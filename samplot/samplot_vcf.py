#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create samplot vcf commands to execute and generate
companion HTML image browser.

Note: additional arguments are passed through to samplot plot
"""
from __future__ import print_function

import argparse
from collections import Counter
import logging
import operator
import os
import random
import sys
import re

import pysam
from jinja2 import Environment, FileSystemLoader, select_autoescape

try:
    from shlex import quote
except ImportError:
    from pipes import quote

from .samplot import add_plot


logger = logging.getLogger(__name__)

cmp_lookup = {
    ">": operator.gt,  # e.g. DHFC < 0.5
    "<": operator.lt,
    "<=": operator.le,
    ">=": operator.ge,
    "==": operator.eq,
    "contains": operator.contains,  # e.g. CSQ contains HIGH
    "exists": lambda a, b: True,  # e.g. exists smoove_gene
}


class Sample(object):
    __slots__ = [
        "family_id",
        "id",
        "paternal_id",
        "maternal_id",
        "mom",
        "dad",
        "kids",
        "i",
    ]

    def __init__(self, line):
        toks = line.rstrip().split()
        self.family_id = toks[0]
        self.id = toks[1]
        self.paternal_id = toks[2]
        self.maternal_id = toks[3]
        self.kids = []
        self.i = -1  # index in the vcf.

    def __repr__(self):
        return "Sample(id:{id},paternal_id:{pid},maternal_id:{mid})".format(
            id=self.id, pid=self.paternal_id, mid=self.maternal_id
        )


def flatten(value, sep=","):
    """
    >>> flatten([1,2,3,4])
    '1,2,3,4'
    >>> flatten((5,6))
    '5,6'
    >>> flatten(0.987654321)
    '0.987654'
    >>> flatten(7)
    '7'
    >>> flatten("flatten")
    'flatten'
    """
    flat = None
    # tuple or list
    if isinstance(value, tuple) or isinstance(value, list):
        flat = sep.join([str(i) for i in value])
    # reformats long float values
    elif isinstance(value, float):
        flat = "%.6f" % (value,)
    # string and int
    else:
        flat = str(value)
    return flat


def get_format_fields(ids, variant):
    """
    args:
        ids (list) - list of FORMAT field IDs, e.g. ['AS', 'AP', 'DHFFC']
        variant (pysam.libcbcf.VariantRecord)

    returns:
        list
    """
    sample_format = []
    for i, sample_fields in enumerate(variant.samples.values()):
        for field_id in ids:
            sample_field_val = flatten(sample_fields.get(field_id, ""))
            if sample_field_val:
                if len(sample_format) < i + 1:
                    sample_format.append("")
                else:
                    sample_format[i] += " "
                sample_format[i] += "{}={}".format(field_id, sample_field_val)
    return sample_format


def get_format_title(samples, ids, variant):
    """
    args:
        samples (list) - list of sample IDs in order of VCF annotations
        ids (list) - list of FORMAT field IDs, e.g. ['AS', 'AP', 'DHFFC']
        variant (pysam.libcbcf.VariantRecord)

    returns:
        dict
    """
    fields = get_format_fields(ids, variant)
    return dict(zip(samples, fields))


def make_plot_titles(samples, attr_values):
    """
    keeping this method separate in the event we add more things to the title

    args:
        samples (list) - list of sample IDs
        attr_values (str) - string of VCF FORMAT values

    returns:
        dict

    >>> make_plot_titles(
        ['s1', 's2', 's3'],
            {
                's1': 'AS=0 AP=0',
                's2': 'AS=0 AP=1',
                's3': 'AS=1 AP=1'
            }
        )
    {
        's1': "'s1 AS=0 AP=0'",
        's2': "'s2 AS=0 AP=1'",
        's3': "'s3 AS=1 AP=1'"
    }
    """
    plot_titles = dict()
    for sample in samples:
        if sample in attr_values:
            plot_titles[sample] = quote("%s %s" % (sample, attr_values[sample]))
    return plot_titles


def get_overlap(
    tabix,
    chrom,
    start,
    end,
    priority=["exon", "gene", "transcript", "cds"],
    no_hit="intergenic",
    fix_chr=True,
):
    """
    args:
        tabix (pysam.libctabix.TabixFile) - open TabixFile
        chrom (str)
        start (int)
        end (int)
        priority (Optional[list]) - order of preferred region annotation
        no_hit (Optional[str]) - use this annotation if no matches among priority
        fix_chr (Optional[bool]) - try to fetch a region using both
                            non-'chr' and 'chr' prefix on failures

    returns:
        str
    """
    overlaps = None
    try:
        overlaps = set(
            [i.split("\t")[2].lower() for i in tabix.fetch(chrom, start, end)]
        )
    except IndexError:
        # probably not a gff3
        logger.warning("Invalid annotation file specified for --gff3")
        overlaps = None
    except ValueError:
        if fix_chr:
            # try removing chr
            if chrom.startswith("chr"):
                overlaps = get_overlap(
                    tabix, chrom[3:], start, end, priority, no_hit, False
                )
            # or adding chr
            else:
                overlaps = get_overlap(
                    tabix,
                    "chr{chrom}".format(chrom=chrom),
                    start,
                    end,
                    priority,
                    no_hit,
                    False,
                )
    except:
        # bad regions
        logger.warning(
            "Error fetching {chrom}:{start}-{end}".format(
                chrom=chrom, start=start, end=end
            )
        )
        overlaps = None

    overlap = ""
    if overlaps:
        for feature in priority:
            if feature in overlaps:
                overlap = feature
                break
    else:
        # fetching overlaps failed
        overlap = "unknown"

    if not overlap and no_hit:
        overlap = no_hit
    return overlap


def parse_ped(path, vcf_samples=None):
    if path is None:
        return {}
    samples = []
    look = {}
    for line in open(path):
        samples.append(Sample(line))
        look[samples[-1].id] = samples[-1]

    for s in samples:
        s.dad = look.get(s.paternal_id)
        if s.dad is not None:
            s.dad.kids.append(s)
        s.mom = look.get(s.maternal_id)
        if s.mom is not None:
            s.mom.kids.append(s)
    # match these samples to the ones in the VCF.
    if vcf_samples is not None:
        result = []
        for i, variant_sample in enumerate(vcf_samples):
            if variant_sample not in look:
                continue
            result.append(next(s for s in samples if s.id == variant_sample))
            result[-1].i = i
        samples = result

    return {s.id: s for s in samples}


def get_names_to_bams(bams, name_list=None):
    """
    get mapping from names (read group samples) to bam paths)
    this is useful because the VCF has the names and we'll want the bam paths
    for those samples
    if name_list is passed in as a parameter those will be used instead
    """
    names = {}
    if name_list:
        if len(name_list) != len(bams):
            logger.error("List of sample IDs does not match list of alignment files.")
            sys.exit(1)

        for i, p in enumerate(bams):
            names[name_list[i]] = p
    else:
        for p in bams:
            b = pysam.AlignmentFile(p)
            # TODO - catch specific exception
            try:
                names[b.header["RG"][0]["SM"]] = p
            except Exception as e:
                logger.error("No RG field in alignment file {}".format(p))
                logger.error("Include ordered list of sample IDs to avoid this error")
                print(e, file=sys.stderr)
                sys.exit(1)
    return names


def tryfloat(v):
    try:
        return float(v)
    except:
        return v


def to_exprs(astr):
    """
    an expr is just a 3-tuple of "name", fn, value"
    e.g. "DHFFC", operator.lt, 0.7"
    >>> to_exprs("DHFFC < 0.5 & SVTYPE == 'DEL'")
    [('DHFFC', <built-in function lt>, 0.5), ('SVTYPE', <built-in function eq>, 'DEL')]

    >>> to_exprs("CSQ contains 'HIGH'")
    [('CSQ', <built-in function contains>, 'HIGH')]
    """
    astr = (x.strip() for x in astr.strip().split("&"))
    result = []
    for a in astr:
        a = [x.strip() for x in a.split()]
        if len(a) == 2:
            assert a[1] == "exists", ("bad expression", a)
            a.append("extra_arg")
        assert len(a) == 3, ("bad expression", a)
        assert a[1] in cmp_lookup, (
            "comparison:"
            + a[1]
            + " not supported. must be one of:"
            + ",".join(cmp_lookup.keys())
        )
        result.append((a[0], cmp_lookup[a[1]], tryfloat(a[2].strip("'").strip('"'))))
    return result


def check_expr(vdict, expr):
    """
    >>> check_expr({"CSQ": "asdfHIGHasdf"},
            to_exprs("CSQ contains 'HIGH'"))
    True

    >>> check_expr({"CSQ": "asdfHIGHasdf", "DHFC": 1.1},
            to_exprs("CSQ contains 'HIGH' & DHFC < 0.5"))
    False

    >>> check_expr({"CSQ": "asdfHIGHasdf", "DHFC": 1.1},
            to_exprs("CSQ contains 'HIGH' & DHFC < 1.5"))
    True

    >>> check_expr({"smoove_gene": "asdf"},
            to_exprs("smoove_gene exists"))
    True

    >>> check_expr({"smooe_gene": "asdf"},
            to_exprs("smoove_gene exists"))
    False

    >>> check_expr({"smoove_gene": ""},
            to_exprs("smoove_gene exists"))
    True
    """

    # a single set of exprs must be "anded"
    for name, fcmp, val in expr:
        # NOTE: asking for a missing annotation will return false.
        if name not in vdict:
            return False
        if not fcmp(vdict[name], val):
            return False
    return True


def make_single(vdict):
    """
    >>> d = {"xx": (1,)}
    >>> make_single(d)
    {'xx': 1}
    """
    for k in vdict.keys():
        if isinstance(vdict[k], tuple) and len(vdict[k]) == 1:
            vdict[k] = vdict[k][0]
    return vdict


def get_dn_row(ped_samples):
    for s in ped_samples.values():
        if s.mom is not None and s.dad is not None:
            return '{title:"de novo", field:"dn"}'
    return ""


def read_important_regions(bedfilename):
    if not bedfilename:
        return None
    important_regions = {}
    with open(bedfilename, "r") as bedfile:
        for line in bedfile:
            pos_fields = line.strip().split()
            region_string = "_".join(pos_fields[1:3])
            if pos_fields[0] not in important_regions:
                important_regions[pos_fields[0]] = []
            important_regions[pos_fields[0]].append(region_string)

    return important_regions


def var_in_important_regions(important_regions, chrom, start, end, svtype):
    if not important_regions:
        # if no important regions are set all locations are valid
        return True

    if chrom in important_regions:
        for region in important_regions[chrom]:
            region_st, region_end = [int(x) for x in region.split("_")]
            if (
                region_st <= start <= region_end
                or region_st <= end <= region_end
                or start <= region_st <= end
            ):
                return True

    logger.debug(
        "Skipping {} at {}:{}-{}, outside important_regions coordinates".format(
            svtype, chrom, start, end
        )
    )
    return False


def cram_input(bams):
    for bam in bams:
        if bam.endswith(".cram"):
            return True
    return False


def above_call_rate(gts, sample_count, min_call_rate, svtype, chrom, start, end):
    """
    skips variants with call rate below min_call_rate if set
    """
    if not min_call_rate:
        return True

    call_rate = (sample_count - sum(None in g for g in gts)) / sample_count
    if min_call_rate and (call_rate < min_call_rate):
        logger.debug(
            (
                "Skipping {} at {}:{}-{}, call rate of variant "
                + "({}) below min_call_rate"
            ).format(svtype, chrom, start, end, call_rate),
        )
        return False
    return True


def below_max_hets(gts, max_hets, svtype, chrom, start, end):
    """
    skips variants with more than max_hets heterozygotes
    if max_hets is set
    """
    if not max_hets:
        return False

    # requisite hets/hom-alts
    het_count = sum(sum(x) >= 1 for x in gts if None not in x)
    if het_count > max_hets:
        logger.debug(
            "Skipping {} at {}:{}-{}, more than max_hets heterozygotes".format(
                svtype, chrom, start, end
            )
        )
        return False
    return True


def no_variant_found(gts, svtype, chrom, start, end):
    """
    skips variants with no non-ref samples
    """
    if not any(sum(x) > 0 for x in gts if None not in x):
        logger.debug(
            "Skipping {} at {}:{}-{}, no samples have non-ref genotypes".format(
                svtype, chrom, start, end
            )
        )
        return True
    return False


def get_plottable_samples(
    gts, variant, plot_all, filters, svtype, chrom, start, end,
):
    """
    gets the samples and indices for all those which need to be plotted,
    which means passing filters and, if not plot_all, having a nonref genotype
    """
    if plot_all:
        test_idxs = [i for i, gt in enumerate(gts)]
        if len(test_idxs) == 0:
            logger.debug(
                "No samples found for {} at {}:{}-{}".format(svtype, chrom, start, end)
            )
    else:
        test_idxs = [i for i, gt in enumerate(gts) if None not in gt and sum(gt) > 0]
        if len(test_idxs) == 0:
            logger.debug(
                "No non-reference samples found for {} at {}:{}-{}".format(
                    svtype, chrom, start, end
                )
            )

    test_samples = [s for i, s in enumerate(variant.samples.values()) if i in test_idxs]

    # apply filters if set
    if len(filters) == 0:
        idxs = test_idxs
    else:
        idxs = []
        odict = make_single(dict(variant.info.items()))
        for i, ts in enumerate(test_samples):
            vdict = odict.copy()
            vdict.update(make_single(dict(ts.items())))

            if any(check_expr(vdict, fs) for fs in filters):
                idxs.append(test_idxs[i])
    if len(idxs) == 0:
        logger.debug(
            "No samples pass filters for {} at {}:{}-{}".format(
                svtype, chrom, start, end
            )
        )
    return idxs, test_samples


def get_variant_samples(
    idxs, vcf_samples, names_to_bams, svtype, chrom, start, end,
):
    """
    gets the samples that need to be plotted and have alignment files assigned
    """
    variant_samples = []
    for i in idxs:
        if vcf_samples[i] in names_to_bams:
            variant_samples.append(vcf_samples[i])
    if len(variant_samples) == 0:
        logger.debug(
            (
                "Skipping {} at {}:{}-{}, no plottable samples "
                + "with matched alignment files"
            ).format(svtype, chrom, start, end),
        )
    return variant_samples


def get_denovos(
    denovo_row,
    test_samples,
    variant_samples,
    ped_samples,
    svtype,
    chrom,
    start,
    end,
    dn_only,
):
    """
    we call it a de novo if the sample passed the filters but the mom and
    dad had homref genotypes before filtering.
    so stringent filtering on the kid and lenient on parents.
    """
    denovo_svs = []
    if denovo_row != "":
        test_sample_names = {s.name for s in test_samples}
        for variant_sample in variant_samples:
            sample = ped_samples[variant_sample]
            if sample.mom is None or sample.dad is None:
                continue
            if (
                sample.mom.id not in test_sample_names
                and sample.dad.id not in test_sample_names
            ):
                denovo_svs.append(sample.id)

    if len(denovo_svs) <= 0 and dn_only:
        logger.debug(
            "Skipping {} at {}:{}-{}, dn_only selected and no de novos found".format(
                svtype, chrom, start, end
            ),
        )
    return denovo_svs


def get_family_controls(
    ped,
    denovo_svs,
    variant_samples,
    ped_samples,
    max_hets,
    bams,
    names_to_bams,
    vcf_samples_set,
):
    """
    tries to find family members to use as controls for putative de novos
    """
    # do DN samples first so we can see parents.
    # TODO also need to do the non-denovos as they seem to have been forgotten
    for variant_sample in denovo_svs + [
        x for x in variant_samples if x not in denovo_svs
    ]:
        sample = ped_samples.get(variant_sample)
        if sample is None:
            continue
        if (
            sample.mom is not None
            and sample.mom.id not in variant_samples
            and sample.mom.id in vcf_samples_set
        ):
            variant_samples.append("mom-of-%s[%s]" % (variant_sample, sample.mom.id))
            bams.append(names_to_bams[sample.mom.id])
        if (
            sample.dad is not None
            and sample.dad.id not in variant_samples
            and sample.dad.id in vcf_samples_set
        ):
            variant_samples.append("dad-of-%s[%s]" % (variant_sample, sample.dad.id))
            bams.append(names_to_bams[sample.dad.id])
        for kid in sample.kids:
            if kid.id not in variant_samples and kid.id in vcf_samples_set:
                variant_samples.append("kid-of-%s[%s]" % (variant_sample, kid.id))
                bams.append(names_to_bams[kid.id])
            if max_hets:
                if len(bams) > 1.5 * max_hets:
                    break
        if max_hets:
            if len(bams) > 1.5 * max_hets:
                break
    return variant_samples, bams


def get_nonfamily_controls(
    gts, vcf_samples, variant_samples, names_to_bams, min_entries, bams
):
    # extend with some controls:
    hom_ref_idxs = [
        i for i, gt in enumerate(gts) if len(gt) == 2 and gt[0] == 0 and gt[1] == 0
    ]

    if len(hom_ref_idxs) > 3:
        random.shuffle(hom_ref_idxs)

    hom_ref_samples = []
    for i in hom_ref_idxs:
        if vcf_samples[i] in names_to_bams:
            hom_ref_samples.append(vcf_samples[i])

    to_add_count = min_entries - len(bams)
    bams.extend(names_to_bams[s] for s in hom_ref_samples[:to_add_count])
    variant_samples += ["control-sample:" + s for s in hom_ref_samples[:to_add_count]]
    return variant_samples, bams


def create_metadata(
    variant,
    translocation_chrom,
    svtype,
    sample_str,
    n_samples,
    annotations,
    denovo_row,
    denovo_svs,
):
    """
    creates a dict with the info about the SV
    that will be used in the website
    """
    data_dict = {
        "chrom": variant.chrom,
        "chrom2": translocation_chrom,
        "start": variant.start,
        "end": variant.stop,
        "svtype": svtype,
        "svlength": variant.stop - variant.start,
        "samples": sample_str,
        "nsamples": n_samples,
    }
    if annotations:
        data_dict["overlaps"] = get_overlap(
            annotations, variant.chrom, variant.start, variant.stop
        )
    if denovo_row != "":
        data_dict["dn"] = ",".join(denovo_svs)
    return data_dict


def format_template(
    variant,
    data_dict,
    max_entries,
    bams,
    variant_samples,
    plot_titles,
    out_dir,
    output_type,
    svtype,
    downsample,
    pass_through_args,
):
    """
    formates the template string for generation of the final command
    """
    if data_dict["chrom2"] is None:
        figname_template = "{svtype}_{chrom}_{start}_{end}.{itype}"
    else:
        figname_template = "{svtype}_{chrom}_{start}_{chrom2}_{end}.{itype}"

    fig_path = os.path.join(
        out_dir, figname_template.format(itype=output_type, **data_dict),
    )

    if "CIPOS" in variant.info:
        v = variant.info["CIPOS"]
        cipos = "--start_ci '%s,%s'" % (abs(int(v[0])), abs(int(v[1])))
    else:
        cipos = ""
    if "CIEND" in variant.info:
        v = variant.info["CIEND"]
        ciend = "--end_ci '%s,%s'" % (abs(int(v[0])), abs(int(v[1])))
    else:
        ciend = ""
    # dynamically set Z to speed drawing and remove noise for larger events
    z = 3
    if variant.stop - variant.start > 2000:
        z = 4
    if variant.stop - variant.start > 10000:
        z = 6
    if variant.stop - variant.start > 20000:
        z = 9
    if data_dict["chrom2"] is None:
        z = 3

    if max_entries:
        bams = bams[:max_entries]
        variant_samples = variant_samples[:max_entries]

    # update titles based on FORMAT fields requested
    title_list = list()
    for variant_sample in variant_samples:
        if variant_sample in plot_titles:
            title_list.append(plot_titles[variant_sample])
        else:
            title_list.append(variant_sample)

    start = variant.start
    stop = variant.stop
    start2 = None
    stop2 = None

    if data_dict["chrom2"] is None:
        template = (
            "samplot plot {extra_args} -z {z} -n {titles} "
            + "{cipos} {ciend} {svtype} -c {chrom} -s {start} "
            + "-e {end} -o {fig_path} -d {downsample} -b {bams}"
        )
    else:
        template = (
            "samplot plot {extra_args} -z {z} -n {titles} "
            + "{cipos} {ciend} {svtype} -c {chrom} -s {start} "
            + "-e {end} -o {fig_path} -d {downsample} -b {bams} "
            + "-c {chrom2} -s {start2} -e {end2}"
        )
        # For interchromosomal variants the 2nd breakpoint position should 
        # not be encoded in INFO/END tag although some callers still do this. 
        # Currently it is unclear if there is a good replacement. Delly uses 
        # INFO/POS2 for this, GATK-SV uses INFO/END2, dysgu uses INFO/CHR2_POS.
        # see:  https://github.com/dellytools/delly/issues/159
        # see: https://gatk.broadinstitute.org/hc/en-us/articles/5334587352219-How-to-interpret-SV-VCFs
        # TODO - if the SV breakpoints are specified in the ALT field one 
        #        could use this info to get the 2nd breakpoint position
        if "POS2" in variant.info:
            start2 = variant.info["POS2"]
        elif "END2" in variant.info:
            start2 = variant.info["END2"]
        elif "CHR2_POS" in variant.info:
            start2 = variant.info["CHR2_POS"]
        else:
            start2 = stop
            # Update stop if INFO/END denotes the 2nd breakpoint
            stop = start + 1

        stop2 = start2 + 1
        

    command = template.format(
        extra_args=" ".join(pass_through_args),
        bams=" ".join(bams),
        titles=" ".join(title_list),
        z=z,
        cipos=cipos,
        ciend=ciend,
        svtype="-t " + svtype if svtype != "SV" else "",
        fig_path=fig_path,
        chrom=variant.chrom,
        start=start,
        end=stop,
        downsample=downsample,
        chrom2=data_dict["chrom2"],
        start2=start2,
        end2=stop2,
    ) + "\n"
    return command


def write_site(table_data, out_dir, output_type, annotations, denovo_row):
    # grab the template
    env = Environment(
        loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), "templates")),
        autoescape=select_autoescape(["html"]),
    )
    html_template = env.get_template("samplot_vcf.html")
    # write index.html
    with open("{out_dir}/index.html".format(out_dir=out_dir), "w") as fh:
        print(
            html_template.render(
                data=table_data,
                plot_type=output_type,
                gff3="true" if annotations else "false",
                denovo="true" if denovo_row else "false",
            ),
            file=fh,
        )


def is_simply_skippable(
    variant,
    vcf_samples,
    gts,
    important_regions,
    max_mb,
    min_bp,
    min_call_rate,
    max_hets,
    plot_all,
    translocation_chrom,
):
    """
    checks several basic terms that could filter this variant out
    specifically, if the variant type is INS,
    or fails the important regions,
    max_mb, min_bp, min_call_rate, or max_hets filters
    """
    svtype = variant.info.get("SVTYPE", "SV")

    # skips variants outside important regions if those are set
    if not var_in_important_regions(
        important_regions, variant.chrom, variant.start, variant.stop, svtype,
    ):
        return True

    # skips insertions
    if svtype in ("INS"):
        logger.debug(
            "Skipping {} at {}:{}-{}, INS type not supported".format(
                svtype, variant.chrom, variant.start, variant.stop
            )
        )
        return True

    # skips variants over max_mb length, if set
    if max_mb and (variant.stop - variant.start > max_mb * 1000000):
        logger.debug(
            "Skipping {} at {}:{}-{}, variant length greater than max_mb".format(
                svtype, variant.chrom, variant.start, variant.stop
            )
        )
        return True
    
    # skips variants under min_bp, if set
    if (variant.stop - variant.start < min_bp) and translocation_chrom is None:
        logger.debug(
            "Skipping {} at {}:{}-{}, variant length shorter than min_bp".format(
                svtype, variant.chrom, variant.start, variant.stop
            )
        )
        return True

    # skips variants if the call rate is below min_call_rate, if set
    if not above_call_rate(
        gts,
        len(vcf_samples),
        min_call_rate,
        svtype,
        variant.chrom,
        variant.start,
        variant.stop,
    ):
        return True

    # skips variants if there are more hets than max_hets, if set
    if below_max_hets(
        gts, max_hets, svtype, variant.chrom, variant.start, variant.stop
    ):
        return True

    # skips variants where no sample is non-ref, if plot_all is not set
    if not plot_all:
        if no_variant_found(
            gts, svtype, variant.chrom, variant.start, variant.stop
        ):
            return True

    return False


def generate_commands(
    vcf,
    plot_all,
    max_mb,
    min_bp,
    min_call_rate,
    max_hets,
    dn_only,
    ped,
    important_regions,
    format_field_ids,
    min_entries,
    max_entries,
    out_dir,
    output_type,
    downsample,
    filters,
    ped_samples,
    denovo_row,
    names_to_bams,
    annotations,
    pass_through_args,
):
    """
    for every variant in vcf, process and output plot
    command - if and only if it passes filters
    """
    commands = []
    table_data = []
    vcf_samples = vcf.header.samples
    vcf_samples_set = set(vcf_samples)
    vcf_samples_list = list(vcf_samples)
    vcf_stats = Counter()

    # Check if VCF samples match BAMs
    if vcf_samples_set != set(names_to_bams):
        missing_vcf_samples = vcf_samples_set - set(names_to_bams)
        missing_bam_samples = set(names_to_bams) - vcf_samples_set
        logger.warning(
            "VCF samples and BAMs do not match. "
            "This may be due to different sample names in the VCF and BAMs."
        )
        if missing_vcf_samples:
            logger.warning(
                "VCF samples missing from BAMs: {}".format(", ".join(missing_vcf_samples))
            )
        if missing_bam_samples:
            logger.warning(
                "BAMs missing from VCF samples: {}".format(", ".join(missing_bam_samples))
            )
    
    var_count = 0 # define it here, in case of zero variants (otherwise it will raise NameError when used below)
    for var_count, variant in enumerate(vcf):
        translocation_chrom = None
        svtype = variant.info.get("SVTYPE", "SV")

        # get genotypes
        gts = [s.get("GT", (None, None)) for s in variant.samples.values()]

        # handle translocations
        if svtype in ["BND", "TRA"]:
            try:
                translocation_chrom = variant.info.get("CHR2")
            except (KeyError, ValueError) as e:
                logger.debug(e)
                logger.info(f"Translocation {svtype} on {variant.chrom}:{variant.start}"
                              "skipped due to missing CHR2 INFO field.")

        if is_simply_skippable(
            variant,
            vcf_samples,
            gts,
            important_regions,
            max_mb,
            min_bp,
            min_call_rate,
            max_hets,
            plot_all,
            translocation_chrom,
        ):
            vcf_stats["Skipped"] += 1
            continue

        # gets the list of samples to plot
        # skips ref samples if plot_all isn't set
        # and applies user-defined filters
        idxs, test_samples = get_plottable_samples(
            gts,
            variant,
            plot_all,
            filters,
            svtype,
            variant.chrom,
            variant.start,
            variant.stop,
        )
        if len(idxs) == 0:
            vcf_stats["No plottable samples"] += 1
            continue

        # matches alignment files to variant samples
        variant_samples = get_variant_samples(
            idxs,
            vcf_samples,
            names_to_bams,
            svtype,
            variant.chrom,
            variant.start,
            variant.stop,
        )
        if len(variant_samples) <= 0:
            vcf_stats["No plottable samples with matched BAM"] += 1
            continue

        bams = [names_to_bams[s] for s in variant_samples]

        # finds putative de novo variants
        denovo_svs = get_denovos(
            denovo_row,
            test_samples,
            variant_samples,
            ped_samples,
            svtype,
            variant.chrom,
            variant.start,
            variant.stop,
            dn_only,
        )
        if dn_only and (len(denovo_svs) <= 0):
            vcf_stats["Non de novo ('--dn_only' specified)"] += 1
            continue

        # save fields for the html.
        n_samples = len(variant_samples)
        # semi-colon delimited eases CSV export from HTML
        sample_str = ";".join(variant_samples)
        # dict holding sample to FORMAT title string
        plot_titles = dict()
        if format_field_ids:
            format_attrs = get_format_title(vcf_samples_list, format_field_ids, variant)
            plot_titles = make_plot_titles(variant_samples, format_attrs)

        # get control samples if possible
        # try to get family members if ped is set
        # and reference samples is ped is not set
        if ped is not None:
            variant_samples, bams = get_family_controls(
                ped,
                denovo_svs,
                variant_samples,
                ped_samples,
                max_hets,
                bams,
                names_to_bams,
                vcf_samples_set,
            )
        elif min_entries and len(bams) < min_entries:
            variant_samples, bams = get_nonfamily_controls(
                gts, vcf_samples, variant_samples, names_to_bams, min_entries, bams
            )

        data_dict = create_metadata(
            variant,
            translocation_chrom,
            svtype,
            sample_str,
            n_samples,
            annotations,
            denovo_row,
            denovo_svs,
        )
        table_data.append(data_dict)

        command = format_template(
            variant,
            data_dict,
            max_entries,
            bams,
            variant_samples,
            plot_titles,
            out_dir,
            output_type,
            svtype,
            downsample,
            pass_through_args,
        )
        commands.append(command)


    logger.debug("VCF entry count: {}".format(var_count + 1))
    if vcf_stats:
        logger.debug("VCF entrys filtered out: {}".format(sum(vcf_stats.values())))
        for reason, count in vcf_stats.items():
            logger.debug(" - {}: {}".format(reason, count))

    return commands, table_data


def run_plot_command(command_string: str):
    # Setup a parser for translating the command_string
    parent_parser = argparse.ArgumentParser()
    sub_parser = parent_parser.add_subparsers(title="[sub-commands]", dest="command")
    add_plot(sub_parser)

    # Convert command_string to list and remove leading 'samplot' argument
    # Taken from https://stackoverflow.com/a/524796.
    # NOTE: If python2 is dropped, `shlex.split` could be used for simpler syntax
    command = [p.strip("'") for p in re.split("( |\\\".*?\\\"|'.*?')", command_string.strip()) if p.strip()]
    command = command[1:]

    # Skipped parse_known_args here since extra_args are not used in `samplot plot`.
    # This means that any fauly extra arguments given to `samplot vcf` will raise
    # and error here
    args = parent_parser.parse_args(command)
    args.func(parent_parser, args)


def vcf(parser, args, pass_through_args):
    """
    Generate commands and html for plotting/reviewing variants from VCF
    """
    if args.debug:
        logger.setLevel(logging.DEBUG)

    if args.dn_only and not args.ped:
        logger.error("Missing --ped, required when using --dn_only")
        sys.exit(1)

    if cram_input(args.bams):
        if "-r" not in pass_through_args and "--reference" not in pass_through_args:
            logger.error(
                "ERROR: missing reference file required for CRAM. "
                + "Use -r option. (Run `samplot.py -h` for more help)"
            )
            sys.exit(1)

    vcf = pysam.VariantFile(args.vcf)
    vcf_samples = vcf.header.samples

    annotations = None
    if args.gff3:
        annotations = pysam.TabixFile(args.gff3)

    filters = [to_exprs(f) for f in args.filter]

    ped_samples = parse_ped(args.ped, vcf_samples)

    # this is empty unless we have a sample with both parents defined.
    denovo_row = get_dn_row(ped_samples)

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # connect the sample IDs to bam files
    names_to_bams = get_names_to_bams(args.bams, args.sample_ids)

    # check that at least one sample is can be plotted  
    if not any(vcf_sample in names_to_bams for vcf_sample in vcf_samples):
        other = "'--sample_ids'" if args.sample_ids else "BAM"
        logger.error("Samples in VCF do not match samples specified in {}".format(other))
        logger.error("VCF samples: {}".format(', '.join(vcf_samples)))
        logger.error("{} samples: {}".format(other, ', '.join(vcf_samples)))
        sys.exit(1)

    # if important regions are included, load those intervals
    # and only show SVs inside them
    important_regions = read_important_regions(args.important_regions)

    # user-requested FORMAT fields to add to plot title
    format_field_ids = None
    if args.format:
        format_field_ids = args.format.split(",")

    # for every variant in vcf, process and output plot
    # command - if and only if it passes filters
    commands, table_data = generate_commands(
        vcf,
        args.plot_all,
        args.max_mb,
        args.min_bp,
        args.min_call_rate,
        args.max_hets,
        args.dn_only,
        args.ped,
        important_regions,
        format_field_ids,
        args.min_entries,
        args.max_entries,
        args.out_dir,
        args.output_type,
        args.downsample,
        filters,
        ped_samples,
        denovo_row,
        names_to_bams,
        annotations,
        pass_through_args,
    )

    write_site(table_data, args.out_dir, args.output_type, annotations, denovo_row)

    if args.manual_run:
        with open(args.command_file, "w") as outfile:
            outfile.writelines(commands)
    else:
        if args.threads == 1:
            for command in commands:
                run_plot_command(command)
        else:
            from multiprocessing import Pool
            with Pool(processes=args.threads) as pool:
                pool.map(run_plot_command, commands)


def add_vcf(parent_parser):
    """Defines allowed arguments for samplot's vcf plotter
    """
    import doctest

    parser = parent_parser.add_parser(
        "vcf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Generates commands to plot images with `samplot plot`,"
        + " using VCF file to define regions",
    )

    if len(sys.argv) > 1 and sys.argv[1] == "test":
        r = doctest.testmod()
        print(r)
        sys.exit(r.failed)

    parser.add_argument("--vcf", "-v", help="VCF file containing structural variants")
    parser.add_argument(
        "-d", "--out-dir", help="path to write output images", default="samplot-out",
    )
    parser.add_argument(
        "--ped", help="path to ped (or .fam) file",
    )
    parser.add_argument(
        "--dn_only",
        help="plots only putative de novo variants (PED file required)",
        action="store_true",
    )
    parser.add_argument(
        "--min_call_rate",
        type=float,
        help="only plot variants with at least this call-rate",
        required=False,
    )
    parser.add_argument(
        "--filter",
        action="append",
        help="simple filter that samples"
        + " must meet. Join multiple filters with '&' "
        + "and specify --filter multiple times for 'or'"
        + " e.g. DHFFC < 0.7 & SVTYPE = 'DEL'",
        default=[],
    )
    parser.add_argument(
        "-O",
        "--output_type",
        choices=("png", "pdf", "eps", "jpg"),
        help="type of output figure",
        default="png",
    )
    parser.add_argument(
        "--max_hets",
        type=int,
        help="only plot variants with at most this many heterozygotes",
        required=False,
    )
    parser.add_argument(
        "--min_entries",
        type=int,
        help="try to include homref samples as controls to get this many samples in plot",
        default=6,
        required=False,
    )
    parser.add_argument(
        "--max_entries",
        type=int,
        help="only plot at most this many heterozygotes",
        default=10,
        required=False,
    )
    parser.add_argument(
        "--max_mb",
        type=int,
        help="skip variants longer than this many megabases",
        required=False,
    )
    parser.add_argument(
        "--min_bp",
        type=int,
        help="skip variants shorter than this many bases",
        default=20,
    )
    parser.add_argument(
        "--important_regions",
        help="only report variants that overlap regions in this bed file",
        required=False,
    )
    parser.add_argument(
        "-b",
        "--bams",
        type=str,
        nargs="+",
        help="Space-delimited list of BAM/CRAM file names",
        required=True,
    )
    parser.add_argument(
        "--sample_ids",
        type=str,
        nargs="+",
        help="Space-delimited list of sample IDs, "
        + "must have same order as BAM/CRAM file names. "
        + "BAM RG tag required if this is omitted.",
        required=False,
    )
    parser.add_argument(
        "--command_file",
        help="store commands in this file.",
        default="samplot_vcf_cmds.tmp",
        required=False,
    )
    parser.add_argument(
        "--format",
        default="AS,AP,DHFFC",
        help="comma separated list of FORMAT fields to include in sample plot title",
        required=False,
    )
    parser.add_argument(
        "--gff3",
        help="genomic regions (.gff with .tbi in same directory) "
        + "used when building HTML table and table filters",
        required=False,
    )
    parser.add_argument(
        "--downsample", help="Number of normal reads/pairs to plot", default=1, type=int
    )
    parser.add_argument(
        "--manual_run",
        help="disables auto-run for the plotting commands",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--plot_all",
        help="plots all samples and all variants - "
        + "limited by any filtering arguments set",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of threads to use to generate plots. Default: %(default)s",
    )
    parser.add_argument(
        "--debug",
        help="prints out the reason for skipping any skipped variant entry",
        default=False,
        action="store_true",
    )

    parser.set_defaults(func=vcf)


if __name__ == "__main__":
    print("Run as samplot module with `samplot vcf`")
