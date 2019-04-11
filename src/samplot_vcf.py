from __future__ import print_function
import json
import pysam
import operator
import random
import gzip
import os

HERE = os.path.dirname(__file__)
def xopen(path): return gzip.open(path) if path.endswith(".gz") else open(path)

html_tmpl = """<!DOCTYPE html>
<html lang="en">
<head>
<link href="https://unpkg.com/tabulator-tables@4.1.3/dist/css/tabulator_modern.min.css" rel="stylesheet">
<script type="text/javascript" src="https://unpkg.com/tabulator-tables@4.1.3/dist/js/tabulator.min.js"></script>
<script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
</head>
<body>
Click a row to see the <a href="https://github.com/ryanlayer/samplot" target="_blank">samplot</a> image.
<div id="xtable"></div>
<div id="ximg"><img src="https://raw.githubusercontent.com/ryanlayer/samplot/master/doc/imgs/samplot_logo_v5.png"/></div>
<script>

var rowClick = %s

var table = new Tabulator("#xtable", {
    height: 400,
    selectable: 1,
    layout:"fitColumns",
    pagination:"local",
    paginationSize:1000,
    columns: [
      {title:"chrom", field:"chrom", width: 100},
      {title:"start", field:"start"},
      {title:"end", field:"end"},
      {title:"size", field:"size"},
      {title:"type", field:"SVTYPE", width: 80},
      {title:"# of samples", field:"n_samples"},
      {title:"samples", field:"samples", width:200},
      %s
    ],
    rowClick: rowClick,
})
table.setData(%s)

</script>
</body>
</html>
"""


class Sample(object):
    __slots__ = 'family_id id paternal_id maternal_id mom dad kids i'.split()
    def __init__(self, line):
        toks = line.rstrip().split("\t")
        self.family_id = toks[0]
        self.id = toks[1]
        self.paternal_id = toks[2]
        self.maternal_id = toks[3]
        self.kids = []
        self.i = -1 # index in the vcf.

    def __repr__(self):
        return "Sample(id:{id})".format(id=self.id)

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
            if not variant_sample in look: continue
            result.append(next(s for s in samples if s.id == variant_sample))
            result[-1].i = i
        samples = result

    return {s.id: s for s in samples}


def get_names_to_bams(bams):
    """
    get mapping from names (read group samples) to bam paths)
    this is useful because the VCF has the names and we'll want the bam paths
    for those samples
    """
    names = {}
    for p in bams:
        b = pysam.AlignmentFile(p)
        names[b.header["RG"][0]['SM']] = p
    return names


cmp_lookup = {
        '>': operator.gt, # e.g. DHFC < 0.5
        '<': operator.lt,
        '<=': operator.le,
        '>=': operator.ge,
        '==': operator.eq,
        'contains': operator.contains, # e.g. CSQ contains HIGH
        'exists': lambda a,b: True, # e.g. exists smoove_gene
        }

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
        assert a[1] in cmp_lookup, ("comparison:" + a[1] + " not supported. must be one of:" + ",".join(cmp_lookup.keys()))
        result.append((a[0], cmp_lookup[a[1]], tryfloat(a[2].strip("'").strip('"'))))
    return result

def check_expr(vdict, expr):
    """
    >>> check_expr({"CSQ": "asdfHIGHasdf"}, to_exprs("CSQ contains 'HIGH'"))
    True

    >>> check_expr({"CSQ": "asdfHIGHasdf", "DHFC": 1.1}, to_exprs("CSQ contains 'HIGH' & DHFC < 0.5"))
    False

    >>> check_expr({"CSQ": "asdfHIGHasdf", "DHFC": 1.1}, to_exprs("CSQ contains 'HIGH' & DHFC < 1.5"))
    True

    >>> check_expr({"smoove_gene": "asdf"}, to_exprs("smoove_gene exists"))
    True

    >>> check_expr({"smooe_gene": "asdf"}, to_exprs("smoove_gene exists"))
    False

    >>> check_expr({"smoove_gene": ""}, to_exprs("smoove_gene exists"))
    True
    """

    # a single set of exprs must be "anded"
    for name, fcmp, val in expr:
        # NOTE: asking for a missing annotation will return false.
        if not name in vdict: return False
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
    return ''


def main(args, pass_through_args):
    vcf = pysam.VariantFile(args.vcf)
    vcf_samples = vcf.header.samples
    vcf_samples_set = set(vcf_samples)

    filters = [to_exprs(f) for f in args.filter]

    ped_samples = parse_ped(args.ped, vcf_samples)

    # this is empty unless we have a sample with both parents defined.
    dn_row = get_dn_row(ped_samples)
    if not os.path.exists(args.out_dir):
      os.makedirs(args.out_dir)

    names_to_bams = get_names_to_bams(args.bams)
    tabledata = []

    for variant in vcf:
        svtype = variant.info.get("SVTYPE", "SV")
        if svtype in ("BND", "INS"): continue
        if variant.stop - variant.start > args.max_mb * 1000000: continue

        gts = [s.get("GT", (None, None)) for s in variant.samples.values()]

        if sum(None in g for g in gts) >= args.min_call_rate * len(vcf_samples): continue
        # requisite hets/hom-alts
        if sum(sum(x) >= 1 for x in gts if not None in x) > args.max_hets: continue
        if not any(sum(x) > 0 for x in gts if not None in x): continue


        test_idxs = [i for i, gt in enumerate(gts) if not None in gt and sum(gt) > 0]
        test_samples = [s for i, s in enumerate(variant.samples.values()) if i in test_idxs]

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
        if len(idxs) == 0: continue
        is_dn = []

        # we call it a de novo if the sample passed the filters but the mom and
        # dad were not even with the genotype before filtering. so stringent
        # filtering on the kid and lenient on parents.
        variant_samples = []
        for i in idxs:
            if vcf_samples[i] in names_to_bams:
                variant_samples.append(vcf_samples[i])
        if len(variant_samples) <= 0: continue

        bams = [names_to_bams[s] for s in variant_samples]
        if dn_row != "":
            test_sample_names = {s.name for s in test_samples}
            for variant_sample in variant_samples:
                sample = ped_samples[variant_sample]
                if sample.mom is None or sample.dad is None: continue
                if not sample.mom.id in test_sample_names and not sample.dad.id in test_sample_names:
                    is_dn.append(sample.id)


        # save these for the html.
        n_samples = len(variant_samples)
        sample_str = ",".join(variant_samples)

        # try to get family members
        if args.ped is not None:
            # do DN samples first so we can see parents.
            for variant_sample in is_dn + [x for x in variant_samples if not x in is_dn]:
                s = ped_samples.get(variant_sample)
                if s is None: continue
                if s.mom is not None and not s.mom.id in variant_samples and s.mom.id in vcf_samples_set:
                    variant_samples.append("mom-of-%s[%s]" % (variant_sample, s.mom.id))
                    bams.append(names_to_bams[s.mom.id])
                if s.dad is not None and not s.dad.id in variant_samples and s.dad.id in vcf_samples_set:
                    variant_samples.append("dad-of-%s[%s]" % (variant_sample, s.dad.id))
                    bams.append(names_to_bams[s.dad.id])
                for kid in s.kids:
                    if not kid.id in variant_samples and kid.id in vcf_samples_set:
                        variant_samples.append("kid-of-%s[%s]" % (variant_sample, kid.id))
                        bams.append(names_to_bams[kid.id])

                    if len(bams) > 1.5 * args.max_hets: break
                if len(bams) > 1.5 * args.max_hets: break

        elif len(bams) < 6:
            # extend with some controls:
            hom_ref_idxs = [i for i, gt in enumerate(gts) if len(gt) == 2 and gt[0] == 0 and gt[1] == 0]
            if len(hom_ref_idxs) > 3:
              random.shuffle(hom_ref_idxs)
              hom_ref_idxs = hom_ref_idxs[:3]

            hom_ref_samples = []
            for i in hom_ref_idxs:
                if vcf_samples[i] in names_to_bams:
                    hom_ref_samples.append(vcf_samples[i])

            bams.extend(names_to_bams[s] for s in hom_ref_samples)
            variant_samples += ["control-sample:" + s for s in hom_ref_samples]

        data_dict = {"chrom": variant.chrom, "start": variant.start, "end": variant.stop,
                "SVTYPE": svtype, "size": variant.stop - variant.start, "samples":
                sample_str, "n_samples": n_samples}
        if dn_row != "":
            data_dict["dn"] = ",".join(is_dn)
        fig_path = os.path.join(args.out_dir, "{SVTYPE}_{chrom}_{start}_{end}.{itype}".format(
                itype=args.output_type,
                **data_dict))
        tabledata.append(data_dict)

        if "CIPOS" in variant.info:
            v = variant.info["CIPOS"]
            cipos = "--start_ci '%s,%s'" % (abs(v[0]), abs(v[1]))
        else:
            cipos = ""
        if "CIEND" in variant.info:
            v = variant.info["CIEND"]
            ciend = "--end_ci '%s,%s'" % (abs(v[0]), abs(v[1]))
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

        print("python {here}/samplot.py {extra_args} -z {z} --minq 0 -n {titles} {cipos} {ciend} {svtype} -c {chrom} -s {start} -e {end} -o {fig_path} -d 1 -b {bams}".format(here=HERE,
            extra_args=" ".join(pass_through_args), bams=" ".join(bams),
            titles=" ".join(variant_samples),
            z=z,
            cipos=cipos, ciend=ciend,
            svtype="-t " + svtype if svtype != "SV" else "",
            fig_path=fig_path,
            chrom=variant.chrom, start=variant.start, end=variant.stop))


    rowFn = """
    function(e, row) {{
        var d = row.getData()
        var img = d['SVTYPE'] + "_" + d['chrom'] + "_" + d['start'] + "_" + d["end"] + ".{itype}"
        var i = jQuery("#ximg")
        i.empty()
        i.append('<img width="100%" src="' + img + '"/>"')
        return true

    }}""".format(out_dir=args.out_dir, itype=args.output_type)

    with open("{out_dir}/index.html".format(out_dir=args.out_dir), "w") as fh:
        fh.write(html_tmpl % (rowFn, dn_row, json.dumps(tabledata)))

if __name__ == "__main__":
    import argparse
    import sys
    import doctest
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        r = doctest.testmod()
        print(r)
        sys.exit(r.failed)

    p = argparse.ArgumentParser("note that additional arguments are passed through to samplot.py")
    p.add_argument("--vcf", "-v", help="VCF file containing structural variants")
    p.add_argument("-d", "--out-dir", help="path to write output PNGs", default="samplot-out")
    p.add_argument("--ped", help="path ped (or .fam) file")
    p.add_argument("--min-call-rate", type=float, help="only plot variants with at least this call-rate", default=0.95)
    p.add_argument("--filter", action="append", help="simple filter that samples" +
            " must meet. Join multiple filters with '&' and specify --filter multiple times for 'or'" +
            " e.g. DHFFC < 0.7 & SVTYPE = 'DEL'" , default=[])
    p.add_argument("-O", "--output-type", choices=("png", "pdf", "eps", "jpg"), help="type of output figure", default="png")
    p.add_argument("--max-hets", type=int, help="only plot variants with at most this many heterozygotes", default=10)
    p.add_argument("--max-mb", type=int, help="skip variants longer than this many megabases", default=1)
    p.add_argument("-b", "--bams", type=str, nargs="+", help="Space-delimited list of BAM/CRAM file names", required=True)

    a, pass_through_args = p.parse_known_args()
    main(a, pass_through_args)
