from __future__ import print_function
import pysam
import random
import gzip
import os

HERE = os.path.dirname(__file__)

def xopen(path): return gzip.open(path) if path.endswith(".gz") else open(path)

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

def main(args, pass_through_args):
    vcf = pysam.VariantFile(args.vcf)
    vcf_samples = vcf.header.samples

    if not os.path.exists(args.out_dir):
      os.makedirs(args.out_dir)

    names_to_bams = get_names_to_bams(args.bams)

    for v in vcf:
        gts = [s.get("GT", (None, None)) for s in v.samples.values()]

        if sum(None in g for g in gts) >= args.min_call_rate * len(vcf_samples): continue
        # requisite hets
        if sum(sum(x) == 1 for x in gts if not None in x) > args.max_hets: continue
        if not any(sum(x) > 0 for x in gts if not None in x): continue

        idxs = [i for i, gt in enumerate(gts) if not None in gt and sum(gt) > 0]

        vsamples = [vcf_samples[i] for i in idxs]
        bams = [names_to_bams[s] for s in vsamples]

        if len(bams) < 6:
            # extend with some controls:
            hom_ref_idxs = [i for i, gt in enumerate(gts) if len(gt) == 2 and gt[0] == 0 and gt[1] == 0]
            if len(hom_ref_idxs) > 3:
              random.shuffle(hom_ref_idxs)
              hom_ref_idxs = hom_ref_idxs[:3]

            hsamples = [vcf_samples[i] for i in hom_ref_idxs]
            bams.extend(names_to_bams[s] for s in hsamples)
            vsamples += ["control-sample: " + s for s in hsamples]

        svtype = v.info.get("SVTYPE", "SV")
        if svtype == "BND": continue

        fig_path = "{out_dir}/{svtype}_{chrom}_{start}_{end}.{itype}".format(
                svtype=svtype,
                out_dir=args.out_dir, chrom=v.chrom, start=v.start, end=v.stop, itype=args.output_type)

        print("python {here}/samplot.py {extra_args} --minq 0 -n {titles} {svtype} -c {chrom} -s {start} -e {end} -o {fig_path} -d 1 -w 2000 -b {bams}".format(here=HERE,
            extra_args=" ".join(pass_through_args), bams=" ".join(bams),
            titles=" ".join(vsamples),
            svtype="-t " + svtype if svtype != "SV" else "",
            fig_path=fig_path,
            chrom=v.chrom, start=v.start, end=v.stop))

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser("note that additional arguments are passed through to samplot.py")
    p.add_argument("--vcf", "-v", help="VCF file containing structural variants")
    p.add_argument("-d", "--out-dir", help="path to write output PNGs", default="samplot-out")
    p.add_argument("--min-call-rate", type=float, help="only plot variants with at least this call-rate", default=0.9)
    p.add_argument("-O", "--output-type", choices=("png", "pdf", "eps", "jpg"), help="type of output figure", default="png")
    p.add_argument("--max-hets", type=int, help="only plot variants with at most this many heterozygotes", default=10)
    p.add_argument("-b", "--bams", type=str, nargs="+", help="Space-delimited list of BAM/CRAM file names", required=True)

    a, pass_through_args = p.parse_known_args()
    main(a, pass_through_args)
