# samplot

Creates image views of genomic intervals from BAM files.

##Example:
```
python samplot/src/samplot.py -c 2 -s 89161083 -e 89185670 \
    -b "samplot/test/data/low_coverage/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.restricted_sv_regions.20121211.bam" \
	"samplot/test/data/low_coverage/NA12889.mapped.ILLUMINA.bwa.CEU.low_coverage.restricted_sv_regions.20130415.bam" \
	"samplot/test/data/low_coverage/NA12890.mapped.ILLUMINA.bwa.CEU.low_coverage.restricted_sv_regions.20130415.bam" -o "chr2.jpg"
```
