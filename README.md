<center><img src="/doc/imgs/samplot_logo_v5.png" width="300"/></center>

<center><img src="/doc/imgs/montage.jpg" width="75%"/></center>

`samplot` is a command line tool for rapid, multi-sample structural variant
visualization. `samplot` takes SV coordinates and bam files and produces
high-quality images that highlight any alignment and depth signals that
substantiate the SV.

## Usage
```
Usage: samplot.py [options]

Options:
  -h, --help            show this help message and exit
  -n TITLES             Plot title (CSV)
  -r REFERENCE          Reference file for CRAM
  -z Z                  Number of stdevs from the mean (default 4)
  -b BAMS               Bam file names (CSV)
  -o OUTPUT_FILE        Output file name
  -s START              Start range
  -e END                End range
  -c CHROM              Chromosome range
  -w WINDOW             Window size (count of bases to include), default(0.5 *
                        len)
  -d MAX_DEPTH          Max number of normal pairs to plot
  -t SV_TYPE            SV type
  -T TRANSCRIPT_FILE    GFF of transcripts
  -A ANNOTATION_FILE    bed.gz tabixed file of transcripts
  -a                    Print commandline arguments
  -H PLOT_HEIGHT        Plot height
  -W PLOT_WIDTH         Plot width
  --common_insert_size  Set common insert size for all plots
```

## Dependencies

* numpy
* matplotlib
* pysam
* statistics

All of these are available from [pip](https://pypi.python.org/pypi/pip).

## Examples: 

Samplot requires either BAM files or CRAM files as primary input. If you use
CRAM, you'll also need a reference genome like one used the the 1000 Genomes
Project
(ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz).

### Basic use case
Using data from NA12878, NA12889, and NA12890 in the 
[1000 Genomes Project](http://www.internationalgenome.org/about), we will
inspect a possible deletion in NA12878 at 4:115928726-115931880 with respect
to that same region in two unrelated samples NA12889 and NA12890.

The following command will create an image of that region:
```
time samplot/src/samplot.py \
    -n NA12878,NA12889,NA12890 \
    -b samplot/test/data/NA12878_restricted.bam,samplot/test/data/NA12889_restricted.bam,samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL

real    0m9.450s
user    0m9.199s
sys     0m0.217s
```

The arguments used above are:

`-n` The names to be shown for each sample in the plot

`-b` The BAM/CRAM files of the samples (comma-separated)

`-o` The name of the output file containing the plot

`-c` The chromosome of the region of interest

`-s` The start location of the region of interest

`-e` The end location of the region of interest

`-t` The type of the variant of interest

This will create an image file named `4_115928726_115931880.png`, shown below:

<img src="/doc/imgs/4_115928726_115931880.png">

### Downsampling "normal" pairs

The runtime of `samplot` can be reduced by only plotting a portion of the concordant 
pair-end reads (+/- strand orientation, within z s.d. of the mean insert size where z 
is a command line option the defaults to 4). If we rerun the prior example, but only plot
a random sampling of 100 normal pairs we get a similar result 3.6X faster.

```
time samplot/src/samplot.py \
    -n NA12878,NA12889,NA12890 \
    -b samplot/test/data/NA12878_restricted.bam,samplot/test/data/NA12889_restricted.bam,samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.d100.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL \
    -d 100

real    0m2.621s
user    0m2.466s
sys     0m0.124s
```

<img src="/doc/imgs/4_115928726_115931880.d100.png">


### Gene and other genomic feature annotations

Gene annotations (tabixed, gff3 file) and genome features (tabixed, bgzipped, bed file) can be 
included in the plots.

Get the gene annotations:
```
wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.gff3.gz
bedtools sort -i Homo_sapiens.GRCh37.82.gff3.gz \
| bgzip -c > Homo_sapiens.GRCh37.82.sort.gff3.gz
tabix Homo_sapiens.GRCh37.82.sort.gff3.gz
```

Get genome annotations, in this case Repeat Masker tracks and a mappability track:
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityUniqueness35bp.bigWig
bigWigToBedGraph wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bed
bgzip wgEncodeDukeMapabilityUniqueness35bp.bed
tabix wgEncodeDukeMapabilityUniqueness35bp.bed.gz

curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
| bgzip -d -c \
| cut -f 6,7,8,13 \
| bedtools sort -i stdin \
| bgzip -c > rmsk.bed.gz
tabix rmsk.bed.gz
```

Plot:
```
samplot/src/samplot.py \
    -n NA12878,NA12889,NA12890 \
    -b samplot/test/data/NA12878_restricted.bam,samplot/test/data/NA12889_restricted.bam,samplot/test/data/NA12890_restricted.bam \
    -o 4_115928726_115931880.d100.genes_reps_map.png \
    -c chr4 \
    -s 115928726 \
    -e 115931880 \
    -t DEL \
    -d 100 \
    -T Homo_sapiens.GRCh37.82.sort.gff3.gz \
    -A rmsk.bed.gz,wgEncodeDukeMapabilityUniqueness35bp.bed.gz

real    0m2.784s
user    0m2.633s
sys 0m0.129s
```

<img src="/doc/imgs/4_115928726_115931880.d100.genes_reps_map.png">

### Generating images from a VCF file
To plot images from all structural variants in a VCF file, use samplot's
`samplot_vcf.sh` script. This accepts a VCF file and the BAM files of samples
you wish to plot, outputting images and related metadata to a directory of your
choosing.

This script is especially useful as part of the 
[SV-plaudit pipeline](https://github.com/jbelyeu/SV-plaudit) and creates
metadata files for
all images which SV-plaudit requires.

```
samplot/src/samplot_vcf.sh \
    -o output_dir \
    -B $HOME/bin/bcftools \
    -S samplot/src/samplot.py \
    -v samplot/test/data/NA12878.trio.svt.subset.vcf \
    samplot/test/data/NA12878_restricted.bam \
    samplot/test/data/NA12889_restricted.bam \
    samplot/test/data/NA12890_restricted.bam
```
The arguments used above are:

`-o` output directory (make this directory before executing)

`-B` Executable file of [bcftools](https://samtools.github.io/bcftools/)

`-S` samplot.py script

`-v` VCF file with variants to plot


#### CRAM inputs
Samplot also support CRAM input, which requires a reference fasta file for
reading as noted above. Notice that the reference file is not included in this
repository due to size. This time we'll plot an interesting duplication at
X:101055330-101067156.

```
samplot/src/samplot.py \
    -n NA12878,NA12889,NA12890 \
    -b samplot/test/data/NA12878_restricted.cram,samplot/test/data/NA12889_restricted.cram,samplot/test/data/NA12890_restricted.cram \
    -o cramX_101055330_101067156.png 
    -c chrX \
    -s 101055330 \
    -e 101067156 \
    -t DUP \
    -r hg19.fa
```

The arguments used above are the same as those used for the basic use case, with the addition of the following:

`-r` The reference file used for reading CRAM files

#### Plotting without the SV 
Samplot can also plot genomic regions that are unrelated to an SV. If you do
not pass the SV type option (`-t`) then the top SV bar will go away and only
the region that is given by `-c` `-s` and `-e` will be displayed.

#### Long read (Oxford nanopore and PacBio) and linked read support
Any alignment that is longer than 1000 bp are treated as a longread, and
the plot design will focus on aligned regions and gaps. Aligned regions are in orange, and gaps follow the same DEL/DUP/INV color code used for short reads. The height of the alignment is based on the size of its largest gap.

<img src="/doc/imgs/longread_del.png">

If the bam file has an MI tag, then the reads will be treated as linked reads.
The plots will be similar to short read plots, but all alignments with the same MI is plotted at the same height according to alignment with the largest gap in the group. A green line connects all alignments in a group.

<img src="/doc/imgs/linkedread_del.png">
