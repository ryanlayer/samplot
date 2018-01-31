# samplot

Creates image views of genomic intervals from read alignment files, optimized for structural variant viewing.

Samplot workflow is simple. Clone this repository and ensure that you have Python 2.7.X (samplot does not currently support Python 3.XX). Also make sure you have the following python libraries:

* numpy
* matplotlib
* pylab
* pysam
* statistics

All of these are available from [pip](https://pypi.python.org/pypi/pip).

Samplot requires either BAM files or CRAM files as primary input (if you use CRAM, you'll also need a reference genome like one used the the 1000 Genomes Project (ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz). Follow the usage examples below to format your commands.

## Usage Examples: 


### Basic use case
We're  using data from NA12878, NA12889, and NA12890 in the [1000 Genomes Project](http://www.internationalgenome.org/about). 

Let's say we have BAM files and want to see what the deletion in NA12878 at 4:115928726-115931880 looks like compared to the parents (NA12889, NA12890). 
The following command will create an image of that region:
```
python src/samplot.py -n NA12878,NA12889,NA12890 -b test/data/alignments/NA12878_restricted.bam,Samplot/test/data/alignments/NA12889_restricted.bam,Samplot/test/data/alignments/NA12890_restricted.bam -o 4_115928726_115931880.png -s 115928726 -e 115931880 -c chr4 -a -t DEL > 4_115928726_115931880.args
```

<img src="/doc/imgs/4_115928726_115931880.png">

### CRAM inputs
Samplot also support CRAM input, which requires a reference fasta file for reading as noted above. Notice that the reference file is not included in this repository due to size.

```
python src/samplot.py -n NA12878,NA12889,NA12890 -b test/data/alignments/NA12878_restricted.cram,Samplot/test/data/alignments/NA12889_restricted.cram,Samplot/test/data/alignments/NA12890_restricted.cram -o cramX_101055330_101067156.png -s 101055330 -e 101067156 -c chrX -a -t DUP -r ~/Research/data/reference/hg19/hg19.fa > cram_X_101055330_101067156.args
```
<img src="doc/imgs/cramX_101055330_101067156.png">
