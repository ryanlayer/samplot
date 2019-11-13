set -e

#download hg19 reference for cram
FILE="hg19.fa.gz"
if [ ! -f $FILE ]; then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip hg19.fa.gz
    bgzip hg19.fa
fi

#images of each type with all technologies
mkdir -p test_imgs
samplot plot -n Illumina PacBio ONT 10X -t DEL -c 1 -s 24804397 -e 24807302 -o test_imgs/DEL_1_24804397_24807302.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam  -r hg19.fa.gz
samplot plot -n Illumina PacBio ONT 10X -t DUP -c 4 -s 99813786 -e 99817098 -o test_imgs/DUP_4_99813786_99817098.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam -r hg19.fa.gz
samplot plot -n Illumina PacBio ONT 10X -t DUP -c 11 -s 67974431 -e 67975639 -o test_imgs/DUP_11_67974431_67975639.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam -r hg19.fa.gz
samplot plot -n Illumina PacBio ONT 10X -t INV -c 12 -s 12544867 -e 12546613 -o test_imgs/INV_12_12544867_12546613.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam -r hg19.fa.gz

#zoom example
samplot plot -n Illumina PacBio ONT 10X -t DUP -c 4 -s 99813786 -e 99817098 -o test_imgs/DUP_4_99813786_99817098_zoom.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam -r hg19.fa.gz --zoom 1000

#trios with no variant
samplot plot -n HG002 HG003 HG004 -c 1 -s 43059290 -e 43059950 -o test_imgs/1_43059290_43059950.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam

#trios of each type
samplot plot -n HG002 HG003 HG004 -t DEL -c 1 -s 24804397 -e 24807302 -o test_imgs/trio_DEL_1_24804397_24807302.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam
samplot plot -n HG002 HG003 HG004 -t DUP -c 4 -s 99813786 -e 99817098 -o test_imgs/trio_DUP_4_99813786_99817098.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam
samplot plot -n HG002 HG003 HG004 -t DUP -c 11 -s 67974431 -e 67975639 -o test_imgs/trio_DUP_11_67974431_67975639.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam
samplot plot -n HG002 HG003 HG004 -t INV -c 12 -s 12544867 -e 12546613 -o test_imgs/trio_INV_12_12544867_12546613.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam

#create a temporary example website
mkdir -p test_site
samplot vcf -d test_site/ --vcf test.vcf --sample_ids HG002 HG003 HG004 -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam > test_site_cmds.sh
