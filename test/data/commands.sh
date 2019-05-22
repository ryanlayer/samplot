set -e

#images of each type with all technologies
mkdir -p test_imgs
python ../../src/samplot.py -n Illumina PacBio ONT 10X -t DEL -c 1 -s 24804397 -e 24807302 -o test_imgs/DEL_1_24804397_24807302.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam  -r hg19.fa.gz
python ../../src/samplot.py -n Illumina PacBio ONT 10X -t DUP -c 4 -s 99813786 -e 99817098 -o test_imgs/DUP_4_99813786_99817098.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam -r hg19.fa.gz
python ../../src/samplot.py -n Illumina PacBio ONT 10X -t DUP -c 11 -s 67974431 -e 67975639 -o test_imgs/DUP_11_67974431_67975639.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam -r hg19.fa.gz
python ../../src/samplot.py -n Illumina PacBio ONT 10X -t INV -c 12 -s 12544867 -e 12546613 -o test_imgs/INV_12_12544867_12546613.png -b HG002_Illumina.bam HG002_PacBio.bam HG002_ONT.cram HG002_10X.bam -r hg19.fa.gz

#trios with no variant
python ../../src/samplot.py -n HG002 HG003 HG004 -c 1 -s 43059290 -e 43059950 -o test_imgs/1_43059290_43059950.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam

#trios of each type
python ../../src/samplot.py -n HG002 HG003 HG004 -t DEL -c 1 -s 24804397 -e 24807302 -o test_imgs/trio_DEL_1_24804397_24807302.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam
python ../../src/samplot.py -n HG002 HG003 HG004 -t DUP -c 4 -s 99813786 -e 99817098 -o test_imgs/trio_DUP_4_99813786_99817098.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam
python ../../src/samplot.py -n HG002 HG003 HG004 -t DUP -c 11 -s 67974431 -e 67975639 -o test_imgs/trio_DUP_11_67974431_67975639.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam
python ../../src/samplot.py -n HG002 HG003 HG004 -t INV -c 12 -s 12544867 -e 12546613 -o test_imgs/trio_INV_12_12544867_12546613.png -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam

#create a temporary example website
mkdir -p test_site
python ../../src/samplot_vcf.py -d test_site/ --vcf test.vcf --sample_ids HG002 HG003 HG004 -b HG002_Illumina.bam HG003_Illumina.bam HG004_Illumina.bam > test_site_cmds.sh
bash test_site_cmds.sh
