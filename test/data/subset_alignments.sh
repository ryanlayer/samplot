#download example regions from GIAB 300X Illumina Ashkenazi Trio
samtools view -h -b -L examples_padded.bed ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam > HG002_Illumina.bam
samtools view -h -b -L examples_padded.bed ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam > HG003_Illumina.bam
samtools view -h -b -L examples_padded.bed ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam > HG004_Illumina.bam

#download example regions from GIAB 300X Illumina Ashkenazi Trio Son
samtools view -h -b -L examples_padded.bed ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/HG002_PB_70x_RG_HP10XtrioRTG.bam > HG002_PacBio.bam
samtools view -h -b -L examples_padded.bed ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh37/NA24385_300G/HG002_10x_84x_RG_HP10xtrioRTG.bam > HG002_10X.bam
samtools view -h -b -C -L examples_padded.bed ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/combined_2018-08-10/HG002_ONTrel2_16x_RG_HP10xtrioRTG.cram > HG002_ONT.bam

#index new alignment files
samtools index HG002_10X.bam
samtools index HG002_Illumina.bam
samtools index HG002_ONT.bam
samtools index HG002_PacBio.bam
samtools index HG003_Illumina.bam
samtools index HG004_Illumina.bam
