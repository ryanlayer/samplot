#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

STOP_ON_FAIL=0

bam_1=../data/NA12878_restricted.bam
bam_2=../data/NA12889_restricted.bam
bam_3=../data/NA12890_restricted.bam

sv_chrm=chr4
sv_start=115928730
sv_end=115931875
sv_type=DEL
out_file_name="test.png"

rm -f $out_file_name
run basic_operation \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1,$bam_2,$bam_3 \
        -o $out_file_name \
        -t $sv_type
assert_exit_code 0
assert_equal $out_file_name $( ls $out_file_name )

sample_out_file_name="sample.png"
run sampling_normal \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1,$bam_2,$bam_3 \
        -o $sample_out_file_name \
        -t $sv_type \
        -d 10
assert_exit_code 0
assert_equal $sample_out_file_name $( ls $sample_out_file_name )

sv_chrm=chrX
sv_start=101055330
sv_end=101067156
sv_type=DUP
out_file_name="dup.png"
rm -f $out_file_name

run common_insert_size_scale \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1,$bam_2,$bam_3 \
        -o $out_file_name \
        -t $sv_type \
        -d 10 \
        --common_insert_size
assert_exit_code 0
assert_equal $out_file_name $( ls $out_file_name )

rm -rf img/
mkdir img
vcf_file=../data/NA12878.trio.svt.subset.vcf
run from_vcf \
    ../../src/samplot_vcf.sh \
    -d 10 \
    -o img \
    -v $vcf_file \
    $bam_1 $bam_2 $bam_3
assert_in_stdout "img/DEL_chr4_115928726-115931880.png"
assert_in_stdout "img/DUP_chrX_101055330-101067156.png"
assert_no_stderr
assert_exit_code 0
