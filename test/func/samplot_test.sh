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
out_file_name="test.jpg"

rm -f $out_file_name
run basic_operation \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1,$bam_2,$bam_3 \
        -o $out_file_name \
        -t $sv_type
assert_exit_code 0
assert_equal $out_file_name $( ls $out_file_name )
