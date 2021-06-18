#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

STOP_ON_FAIL=0
data_path="test/data/"
func_path="test/func/"

bam_1=$data_path"NA12878_restricted.bam"
bam_2=$data_path"NA12889_restricted.bam"
bam_3=$data_path"NA12890_restricted.bam"


vcf_file=$data_path"test.vcf"
cmd_file=$func_path"test.cmd"
test_dir=$func_path"test_vcf_dir"
rm -f $cmd_file
rm -rf $test_dir
run from_vcf \
    samplot vcf \
        -d $test_dir \
        --vcf $vcf_file \
        --sample_ids HG002 HG003 HG004 \
        -b $data_path"HG002_Illumina.bam" \
        $data_path"HG003_Illumina.bam" \
        $data_path"HG004_Illumina.bam" \
        --manual_run\
        --command_file $cmd_file
if [ $from_vcf ]; then
    assert_no_stderr
    assert_exit_code 0
    assert_equal $test_dir/index.html $( ls $test_dir/index.html )
    assert_equal $cmd_file $( ls $cmd_file )
fi
rm -f $cmd_file
rm -rf $test_dir

vcf_file=$data_path"test.vcf"
cmd_file=$func_path"test.cmd"
test_dir=$func_path"test_vcf_gff3_dir"
rm -f $cmd_file
rm -rf $test_dir
run from_vcf_gff3 \
    samplot vcf \
        -d $test_dir \
        --vcf $vcf_file \
        --sample_ids HG002 HG003 HG004 \
        -b $data_path"HG002_Illumina.bam" \
        $data_path"HG003_Illumina.bam" \
        $data_path"HG004_Illumina.bam" \
        --gff3 $data_path"Homo_sapiens.GRCh37.82.sort.2_X.gff3.gz"\
        --manual_run\
        --command_file $cmd_file
if [ $from_vcf_gff3 ]; then
    assert_no_stderr
    assert_exit_code 0
    assert_equal $test_dir/index.html $( ls $test_dir/index.html )
    assert_equal $cmd_file $( ls $cmd_file )
fi
rm -f $cmd_file
rm -rf $test_dir


vcf_file=$data_path"test.vcf"
cmd_file=$func_path"test.cmd"
test_dir=$func_path"test_vcf_gff3_dir"
rm -f $cmd_file
rm -rf $test_dir
run from_vcf_annotated \
    samplot vcf \
        -d $test_dir \
        --vcf $vcf_file \
        --sample_ids HG002 HG003 HG004 \
        -b $data_path"HG002_Illumina.bam" \
        $data_path"HG003_Illumina.bam" \
        $data_path"HG004_Illumina.bam" \
        -T $data_path"Homo_sapiens.GRCh37.82.sort.2_X.gff3.gz"\
        -A $data_path"Alu.2_X.bed.gz" \
        --manual_run\
        --command_file $cmd_file
if [ $from_vcf_annotated ]; then
    assert_no_stderr
    assert_exit_code 0
    assert_equal $test_dir/index.html $( ls $test_dir/index.html )
    assert_equal $cmd_file $( ls $cmd_file )
fi
rm -f $cmd_file
rm -rf $test_dir


vcf_file=$data_path"test.vcf"
cmd_file=$func_path"test.cmd"
test_dir=$func_path"test_vcf_auto_dir"
rm -rf $test_dir
run from_vcf_auto \
    samplot vcf \
        -d $test_dir \
        --vcf $vcf_file \
        --sample_ids HG002 HG003 HG004 \
        -b $data_path"HG002_Illumina.bam" \
        $data_path"HG003_Illumina.bam" \
        $data_path"HG004_Illumina.bam" 
if [ $from_vcf_auto ]; then
    assert_in_stderr "Window size is under 1.5x the estimated fragment length and will be resized to 847. Rerun with -w 604 to override"
    assert_exit_code 0
    assert_equal $test_dir/index.html $( ls $test_dir/index.html )
    assert_equal $test_dir/DEL_1_24804397_24807302.png $( ls $test_dir/DEL_1_24804397_24807302.png )
    assert_equal $test_dir/DUP_4_99813786_99817098.png $( ls $test_dir/DUP_4_99813786_99817098.png )
    assert_equal $test_dir/DUP_11_67974431_67975639.png $( ls $test_dir/DUP_11_67974431_67975639.png )
    assert_equal $test_dir/INV_12_12544867_12546613.png $( ls $test_dir/INV_12_12544867_12546613.png )
    assert_equal $test_dir/DEL_19_12694866_12698924.png $( ls $test_dir/DEL_19_12694866_12698924.png )
    assert_equal $test_dir/TRA_1_24804398_43059290.png $( ls $test_dir/TRA_1_24804398_43059290.png )
    assert_equal $test_dir/TRA_1_24804399_99813787.png $( ls $test_dir/TRA_1_24804399_99813787.png )
fi
rm -f $cmd_file
rm -rf $test_dir

vcf_file=$data_path"test.vcf"
cmd_file=$func_path"test.cmd"
test_dir=$func_path"test_plotall_dir"
rm -f $cmd_file
rm -rf $test_dir
run plot_all \
    samplot vcf \
        -d $test_dir \
        --vcf $vcf_file \
        --sample_ids HG002 HG003 HG004 \
        -b $data_path"HG002_Illumina.bam" \
        $data_path"HG003_Illumina.bam" \
        $data_path"HG004_Illumina.bam" \
        --plot_all
if [ $plot_all ]; then
    assert_in_stderr "Window size is under 1.5x the estimated fragment length and will be resized to 847. Rerun with -w 604 to override"
    assert_exit_code 0
    assert_equal "$test_dir/index.html" $( ls $test_dir/index.html )
    assert_equal "$test_dir/DEL_19_12694866_12698924.png" $( ls "$test_dir/DEL_19_12694866_12698924.png" )
    assert_equal "$test_dir/DUP_4_99813786_99817098.png" $( ls "$test_dir/DUP_4_99813786_99817098.png" )
    assert_equal "$test_dir/DUP_4_99813786_99817098.png" $( ls "$test_dir/DUP_4_99813786_99817098.png" )
    assert_equal "$test_dir/TRA_1_24804398_43059290.png" $( ls $test_dir/TRA_1_24804398_43059290.png )
    assert_equal "$test_dir/TRA_1_24804399_99813787.png" $( ls $test_dir/TRA_1_24804399_99813787.png )
    assert_equal "$test_dir/DEL_1_24804397_24807302.png" $( ls "$test_dir/DEL_1_24804397_24807302.png" )
    assert_equal "$test_dir/DUP_11_67974431_67975639.png" $( ls "$test_dir/DUP_11_67974431_67975639.png" )
    assert_equal "$test_dir/INV_12_12544867_12546613.png" $( ls "$test_dir/INV_12_12544867_12546613.png" )

fi
rm -f $cmd_file
rm -rf $test_dir

vcf_file=$data_path"test.vcf"
cmd_file="test.cmd"
test_dir="test_vcf_dir"
ped_file=$data_path"test.ped"

run denovo_only_noped \
    samplot vcf \
        -d $test_dir \
        --vcf $vcf_file \
        --sample_ids HG002 HG003 HG004 \
        -b $data_path"HG002_Illumina.bam" \
        $data_path"HG003_Illumina.bam" \
        $data_path"HG004_Illumina.bam" \
        --dn_only\
if [ $denovo_only_noped ]; then
    assert_in_stderr "Missing --ped, required when using --dn_only"
fi
rm -f $cmd_file
rm -rf $test_dir

vcf_file=$data_path"test.vcf"
cmd_file="test.cmd"
test_dir="test_vcf_dir"
ped_file=$data_path"test.ped"

run denovo_only \
    samplot vcf \
        -d $test_dir \
        --sample_ids HG002 HG003 HG004 \
        --vcf $vcf_file \
        -b $data_path"HG002_Illumina.bam" \
        $data_path"HG003_Illumina.bam" \
        $data_path"HG004_Illumina.bam" \
        --dn_only\
        --ped $data_path"test.ped"
if [ $denovo_only ]; then
    assert_no_stderr
    assert_exit_code 0
    assert_equal "$test_dir/DEL_19_12694867_12698924.png" $( ls "$test_dir/DEL_19_12694866_12698924.png" )
    assert_equal "ls: test_vcf_dir/DUP_4_99813786_99817098.png: No such file or directory" $( ls "$test_dir/DUP_4_99813786_99817098.png" )
    assert_equal "" $( ls "$test_dir/TRA_1_24804399_43059290.png" )
    assert_equal "" $( ls "$test_dir/TRA_1_24804398_99813787.png" )
    assert_equal "" $( ls "$test_dir/DEL_1_24804397_24807302.png" )
fi

rm -rf ssshtest
