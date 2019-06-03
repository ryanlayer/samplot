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
        -b $bam_1 $bam_2 $bam_3 \
        -o $out_file_name \
        -t $sv_type
if [ $basic_operation ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name="test_zoom.png"
rm -f $out_file_name
run basic_operation_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3 \
        -o $out_file_name \
        -t $sv_type \
        --zoom 500
if [ $basic_operation_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sample_out_file_name="sample.png"
run sampling_normal \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3 \
        -o $sample_out_file_name \
        -t $sv_type \
        -d 10
if [ $sampling_normal ]; then
    assert_exit_code 0
    assert_equal $sample_out_file_name $( ls $sample_out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sample_out_file_name="sample_zoom.png"
run sampling_normal_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3 \
        -o $sample_out_file_name \
        -t $sv_type \
        -d 10 \
        --zoom 500
if [ $sampling_normal_zoom ]; then
    assert_exit_code 0
    assert_equal $sample_out_file_name $( ls $sample_out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sv_chrm=chrX
sv_start=101055330
sv_end=101067156
sv_type=DUP
out_file_name="dup.png"
rm -f $out_file_name

run common_insert_size_scale \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3 \
        -o $out_file_name \
        -t $sv_type \
        -d 10 \
        --common_insert_size
if [ $common_insert_size_scale ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name="dup_zoom.png"
rm -f $out_file_name
run common_insert_size_scale_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3 \
        -o $out_file_name \
        -t $sv_type \
        -d 10 \
        --zoom 500 \
        --common_insert_size
if [ $common_insert_size_scale_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi


out_file_name="no_sv_type.png"
rm -f $out_file_name

run no_sv_type \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3 \
        -o $out_file_name \
        -d 10 \
        --common_insert_size
if [ $no_sv_type ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

rm -rf img/
mkdir img
vcf_file=../data/NA12878.trio.svt.subset.vcf
run from_vcf \
    ../../src/samplot_vcf.sh \
    -d 10 \
    -o img \
    -v $vcf_file \
    $bam_1 $bam_2 $bam_3
if [ $from_vcf ]; then
    assert_in_stdout "img/DEL_chr4_115928726-115931880.png"
    assert_in_stdout "img/DUP_chrX_101055330-101067156.png"
    assert_no_stderr
    assert_exit_code 0
fi

sv_chrm=X
sv_start=101055330
sv_end=101067156
sv_type=DUP
out_file_name="longread_nanopore_dup.png"
bam=../data/nanopore-NA12878.bam
rm -f $out_file_name

run nanopore_dup \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $nanopore_dup ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name="longread_nanopore_dup_zoom.png"
rm -f $out_file_name
run nanopore_dup_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10  \
        --zoom 1000
if [ $nanopore_dup_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sv_chrm=4
sv_start=115928730
sv_end=115931875
sv_type=DEL
out_file_name="longread_nanopore_del.png"
bam=../data/nanopore-NA12878.bam
rm -f $out_file_name
run nanopore_del \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $nanopore_del ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name="longread_nanopore_del_zoom.png"
rm -f $out_file_name
run nanopore_del_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10  \
        --zoom 500
if [ $nanopore_del_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sv_chrm=chr1
sv_start=58343117
sv_end=58343622
sv_type=DEL
out_file_name="longread_del.png"
bam=../data/hg19_chr1_58343117_58343622_deletion.bam
rm -f $out_file_name
run longread_del \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $longread_del ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name="longread_del_zoom_big_zoom.png"
rm -f $out_file_name
run longread_del_zoom_big_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10  \
        --zoom 500
if [ $longread_del_zoom_big_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_in_stderr "Ignoring zoom command."
    assert_no_stdout
fi


out_file_name="longread_del_zoom_zoom.png"
rm -f $out_file_name
run longread_del_zoom_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10  \
        --zoom 200
if [ $longread_del_zoom_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sv_chrm=chr21
sv_start=27373431
sv_end=27375410
sv_type=INV
out_file_name="longread_inv.png"
bam=../data/hg19_chr21_27373431_27375410_inversion.bam
rm -f $out_file_name
run longread_inv \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $longread_inv ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name="longread_inv_zoom.png"
rm -f $out_file_name
run longread_inv_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10  \
        --zoom 750
if [ $longread_inv_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sv_chrm=1
sv_start=89475845
sv_end=89478561
sv_type=DEL
out_file_name="linkedread_del.png"
bam=../data/HG002_1_89475845-89478561_DEL.tenx.bam
rm -f $out_file_name
run linkedread_del \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $linkedread_del ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name="linkedread_del_zoom.png"
rm -f $out_file_name
run linkedread_del_zoom \
    python ../../src/samplot.py \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 \
        --zoom 500
if [ $linkedread_del_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi


sv_chrm_1=2
sv_start_1=59405943
sv_end_1=59405943
sv_chrm_2=X
sv_start_2=151118533
sv_end_2=151118533
sv_type=BND
out_file_name="translocation.png"
bam=../data/2_59305747-59505747_X_151018513-151218513.BND.bam
run translocation \
    python ../../src/samplot.py \
        -c $sv_chrm_1 -s $sv_start_1 -e $sv_end_1 \
        -c $sv_chrm_2 -s $sv_start_2 -e $sv_end_2 \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -A ../data/Alu.2_X.bed.gz \
        -T ../data/Homo_sapiens.GRCh37.82.sort.2_X.gff3.gz \
        --zoom 10000
if [ $translocation ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi
