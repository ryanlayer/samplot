#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

STOP_ON_FAIL=0
data_path="test/data/"
func_path="test/func/"

bam_1=$data_path"NA12878_restricted.bam"
bam_2=$data_path"NA12889_restricted.bam"
bam_3=$data_path"NA12890_restricted.bam"

sv_chrm=chr4
sv_start=115928730
sv_end=115931875
sv_type=DEL
out_file_name=$func_path"test_del.png"

rm -f $out_file_name
run basic_operation \
    samplot plot \
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

sv_chrm=chr4
sv_start=115928730
sv_end=115931875
sv_type=DEL
out_file_name=$func_path"test_max_coverage.png"

rm -f $out_file_name
run max_coverage \
    samplot plot \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3 \
        -o $out_file_name \
        --max_coverage 50\
        -t $sv_type
if [ $max_coverage ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sv_chrm=chr4
sv_start=115928730
sv_end=115931875
sv_type=DEL
out_file_name=$func_path"test_coverage_only.png"

rm -f $out_file_name
run coverage_only \
    samplot plot \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 \
        -o $out_file_name \
        --coverage_only \
        -t $sv_type
if [ $coverage_only ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

out_file_name=$func_path"test_same_yaxis.png"

sv_chrm=chrX
sv_start=101055330
sv_end=101067156
sv_type=DUP
rm -f $out_file_name
run same_yaxis \
    samplot plot \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam_1 $bam_2 $bam_3\
        -o $out_file_name \
        -t $sv_type \
        --same_yaxis_scales
if [ $basic_operation ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi

sv_chrm=chr4
sv_start=115928730
sv_end=115931875
sv_type=DEL
out_file_name=$func_path"test_zoom.png"
rm -f $out_file_name
run basic_operation_zoom \
    samplot plot \
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

sample_out_file_name=$func_path"sample.png"
run sampling_normal \
    samplot plot\
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

sample_out_file_name=$func_path"sample_zoom.png"
run sampling_normal_zoom \
    samplot plot \
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
out_file_name=$func_path"dup.png"
rm -f $out_file_name

run common_insert_size_scale \
    samplot plot\
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

out_file_name=$func_path"dup_zoom.png"
rm -f $out_file_name
run common_insert_size_scale_zoom \
    samplot plot \
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


out_file_name=$func_path"no_sv_type.png"
rm -f $out_file_name

run no_sv_type \
    samplot plot \
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

sv_chrm=X
sv_start=101055330
sv_end=101067156
sv_type=DUP
out_file_name=$func_path"longread_nanopore_dup.png"
bam=$data_path"nanopore-NA12878.bam"
rm -f $out_file_name

run nanopore_dup \
    samplot plot \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $nanopore_dup ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

out_file_name=$func_path"longread_nanopore_dup_zoom.png"
rm -f $out_file_name
run nanopore_dup_zoom \
    samplot plot\
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
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

sv_chrm=4
sv_start=115928730
sv_end=115931875
sv_type=DEL
out_file_name=$func_path"longread_nanopore_del.png"
bam=$data_path"nanopore-NA12878.bam"
rm -f $out_file_name
run nanopore_del \
    samplot plot\
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $nanopore_del ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

out_file_name=$func_path"longread_nanopore_del_zoom.png"
rm -f $out_file_name
run nanopore_del_zoom \
    samplot plot\
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
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

sv_chrm=chr1
sv_start=58343117
sv_end=58343622
sv_type=DEL
out_file_name=$func_path"longread_del.png"
bam=$data_path"hg19_chr1_58343117_58343622_deletion.bam"
rm -f $out_file_name
run longread_del \
    samplot plot \
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $longread_del ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

out_file_name=$func_path"longread_del_zoom_big_zoom.png"
rm -f $out_file_name
run longread_del_zoom_big_zoom \
    samplot plot\
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10  \
        --zoom 500
if [ $longread_del_zoom_big_zoom ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_in_stderr "Insufficient reads for fragment length estimate."
    assert_no_stdout
fi


out_file_name=$func_path"longread_del_zoom_zoom.png"
rm -f $out_file_name
run longread_del_zoom_zoom \
    samplot plot\
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
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

sv_chrm=chr21
sv_start=27373431
sv_end=27375410
sv_type=INV
out_file_name=$func_path"longread_inv.png"
bam=$data_path"hg19_chr21_27373431_27375410_inversion.bam"
rm -f $out_file_name
run longread_inv \
    samplot plot\
        -c $sv_chrm -s $sv_start -e $sv_end \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -d 10 
if [ $longread_inv ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

out_file_name=$func_path"longread_inv_zoom.png"
rm -f $out_file_name
run longread_inv_zoom \
    samplot plot\
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
    assert_in_stderr "Insufficient reads for fragment length estimate."
fi

sv_chrm=1
sv_start=89475845
sv_end=89478561
sv_type=DEL
out_file_name=$func_path"linkedread_del.png"
bam=$data_path"HG002_1_89475845-89478561_DEL.tenx.bam"
rm -f $out_file_name
run linkedread_del \
    samplot plot\
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

out_file_name=$func_path"linkedread_del_zoom.png"
rm -f $out_file_name
run linkedread_del_zoom \
    samplot plot\
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
out_file_name=$func_path"translocation.png"
bam=$data_path"2_59305747-59505747_X_151018513-151218513.BND.bam"
run translocation \
    samplot plot\
        -c $sv_chrm_1 -s $sv_start_1 -e $sv_end_1 \
        -c $sv_chrm_2 -s $sv_start_2 -e $sv_end_2 \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -A $data_path"Alu.2_X.bed.gz" \
        -T $data_path"Homo_sapiens.GRCh37.82.sort.2_X.gff3.gz" \
        --zoom 10000
if [ $translocation ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi


out_file_name=$func_path"csi-annotation.png"
bam=$data_path"2_59305747-59505747_X_151018513-151218513.BND.bam"
run translocation \
    samplot plot\
        -c $sv_chrm_1 -s $sv_start_1 -e $sv_end_1 \
        -c $sv_chrm_2 -s $sv_start_2 -e $sv_end_2 \
        -b $bam \
        -o $out_file_name \
        -t $sv_type \
        -A $data_path"Alu.2_X.csionly.bed.gz" \
        -T $data_path"Homo_sapiens.GRCh37.csionly.2_X.gff3.gz" \
        --zoom 10000
if [ $translocation ]; then
    assert_exit_code 0
    assert_equal $out_file_name $( ls $out_file_name )
    assert_no_stdout
    assert_no_stderr
fi


rm -rf $func_path"img/" ssshtest
