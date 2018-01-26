#!/bin/bash

SAMPLOT=`which samplot 2> /dev/null`
BCFTOOLS=`which bcftools 2> /dev/null`
outdir=$( pwd )
vcf=

set -eu

usage()
{
    cat << EOF
    usage: `basename $0` OPTIONS -v svs.vcf sample_1.bam sample_2.bam ...

    General options:
    -h      Show this message
    -o      Output directory ($outdir)

    Path options:
    -B      BCFTOOLS path ($BCFTOOLS)
    -S      SAMPLOT path ($SAMPLOT)
EOF
}

while getopts "h o:B:S:v:" OPTION; do
case $OPTION in
    h)
        usage
        exit 1
        ;;
    o)
        outdir=$OPTARG
        ;;
    v)
        vcf=$OPTARG
        ;;
    B)
        BCFTOOLS=$OPTARG
        ;;
    S)
        SAMPLOT=$OPTARG
        ;;
    ?)
        usage
        exit
        ;;
    esac
done

if [ -z "$SAMPLOT" ]; then
    echo "ERROR: samplot not found. Set with -S"
    usage
    exit 1
fi

if [ -z "$BCFTOOLS" ]; then
    echo "ERROR: bcftools not found. Set with -B"
    usage
    exit 1
fi

if [ -z "$vcf" ]; then
    echo "ERROR: no vcf given. Set with -v"
    usage
    exit 1
fi

bams=""
tmp=$@
opts=(${tmp// / })
for bam in ${opts[@]:$((OPTIND - 1))}; do
    if [ -z "$bams" ]; then
        bams="$bam"
    else
        bams+=",$bam"
    fi
done

if [ -z "$bams" ]; then
    echo "ERROR: no bams given"
    usage
    exit 1
fi

IFS=$'\n' 
for sv in `$BCFTOOLS view -i 'SVTYPE="DEL" || SVTYPE="DUP" || SVTYPE="INV" || SVTYPE="INS"' $vcf | $BCFTOOLS query -f "%CHROM %POS %INFO/END %INFO/SVTYPE\n"`; do
        IFS=$' '
        arr=($sv)

        $SAMPLOT -c ${arr[0]} -s ${arr[1]} -e ${arr[2]} -t ${arr[3]} -o ${outdir}/${arr[3]}\_${arr[0]}\_${arr[1]}\-${arr[2]}\.png  -b $bams -a \
        > ${outdir}/${arr[3]}\_${arr[0]}\_${arr[1]}\-${arr[2]}\.args

        echo ${outdir}/${arr[3]}\_${arr[0]}\_${arr[1]}\-${arr[2]}\.png

        IFS=$'\n' 
done
