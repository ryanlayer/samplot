#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

STOP_ON_FAIL=0

run basic_operation python ../../src/samplot.py -c 2 -s 89161083 -e 89185670 \
    -b "../data/low_coverage/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.restricted_sv_regions.20121211.bam,../data/low_coverage/NA12889.mapped.ILLUMINA.bwa.CEU.low_coverage.restricted_sv_regions.20130415.bam,../data/low_coverage/NA12890.mapped.ILLUMINA.bwa.CEU.low_coverage.restricted_sv_regions.20130415.bam" -o "test.jpg"
assert_exit_code 0
