#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse

parser = argparse.ArgumentParser(description="uses a bed file to build commands for samplot to create images")
parser.add_argument("-b", "--bed_file", help="bed file input with regions to show in images",required=True)
args = parser.parse_args()

with open (args.bed_file, 'r') as bed: 
    for line in bed:
        fields = line.strip().split()[:3]
        print (fields)
        sys.exit(0)
