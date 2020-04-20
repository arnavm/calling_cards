#!/usr/bin/env python
# coding: utf-8

# This program takes as input an aligned BAM file of calling card
# sequencing data. It reads the XI and GS tags and outputs data into 
# calling card format, which is suitable for visualization.

import argparse
from collections import Counter
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--barcode", type=str, nargs='+', required=False)
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

def getBarcodeString(read, tags):
    return '/'.join([read.get_tag(tag) for tag in tags])

if __name__ == "__main__":
    # Open input and output files
    bamfile = pysam.AlignmentFile(args.input, "rb")
    outfile = open(args.output, "w")
    insertions = Counter()
    barcodeTags = args.barcode
    for read in bamfile:
        if barcodeTags:
            try:
                barcode = getBarcodeString(read, barcodeTags)
                if read.get_tag("XI").endswith('-') or read.get_tag("XI").endswith('+'):
                    insertions.update(['|'.join([read.get_tag("XI"), barcode])])
                else:
                    insertions.update(['|'.join([read.get_tag("XI"), read.get_tag("GS"), barcode])])
            except KeyError:
                continue
        else:
            if read.get_tag("XI").endswith('-') or read.get_tag("XI").endswith('+'):
                insertions.update([read.get_tag("XI")])
            else:
                insertions.update(['|'.join([read.get_tag("XI"), read.get_tag("GS")])])
    for k, v in insertions.items():
        line = k.split('|')
        line.insert(3, str(v))
        print('\t'.join(line), file=outfile)
    # Close files
    outfile.close()
    bamfile.close()
