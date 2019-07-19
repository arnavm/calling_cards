#!/usr/bin/env python
# coding: utf-8

# This program takes as input a BAM file of aligned reads.
# It then tags each read with the supplied tag and writes a new BAM file.

import argparse
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", type = str)
parser.add_argument("output", type = str)
parser.add_argument("-t", "--tag", type = str, required = True)
args = parser.parse_args()

if __name__ == "__main__":
    # Open files
    bamfile = pysam.AlignmentFile(args.input, "rb")
    outfile = pysam.AlignmentFile(args.output, "wb", template = bamfile)
    tag = args.tag.split(':')
    for read in bamfile:
        # Add the species tag 
        read.set_tag(tag[0], tag[2], tag[1])
        outfile.write(read)
    # Close files
    outfile.close()
    bamfile.close()
