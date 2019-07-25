#!/usr/bin/env python

# This script performs a two-hit filter on single cell calling card reads.
# It takes as one BAM file, the circularized Illumina scCC library. It outputs 
# a BAM file of reads filtered by the following criterion: for a given 
# cell barcode and insertion site, there must be at least two different UMIs
# corresponding to that insertion.

import argparse
from collections import defaultdict
from datetime import datetime
import pybedtools as pb
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--platform", type = str, required = True, choices = ["DS", "10x"])
parser.add_argument("-i", "--input", type = str, required = True, help = "Reads from scCC library")
parser.add_argument("-o", "--output", type = str, required = True, help = "Output filename")
parser.add_argument("-v", "--verbose", action = "store_true")
args = parser.parse_args()

if args.platform == "DS":
    barcodeTag = "XC"
    UMITag = "XM"
elif args.platform == "10x":
    barcodeTag = "CB"
    UMITag = "UB"
else:
    print("Unknown platform")
    sys.exit()

totalReads = 0
numSkippedCellBarcode = 0
numSkippedUMI = 0
numSkippedIns = 0
numSkippedSingleUMI = 0

def getTimeString():
    return datetime.now().strftime("%H:%M:%S")

if __name__ == "__main__":
    # Open bamfile
    bamfile = pysam.AlignmentFile(args.input, 'rb')
    # Master dictionary of barcode-insertion site-UMI connections
    cb_ins_umi = defaultdict(lambda: defaultdict(set))
    # Make a first pass through the bamfile
    for read in bamfile:
        totalReads += 1
        # Try to get the read's cell barcode
        try:
            barcode = read.get_tag(barcodeTag)
        except KeyError:
            numSkippedCellBarcode += 1
            continue
        
        # Try to get the read's insertion site
        try:
            insertion = read.get_tag("XI")
        except KeyError:
            numSkippedIns += 1
            continue
        
        # Try to get the read's UMI
        try:
            UMI = read.get_tag(UMITag)
        except KeyError:
            numSkippedUMI += 1
            continue

        # If we've made it this far, the read has all the right info
        cb_ins_umi[barcode][insertion].add(UMI)
    # Close the bamfile
    bamfile.close()

    # Create a dictionary to keep track of which barcode-insertion 
    # pairings are valid (i.e. have at least 2 UMIs per insertion)
    cb_ins = defaultdict(set)
    for barcode, insertions in cb_ins_umi.items():
        for insertion, UMIs in insertions.items():
            if len(UMIs) > 1:
                cb_ins[barcode].add(insertion)

    # Now, make a second pass through the bamfile, writing to the
    # output file only those reads that have a cell barcode 
    # and insertion in cb_ins
    bamfile = pysam.AlignmentFile(args.input, 'rb')
    outfile = pysam.AlignmentFile(args.output, 'wb', template = bamfile)
    for read in bamfile:
        try:
            barcode = read.get_tag(barcodeTag)
        except KeyError:
            continue
        try:
            insertion = read.get_tag("XI")
        except KeyError:
            continue
        try:
            UMI = read.get_tag(UMITag)
        except KeyError:
            continue

        if insertion in cb_ins[barcode]:
            outfile.write(read)
        else:
            numSkippedSingleUMI += 1 
    outfile.close()
    bamfile.close()
    
    # Write summary statistics
    numPassed = totalReads - numSkippedCellBarcode - numSkippedUMI - numSkippedIns - numSkippedSingleUMI
    print("{} reads processed".format(totalReads))
    print("-- {} reads ({:04.2f}%) skipped due to lack of cell barcode".format(numSkippedCellBarcode, numSkippedCellBarcode/totalReads * 100))
    print("-- {} reads ({:04.2f}%) skipped due to lack of UMI".format(numSkippedUMI, numSkippedUMI/totalReads * 100))
    print("-- {} reads ({:04.2f}%) skipped due to lack of insertion site".format(numSkippedIns, numSkippedIns/totalReads * 100))
    print("-- {} reads ({:04.2f}%) skipped due to single UMI".format(numSkippedSingleUMI, numSkippedSingleUMI / totalReads * 100))
    print("-- {} reads ({:04.2f}%) written to output file".format(numPassed, numPassed / totalReads * 100))
