#!/usr/bin/env python
# coding: utf-8

# This program takes as input a list of barcode files (one per line)
# as well as a list of files to write output to. For each
# input file, the program will identify which barcodes are unique to that
# file, and those will be written to the respective output file.

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type = str, nargs = '+', required = True)
parser.add_argument("-o", "--output", type = str, nargs = '+', required = True)
args = parser.parse_args()

if __name__ == "__main__":
    # Check that the number of output files matches the number of input files
    assert(len(args.input) == len(args.output))
    # Open input files, store each file as a list of barcodes
    barcodes = [[line.strip() for line in open(f, 'r').readlines()] for f in args.input]
    totalInitial = 0
    for s in barcodes:
        totalInitial += len(s)
    # Create the set of observed barcodes that only appear once across all input files
    b, c = np.unique(np.concatenate(barcodes), return_counts = True)
    uniqueBarcodes = set(b[c == 1])
    # For each set of input barcodes, find those barcodes that are in uniqueBarcodes
    uniques = [set(s).intersection(uniqueBarcodes) for s in barcodes]
    # Write each set of unique barcodes to their respective output files
    assert(len(uniques) == len(args.output)) # Sanity check
    totalFinal = 0
    for u, f in zip(uniques, args.output):
        with open(f, 'w') as out:
            out.write('\n'.join(sorted(list(u))))
            out.write('\n')
        totalFinal += len(u)
    # Report what percent of input barcodes written to file
    print("Wrote {} out of {} barcodes ({:.2f}%)".format(totalFinal, totalInitial, totalFinal / totalInitial * 100))
