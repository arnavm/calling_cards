import argparse
import pandas as pd
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type = str, required = True, help = "Input BAM file")
parser.add_argument("-b", "--barcodes", type = str, required = True, help = "List of barcodes (one per line)")
parser.add_argument("-o", "--output", type = str, required = True, help = ("Output BAM file"))
args = parser.parse_args()


if __name__ == "__main__":
    bamfile = pysam.AlignmentFile(args.input, 'rb')
    with open(args.barcodes, 'r') as f:
        barcodes = set([line.strip() for line in f])
    outfile = pysam.AlignmentFile(args.output, 'wb', template = bamfile)
    for read in bamfile:
        try:
            if read.get_tag("CB") in barcodes:
                outfile.write(read)
        except KeyError:
            continue
    outfile.close()
