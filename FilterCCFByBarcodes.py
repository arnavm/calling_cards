import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type = str, required = True, help = "Input CCF file")
parser.add_argument("-b", "--barcodes", type = str, required = True, help = "List of barcode (one per line)")
parser.add_argument("-o", "--output", type = str, required = True, help = "Output CCF file")
args = parser.parse_args()

if __name__ == "__main__":
    ccf = pd.read_table(args.input, header = None)
    with open(args.barcodes, 'r') as f:
        barcodes = set([line.strip() for line in f])
    filtered = ccf[ccf[5].isin(barcodes)]
    filtered.to_csv(args.output, sep = '\t', header = False, index = False)

