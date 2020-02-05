#!/bin/env python
# Generates a normalized begdraph of insertion rate
# for calling card data. The output file can be 
# visualized on genome browsers and resembles
# chromosome ideograms.

import argparse
from functools import reduce
import numpy as np
import pandas as pd
from pybedtools import BedTool
from scipy.stats import poisson

# Set up argument parser
parser = argparse.ArgumentParser(description = "Calculate normalized rates of insertions in a calling card experiment", usage = "%(prog)s EXPERIMENT BLOCKS OUTPUT")
parser.add_argument("experiment", type = str, help = "Experiment CCF file")
parser.add_argument("blocks", type = str, help = "Experiment blocks file")
parser.add_argument("output", type = str, help = "File to write output to")

# Fast method for getting number of lines in a file
# For BED files, much faster than calling len() on file
# From https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == "__main__":
    args = parser.parse_args()
    experiment = BedTool(args.experiment)
    blocks = BedTool(args.blocks)
    hopsPerMillion = file_len(args.experiment) / 10**6 # Insertions per million mapped insertions
    intersect = blocks.intersect(experiment, c = True).to_dataframe()
    # For each block, calculate normalized hops per million mapped insertions, per kilobase
    intersect["normRate"] = (intersect["name"] / hopsPerMillion) / ((intersect["end"] - intersect["start"]) / 1000)
    # Convert DataFrame to BED
    outbed = BedTool.from_dataframe(intersect[["chrom", "start", "end", "normRate"]])
    outbed.saveas(args.output)
