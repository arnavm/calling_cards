#!/bin/env python
# A calling card peak caller using
# Bayesian Blocks.

import argparse
import numpy as np
import pandas as pd
from pybedtools import BedTool
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
import sys

# numpy log10 float max
log10FloatMax = np.log10(np.finfo(np.float64).max)

# Set up argument parser
parser = argparse.ArgumentParser(
        description = "Call peaks from calling card data using Bayesian blocks.",
        usage = "%(prog)s [-c] [-r] [-d distance] [-n minSize] [-x maxSize] (-p pValue | -a alpha -m method) [-i intermediate_file] EXPERIMENT BLOCKS BACKGROUND OUTPUT")
parser.add_argument("experiment", type = str, help = "Experiment CCF file")
parser.add_argument("blocks", type = str, help = "Experiment blocks file")
parser.add_argument("background", type = str, help = "Background CCF file")
parser.add_argument("output", type = str, help = "File to write output to")
parser.add_argument("-c", "--pseudocount", type = float, default = 1, help = "Pseudocount for background regions (default: %(default)s)")
group = parser.add_mutually_exclusive_group(required = True)
group.add_argument("-a", "--alpha", type = float, help = "Alpha for multiple hypothesis correction (must be between 0 and 1)")
group.add_argument("-p", "--pValueCutoff", type = float, help = "–log10-transformed p-value cutoff (NOTE: This is a straight cutoff and will not take into account multiple hypothesis correction!)")
parser.add_argument("-m", "--method", type = str, required = "-a" in sys.argv or "--alpha" in sys.argv, choices = ["bonferroni", "sidak", "holm-sidak", "holm", "simes-hochberg", "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"], default = "bonferroni", help = "If alpha provided, need to specificy method of multiple hypothesis correction.")
parser.add_argument("-d", "--distance", type = int, required = False, help = "Merge features closer than this distance (bp)")
parser.add_argument("-n", "--minSize", type = int, required = False, help = "Report peaks larger than this cutoff (bp)")
parser.add_argument("-x", "--maxSize", type = int, required = False, help = "Report peaks smaller than this cutoff (bp)")
parser.add_argument("-r", "--refine", action = "store_true", default = False)
parser.add_argument("-i", "--intermediate", type = str, required = False, help = "Intermediate file to write verbose output (CSV format)")

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
    # Check whether alpha is valid
    if args.alpha:
        assert(0.0 <= args.alpha <= 1.0)
    experiment = BedTool(args.experiment)
    blocks = BedTool(args.blocks)
    background = BedTool(args.background)
    scaleFactor = file_len(args.experiment) / file_len(args.background)
    
    if args.refine:
        data = blocks.intersect(experiment, wa = True, wb = True)
        df = data.to_dataframe()
        df = df.iloc[:, :6]
        df = df.rename(index = str, columns = {"name": "TTAA_chrom", "score": "TTAA_start", "strand": "TTAA_end"})
        df = df.sort_values(["chrom", "start", "end", "TTAA_chrom", "TTAA_start", "TTAA_end"])
        groups = df.groupby(["chrom", "start", "end"])
        first = groups.nth(0)["TTAA_start"]
        last = groups.nth(-1)["TTAA_end"]
        joined = pd.concat([first, last], axis = 1).reset_index()
        refined = joined[["chrom", "TTAA_start", "TTAA_end"]]
        blocks = BedTool.from_dataframe(refined)

    data = blocks.intersect(experiment, c = True).intersect(background, c = True)
    df = data.to_dataframe()
    df = df.rename(index = str, columns = {"name": "expHops", "score": "bgHops"})
    # Add pseudocount to norm_bgHops and floor of pseudocount to expHops.
    # This preserves log-fold change if bgHops is 0, and keeps the value
    # added to expHops an integer.
    df["expHops"] += np.floor(args.pseudocount)
    df["norm_bgHops"] = df["bgHops"] * scaleFactor + args.pseudocount
    # df["norm_bgHops"] = (df["bgHops"] + args.pseudocount) * scaleFactor
    df["pValue"] = poisson.sf(df["expHops"] - 1, df["norm_bgHops"])
    df["negLog10pValue"] = -np.log10(df["pValue"])
    if args.pValueCutoff:
        outdf = df[df["negLog10pValue"] >= args.pValueCutoff]
    else:
        df["rejected"], df["corrected_pValue"], alphacSidak, alphacBonf = multipletests(df["pValue"], alpha = args.alpha, method = args.method, is_sorted = False, returnsorted = False)
        df["negLog10corrected"] = -np.log10(df["corrected_pValue"])
        outdf = df[df["rejected"]]
    outbed = BedTool.from_dataframe(outdf)
    
    # Write intermediate file
    if args.intermediate:
        df.to_csv(args.intermediate)

    if args.distance:
        outbed = outbed.merge(d = args.distance)
    
    if args.minSize:
        minSize = args.minSize
    else:
        minSize = 0

    if args.maxSize:
        maxSize = args.maxSize
    else:
        maxSize = np.inf
    
    sizeFilterDF = outbed.to_dataframe()
    sizeFilterDF["size"] = sizeFilterDF["end"] - sizeFilterDF["start"]
    sizeFilterDF = sizeFilterDF[sizeFilterDF["size"] <= maxSize]
    sizeFilterDF = sizeFilterDF[sizeFilterDF["size"] >= minSize]
    outbed = BedTool.from_dataframe(sizeFilterDF[["chrom", "start", "end"]])
    outbed.saveas(args.output)
