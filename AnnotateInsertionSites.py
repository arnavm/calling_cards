#!/usr/bin/env python
# coding: utf-8

# This program takes as input an aligned BAM file of PB calling card reads.
# It then calculates the insertion site of each read and adds that information
# to each read with the tag "XI." Specifically, this program calculates and 
# the coordinates of the bases immediately preceding the read and stores
# that data in the "XI" tag. The program also includes the sequence of those
# preceding bases in the "XZ" tag. This program is compatible with
# piggyBac and SleepingBeauty transpositions.

import argparse
from functools import lru_cache
import pysam
import sys
from twobitreader import TwoBitFile

transposaseMotifs = {}
transposaseMotifs["PB"] = "TTAA"
transposaseMotifs["SB"] = "TA"

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--transposase", type = str, required = True, choices = ["PB", "SB", "HelR"])
parser.add_argument("-f", "--filter", action = "store_true")
parser.add_argument("input", type = str)
parser.add_argument("reference", type = str)
parser.add_argument("output", type = str)
args = parser.parse_args()

ref = TwoBitFile(args.reference)

# This is to cache the most recent 1024 genomic loci looked up. This saves
# overhead by not having to needlessly query the reference
@lru_cache(maxsize = 1024)
def fetchGenomicSequence(chromosome, start, end):
    return ref[chromosome][start:end].upper()

# This function returns a function specifically tailored to the transposase
# specified by the user. This practice is known as "currying" and prevents
# having to re-initialize the insertSiteLength parameter during each iteration
def makeInsertionSiteFunction(transposase):
    insertSiteLength = len(transposaseMotifs[transposase])
    def insertionSiteTags(read):
        tags = {}
        # Check the orientation of the read:
        if read.flag & 0x10:
            # Read is reverse complemented
            # The insertion point is at the end of the alignment
            insert = read.reference_end
            # If the read is soft clipped at the end, correct the insertion coordinate
            if read.cigartuples[-1][0] == 4:
                insert += read.cigartuples[-1][1]
            # The span of the 4 bases preceding the insertion are:
            # [insert, insert + 4)
            tags["XI"] = "{0}|{1}|{2}|-".format(read.reference_name, insert, insert + insertSiteLength)
            tags["XZ"] = fetchGenomicSequence(read.reference_name, insert, insert + insertSiteLength)
        else:
            # Read is in the forward orientation
            # The insertion point is at the start of the alignment
            insert = read.reference_start
            # If the read is soft clipped at the beginning, correct the insertion coordinate
            if read.cigartuples[0][0] == 4:
                insert -= read.cigartuples[0][1]
            # The span of the 4 bases preceding the insertions are:
            # [insert - 4, insert)
            tags["XI"] = "{0}|{1}|{2}|+".format(read.reference_name, insert - insertSiteLength, insert)
            tags["XZ"] = fetchGenomicSequence(read.reference_name, insert - insertSiteLength, insert)
        return tags
    return insertionSiteTags


if __name__ == "__main__":
    # Open files
    bamfile = pysam.AlignmentFile(args.input, "rb")
    outfile = pysam.AlignmentFile(args.output, "wb", template = bamfile)
    # Create the function for analyzing insertion sites
    getInsertionSiteTags = makeInsertionSiteFunction(args.transposase)
    insertionMotif = transposaseMotifs[args.transposase]
    for read in bamfile:
        # Check that the read is mapped
        if read.flag & 0x4:
            # Unmapped read, do nothing
            pass
            # outfile.write(read)
        else:
            # Mapped read
            # Calculate the insertion site
            tags = getInsertionSiteTags(read)
            # Add the insertion site to the read
            read.set_tag("XI", tags["XI"], "Z")
            read.set_tag("XZ", tags["XZ"], "Z")
            # Set the strand tag
            if read.flag & 0x10:
                read.set_tag("GS", "-", "Z")
            else:
                read.set_tag("GS", "+", "Z")
            # If we are filtering reads, write read only if the bases 
            # preceding the read match the transposase insertion site motif
            if args.filter:
                if tags["XZ"] == insertionMotif:
                    outfile.write(read)
            # Otherwise, write all reads to the output file
            else:
                outfile.write(read)
            
    # Close files
    outfile.close()
    bamfile.close()
