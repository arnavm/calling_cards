#!/usr/bin/env python
# coding: utf-8

# This program takes as input an aligned BAM file of PB calling card reads.
# It then calculates the insertion site of each read and adds that information
# to each read with the tag "XI." Specifically, this program calculates and 
# the coordinates of the bases immediately preceding the read and stores
# that data in the "XI" tag. The program also includes the sequence of those
# preceding bases in the "XZ" tag. This program is compatible with
# piggyBac and Sleeping Beauty transpositions from calling card experiments,
# as well as Tn5 insertions from ATAC-seq data.

import argparse
from functools import lru_cache
import pysam
import sys
from twobitreader import TwoBitFile
import warnings

transposaseMotifs = {}
transposaseMotifs["PB"] = "TTAA"
transposaseMotifs["SB"] = "TA"
transposaseMotifs["Tn5"] = True

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--transposase", type=str, required=True, choices=["PB", "SB", "Tn5"])
parser.add_argument("-f", "--filter", action="store_true")
parser.add_argument("input", type=str)
parser.add_argument("reference", type=str)
parser.add_argument("output", type=str)
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
    # Special consideration for transposase like Tn5
    if transposase == "Tn5":
        insertSiteLength = 9
    else:
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
            # The span of the bases preceding the insertion are:
            # [insert, insert + insertSiteLength)
            # Special consideration for Tn5, where we offset from the
            # read end to map the transposase binding site.
            if transposase == "Tn5":
                tags["XI"] = "{0}|{1}|{2}|-".format(read.reference_name, insert - insertSiteLength, insert)
                tags["XZ"] = fetchGenomicSequence(read.reference_name, insert - insertSiteLength, insert)
            else:
                tags["XI"] = "{0}|{1}|{2}|-".format(read.reference_name, insert, insert + insertSiteLength)
                tags["XZ"] = fetchGenomicSequence(read.reference_name, insert, insert + insertSiteLength)
        else:
            # Read is in the forward orientation
            # The insertion point is at the start of the alignment
            insert = read.reference_start
            # If the read is soft clipped at the beginning, correct the insertion coordinate
            if read.cigartuples[0][0] == 4:
                insert -= read.cigartuples[0][1]
            # The span of the bases preceding the insertions are:
            # [insert - insertSiteLength, insert)
            # Special consideration for Tn5, as above
            if transposase == "Tn5":
                tags["XI"] = "{0}|{1}|{2}|+".format(read.reference_name, insert, insert + insertSiteLength)
                tags["XZ"] = fetchGenomicSequence(read.reference_name, insert, insert + insertSiteLength)
            else:
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
            try:
                tags = getInsertionSiteTags(read)
            except ValueError:
                # This catches alignments at the end of the mitochondrial genome
                continue
            # Add the insertion site to the read
            read.set_tag("XI", tags["XI"], "Z")
            read.set_tag("XZ", tags["XZ"], "Z")
            # Set the strand tag
            if read.flag & 0x10:
                read.set_tag("GS", "-", "Z")
            else:
                read.set_tag("GS", "+", "Z")
            # Set the query name tag (first field of the FASTQ record)
            read.set_tag("QN", read.query_name.split(' ')[0], "Z")
            # If we are filtering reads, write read only if the bases 
            # preceding the read match the transposase insertion site motif
            if args.filter:
                # Filtering should not be used with a motif-agnostic transposase
                # like Tn5, as there is no way to validate the insertion event
                if args.transposase == "Tn5":
                    warnings.warn("The filter flag should not be used with Tn5/ATAC-seq data", UserWarning)
                if tags["XZ"] == insertionMotif:
                    outfile.write(read)
            # Otherwise, write all reads to the output file
            else:
                outfile.write(read)
            
    # Close files
    outfile.close()
    bamfile.close()
