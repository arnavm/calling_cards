#!/usr/bin/env python
# coding: utf-8

# This program takes as input an aligned BAM file of calling card
# sequencing data. It reads the XI and GS tags and outputs data into 
# qBED format, which is suitable for visualization.

import argparse
from collections import Counter
import os
import pysam
import random
import string
import sqlite3
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--barcode", type=str, nargs='+', required=False)
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
parser.add_argument("-db", "--database", action="store_true") 
args = parser.parse_args()

class tempInsertionDB:
    def __init__(self):
        self.name = "tmp" + randomAlphaNumericString() + ".db"
        self.connection = sqlite3.connect(self.name)
        self.cursor = self.connection.cursor()
        self.cursor.execute("CREATE TABLE insertions (key text, reads int)")
    
    def __del__(self):
        self.connection.commit()
        self.connection.close()
        os.remove(self.name)

    def update(self, key):
        assert type(key) == list
        assert len(key) == 1
        key = key[0]
        self.cursor.execute("SELECT * FROM insertions WHERE key = '{}'".format(key))
        record = self.cursor.fetchone()
        if record:
            self.cursor.execute("UPDATE insertions SET reads = {} WHERE key = '{}'".format(record[1] + 1, record[0]))
        else:
            self.cursor.execute("INSERT INTO insertions VALUES ('{}', {})".format(key, 1))

    def items(self):
        self.cursor.execute("SELECT * FROM insertions")
        return self.cursor.fetchall()


# from https://pynative.com/python-generate-random-string/
def randomAlphaNumericString(stringLength=16):
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join((random.choice(lettersAndDigits) for i in range(stringLength)))

def getBarcodeString(read, tags):
    return '/'.join([read.get_tag(tag) for tag in tags])

if __name__ == "__main__":
    # Open input and output files
    bamfile = pysam.AlignmentFile(args.input, "rb")
    outfile = open(args.output, "w")
    if args.database:
        insertions = tempInsertionDB()
    else:
        insertions = Counter()
    barcodeTags = args.barcode
    for read in bamfile:
        insertion_string = read.get_tag("XI")
        if barcodeTags:
            try:
                barcode = getBarcodeString(read, barcodeTags)
                if insertion_string.endswith('-') or insertion_string.endswith('+'):
                    insertion_key = '|'.join([insertion_string, barcode])
                    insertions.update([insertion_key])
                else:
                    insertion_key = '|'.join([insertion_string, read.get_tag("GS"), barcode])
                    insertions.update([insertion_key])
            except KeyError:
                continue
        else:
            if insertion_string.endswith('-') or insertion_string.endswith('+'):
                insertions.update([insertion_string])
            else:
                insertion_key = '|'.join([insertion_string, raed.get_tag("GS")])
                insertions.update([insertion_key])
        print(sys.getsizeof(insertions))
    data = insertions.items()
    del insertions
    for k, v in data:
        line = k.split('|')
        line.insert(3, str(v))
        print('\t'.join(line), file=outfile)
    # Close files
    outfile.close()
    bamfile.close()
