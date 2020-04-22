#!/usr/bin/env python
# coding: utf-8

# This program takes as input an aligned BAM file of calling card
# sequencing data. It reads the XI and GS tags and outputs data into 
# qBED format, which is suitable for visualization.

import argparse
from collections import Counter
import os
import pickle
import pysam
import random
import string
import sqlite3
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--barcode", type=str, nargs='+', required=False)
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-db", "--database", action="store_true")
parser.add_argument("--cache", action="store_true")
parser.add_argument("--cache-size", type=int, required=False, default=1e7)
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

def writeCallingCardOutput(data):#, filename, mode):
    for k, v in data:
        line = k.split('|')
        line.insert(3, str(v))
        print('\t'.join(line))

if __name__ == "__main__":
    # Open input and output files
    bamfile = pysam.AlignmentFile(args.input, "rb")
    # Instantiate DB, if necessary
    if args.database:
        insertions = tempInsertionDB()
    else:
        insertions = Counter()
    # Instantiate caching arguments, if necessary
    if args.cache:
        cache_list = []
        insertion_counter = 0
        previous_insertion_string = ""
    barcodeTags = args.barcode
    
    for read in bamfile:
        insertion_string = read.get_tag("XI")

        # Is this a new insertion string?
        if args.cache and insertion_string != previous_insertion_string:
            # Have we hit the cache limit?
            if insertion_counter == args.cache_size:
                cache_list.append("tmp_" + randomAlphaNumericString() + ".pkl")
                with open(cache_list[-1], 'wb') as f:
                    pickle.dump(insertions, f, pickle.HIGHEST_PROTOCOL)
                insertion_counter = 0
                insertions = Counter()
            insertion_counter += 1
            previous_insertion_string = insertion_string

        if barcodeTags:
            try:
                barcode = getBarcodeString(read, barcodeTags)
                if insertion_string.endswith('-') or insertion_string.endswith('+'):
                    insertion_key = '|'.join([insertion_string, barcode])
                else:
                    insertion_key = '|'.join([insertion_string, read.get_tag("GS"), barcode])
            except KeyError:
                continue
        else:
            if insertion_string.endswith('-') or insertion_string.endswith('+'):
                insertion_key = insertion_string
            else:
                insertion_key = '|'.join([insertion_string, read.get_tag("GS")])
        
        insertions.update([insertion_key])

    if args.cache:
        for _ in cache_list:
            with open(_, 'rb') as f:
                data = pickle.load(f).items()
            writeCallingCardOutput(data)
    data = insertions.items()
    del insertions
    writeCallingCardOutput(data)

    # Clean up
    bamfile.close()
    for _ in cache_list:
        os.remove(_)