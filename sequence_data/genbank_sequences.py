#!/usr/local/bin/python3
# genbank_sequences.py
# Extract sequences from GenBank file
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 31 Aug 2017
# Updated on 31 Aug 2017

from __future__ import absolute_import, print_function
import os
import sys
import argparse
from Bio import SeqIO

desc = "Extract sequences from a GenBank file"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("gb_file", help="GenBank file location")
parser.add_argument("type", help="Type of sequences to extract",
                    choices=["nt", "aa"])
parser.add_argument("out_file", help="Output file path")
parser.add_argument("-i", "--id", help="Additional information to add to "
                    "sequence identifier")
args = parser.parse_args()

# Check that files exist
if not os.path.isfile(args.gb_file):
    sys.exit("GenBank file '" + args.gb_file + "' does not exist.")

info = args.id if args.id else ""

records = SeqIO.parse(args.gb_file, "genbank")
num_recs = len(list(records))
print("{} records in {}".format(num_recs, args.gb_file), file=sys.stderr)

# Re read records
records = SeqIO.parse(args.gb_file, "genbank")
# Iterate through each GenBank record
with open(args.out_file, "w") as fo:
    for idx, rec in enumerate(records, start=1):
        print("Processing record {} out of {}".format(idx, num_recs),
              end="\r", file=sys.stderr)
        # Iterate through each feature in the record
        for feature in rec.features:
            if "product" not in feature.qualifiers:
                continue
            prod = feature.qualifiers["product"][0]

            if args.type == "aa":
                seq = str(feature.extract(rec.seq).translate())
                type = "aa"
            else:
                seq = str(feature.extract(rec.seq))
                type = "bp"

            # Print out DNA sequence
            if info:
                print(">", prod," [", info,  "]", "\n",
                      seq, sep="", file=fo)
            else:
                print(">", prod," [", len(seq), type,  "]", "\n",
                      seq, sep="", file=fo)
print("", file=sys.stderr)

print("Script complete!", file=sys.stderr)