# add_rast_anno.py
# Add RAST annotation information to sequence headers, output to stdout
#
# Author Daniel A. Cuevas (dcuevas08@gmail.com)
# Created on 13 Feb 2018
# Updated on 13 Feb 2018


from __future__ import print_function, absolute_import, division
import sys
import os
import argparse


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
def check_file(file_path, directory_path=None):
    """
    Check if file exists.
    """
    if directory_path:
        file_path = os.path.join(directory_path, file_path)
    return os.path.isfile(file_path)


###############################################################################
# ARGUMENT PARSING
###############################################################################
desc = "Add RAST annotation information to sequence headers, output to stdout"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("fasta", help="FASTA file")
parser.add_argument("rast", help="RAST protein annotation file")
parser.add_argument("--no_header", action="store_true",
                    help="Flag if RAST file has not header")
parser.add_argument("--skip_seqs", action="store_true",
                    help="Flag to not print out sequences not found in RAST")

args = parser.parse_args()

# Check if files exists
if not check_file(args.fasta):
    sys.exit("FASTA file (" + args.fasta + ") does not exist")
if not check_file(args.rast):
    sys.exit("RAST file (" + args.rast + ") does not exist")

# Read in annotation file and make PEG -> annotation dictionary
rastData = {}
with open(args.rast, "r") as f:
    if not args.no_header:
        # Skip header
        header = f.readline()
    for l in f:
        data = l.rstrip("\n").split("\t")
        pegID = data[0]
        function = data[8]
        if pegID in rastData:
            print("Warning:", pegID, "already seen. Skipping.",
                  file=sys.stderr)
            continue
        rastData[pegID] = function

# Iterate through FASTA file
# Keep track of number of sequences and number of sequences not found
numNotFound = 0
numSequences = 0
skip = False  # Flag to skip sequences if the ID was not found in the RAST data
with open(args.fasta, "r") as f:
    for l in f:
        l = l.rstrip("\n")
        if l.startswith(">"):
            skip = False
            numSequences += 1

            # Extract the sequence ID
            seqID = l[1:]

            # Check if sequence ID exists in RAST file
            if seqID in rastData:
                print(">{}_{}".format(seqID, rastData[seqID]))

            # Print out warning message if not found
            else:
                print("Warning: sequence", seqID, "not found in RAST file",
                      file=sys.stderr)
                numNotFound += 1

                # Check if we want to print out the sequence or not
                if args.skip_seqs:
                    skip = True
                else:
                    print(l)
        elif not skip:
            print(l.upper())

if numNotFound > 0:
    print("There were", numNotFound, "out of", numSequences,
          "sequences not found in RAST file", file=sys.stderr)

print("Script complete", file=sys.stderr)
