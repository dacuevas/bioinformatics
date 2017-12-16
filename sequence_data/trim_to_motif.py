# trim_to_motif.py
# Trim off sequences up to a given motif
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Date modified: 16 December 2017


from __future__ import print_function, absolute_import, division
import sys
import os
import argparse


# Function definitions
def check_file(file_path, directory_path=None):
    """
    Check if file exists.
    """
    if directory_path:
        file_path = os.path.join(directory_path, file_path)
    return os.path.isfile(file_path)


def check_for_motif(motifs, seq):
    """
    Check if motif exists.
    """
    for motif in motifs:
        idx = seq.find(motif)
        if idx > -1:
            return True, seq[idx:] + seq[:idx]
    # Did not find any motifs
    return False, seq


# Argument parsing
desc = "Trim off sequences up to a given motif"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("seq_file", help="Sequence input file")
parser.add_argument("motif_file", help="Motif to find")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Verbose output")

args = parser.parse_args()
if not check_file(args.seq_file):
    sys.exit("Sequence file not found: " + args.seq_file)
if not check_file(args.motif_file):
    sys.exit("Motif file not found: " + args.motif_file)

# Read in motif file
motifs = set()
with open(args.motif_file, "r") as f:
    for l in f:
        motif = l.rstrip("\n").upper()
        if motif == "":
            sys.exit("Sequence motif is empty")
        for ci, c in enumerate(motif, start=1):
            if c != "A" and c != "T" and c != "G" and c != "C":
                sys.exit("Nucleotide " + str(ci) + "(" + c + ") in motif "
                                                             "is not "
                                                             "A, T, C, or G")
        motifs.add(motif)

# Process file
no_motif = {}
with_motif = {}
with open(args.seq_file, "r") as f:
    seqID = ""
    seq = ""
    seqCounter = 0
    numFound = 0
    numNotFound = 0
    for li, l in enumerate(f, start=1):
        l = l.rstrip("\n")

        # At sequence ID
        if l.startswith(">"):
            if li > 1:
                # At a new sequence, process previous one
                found, trimseq = check_for_motif(motifs, seq)
                if found:
                    with_motif[seqID] = trimseq
                    numFound += 1
                else:
                    no_motif[seqID] = seq
                    numNotFound += 1

            # Reinitialize info
            seq = ""
            seqID = l[1:]
            seqCounter += 1
            print("Processing sequence", seqCounter, end="\r", file=sys.stderr)
            continue

        # At sequence
        seq += l.upper()

print("", file=sys.stderr)
# Process final sequence
found, trimseq = check_for_motif(motifs, seq)
if found:
    with_motif[seqID] = trimseq
    numFound += 1
else:
    no_motif[seqID] = seq
    numNotFound += 1
print("Found motif in", numFound, "out of", seqCounter, "sequences",
      file=sys.stderr)


# Print out to new files
if len(with_motif) != 0:
    with open(args.seq_file + ".trimmed.fna", "w") as fw:
        for seqID, seq in with_motif.items():
            fw.write(">" + seqID + "\n" + seq + "\n")
else:
    print("No sequences with motif were found!", file=sys.stderr)

if len(no_motif) != 0:
    with open(args.seq_file + ".untrimmed.fna", "w") as fw:
        for seqID, seq in no_motif.items():
            fw.write(">" + seqID + "\n" + seq + "\n")
else:
    print("All sequences found a motif!", file=sys.stderr)
print("Script complete!", file=sys.stderr)