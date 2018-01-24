#!/usr/local/bin/python3
# mgrast-get-function-organism.py
# Read in file of MG-RAST IDs and obtain function and organism information
# from the MG-RAST API "Annotation/similarity" resource
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 23 Jan 2018
# Updated on 23 Jan 2018

from __future__ import absolute_import, print_function
import requests
import os
import sys
import argparse
import re
import time
import datetime


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
def timestamp():
    """
    Return time stamp.
    """
    t = time.time()
    fmt = "[%Y-%m-%d %H:%M:%S]"
    return datetime.datetime.fromtimestamp(t).strftime(fmt)


def print_status(msg, end='\n'):
    """
    Print status message.
    """
    print("{}    {}".format(timestamp(), msg), file=sys.stderr, end=end)
    sys.stderr.flush()


###############################################################################
# ARGUMENT PARSING
###############################################################################
# Parse command line arguments
desc = "Read in file of MG-RAST IDs and obtain function and organism " \
       "information from the MG-RAST API 'Annotation/similarity resource"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("mgfile", help="Metagenome ID file")
parser.add_argument("outfile", help="Output file location")
parser.add_argument("--min_identity", help="Minimum identity threshold",
                    type=int)
parser.add_argument("--min_aln_length",
                    help="Minimum alignment length threshold", type=int)
parser.add_argument("--auth", help="MG-RAST authentication key")
parser.add_argument("--source", help="Database source to use",
                    choices=["RefSeq", "GenBank", "SEED", "PATRIC", "KEGG",
                             "SwissProt"], default="SEED")
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()

# Check if metagenome file exists
if not os.path.isfile(args.mgfile):
    print("ERROR: Could not find meteagenome ID file:", args.mgfile,
          file=sys.stderr)
    sys.exit()

auth = args.auth if args.auth else None
minID = args.min_identity if args.min_identity else None
minLen = args.min_aln_length if args.min_aln_length else None
source = args.source

# Set base URL
API_BASE_URL = "http://api.metagenomics.anl.gov/annotation/similarity"
# Declare resource path and set parameters
myParams = {"type": "all", "source": args.source}
if minID:
    myParams["identity"] = minID
if auth:
    myParams["auth"] = auth
if minLen:
    myParams["length"] = minLen

# Print out options
if args.verbose:
    print("***********************************************", file=sys.stderr)
    print("Running mgrast-get-function-organism.py", file=sys.stderr)
    print("Making queries from", API_BASE_URL, file=sys.stderr)
    print("Options supplied:", file=sys.stderr)
    print("     Input file:                        ", args.mgfile, file=sys.stderr)
    print("     Output file:                       ", args.outfile, file=sys.stderr)
    print("     Source database:                   ", source, file=sys.stderr)
    print("     Minimum identity threshold:        ", minID, file=sys.stderr)
    print("     Minimum alignment length threshold:", minLen, file=sys.stderr)
    print("     Authentication key:                ", auth, file=sys.stderr)
    print("***********************************************\n", file=sys.stderr)
    sys.stderr.flush()

# Open file for output
fout = open(args.outfile, "w")

# Loop through file
with open(args.mgfile, "r") as mgf:
    for l in mgf:
        mgID = l.rstrip()
        if args.verbose:
            print_status("Retrieving metagenome ID: " + mgID)

        # Issue request
        full_url = os.path.join(API_BASE_URL, mgID)
        response = requests.get(full_url, params=myParams)

        # Check that status code is 200 = good
        if response.status_code != 200:
            print_status("ERROR: Connection status code: " +
                         str(response.status_code))
            print_status("Skipping " + mgID)
            continue
        elif args.verbose:
            print_status("Data parsing for metagenome ID:" + mgID)

        # Load response as a json dictionary
        # data = response.json()
        # Iterate through tab-delimited data
        for lineNum, line in enumerate(response.text.split("\n")):
            # First line is the header
            if lineNum == 0:
                continue

            # Download complete terminates data
            if "Download complete" in line:
                break

            if args.verbose:
                print_status("Processing line " + str(lineNum), end="\r")
            contents = line.rstrip("\n").split("\t")

            # Extract information
            queryID = str(contents[0])
            percID = str(contents[2])
            alen = str(contents[3])
            mismatches = str(contents[4])
            gaps = str(contents[5])
            qstart = str(contents[6])
            qend = str(contents[7])
            eval = str(contents[10])

            # Find mgmID
            if re.match(r"(mgm\d+\.\d)", queryID):
                mgmID = re.match(r"(mgm\d+\.\d)", queryID).group(1)
            else:
                mgmID = mgID

            # Loop through annotation data
            for info in contents[-1].split(";"):
                acc = func = org = "None"
                if re.search(r"accession=\[(.*)\]", info):
                    acc = re.search(r"accession=\[(.*?)\]", info).group(1)
                if re.search(r"function=\[(.*)\]", info):
                    func = re.search(r"function=\[(.*?)\]", info).group(1)
                if re.search(r"organism=\[(.*)\]", info):
                    org = re.search(r"organism=\[(.*?)\]", info).group(1)

                toPrint = "\t".join([mgmID, queryID, percID, alen, mismatches,
                                     gaps, qstart, qend, eval, acc, func, org])
                fout.write(toPrint + "\n")
        if args.verbose:
            print("", file=sys.stderr)
            print_status(mgID + " parsing complete")

fout.close()
print_status("Script complete.")
