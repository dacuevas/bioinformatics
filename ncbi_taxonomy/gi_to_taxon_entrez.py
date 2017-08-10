#!/usr/local/bin/python3
# gi_to_taxon_entrez.py
# Collect taxonomy information for given GI numbers by querying Entrez db
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 09 Aug 2017
# Updated on 09 Aug 2017


from __future__ import print_function, absolute_import, division
import sys
import os
import time
import datetime
import argparse
from Bio import Entrez
import re


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
def timestamp():
    """
    Return time stamp.
    """
    t = time.time()
    fmt = '[%Y-%m-%d %H:%M:%S]'
    return datetime.datetime.fromtimestamp(t).strftime(fmt)


def print_status(msg, end='\n'):
    """
    Print status message.
    """
    print('{}    {}'.format(timestamp(), msg), file=sys.stderr, end=end)
    sys.stderr.flush()


def reprint(msg, fh=sys.stderr):
    """
    Print status message with carriage return
    """
    print('\033[1A\033[K{}\n\033[1B'.format(msg), end='', file=fh)
    fh.flush()


def exit_script(num=1):
    """
    Exit script.
    """
    sys.exit(num)


def get_summary(gi, log):
    """
    Get summary info from Entrez

    :param gi: GI number
    :type gi: str
    :param log: Log file handle
    :type log: File
    :return: Tax ID and Name
    :rtype: str, str
    """
    tid = None
    name = None
    try:
        count = 0
        handle = Entrez.esummary(db='nucleotide', id=gi)
        for rec in Entrez.read(handle):
            try:
                tid = rec['TaxId']
                count += 1
            except KeyError:
                log.write(gi + '\tTaxId not found\n')
                continue
            try:
                name = rec['Title']
            except KeyError:
                log.write(gi + '\tTitle not found\n')
    except:
        raise

    if count > 1:
        log.write(gi + '\tMore than one TaxId found: ' + str(count) + '\n')
    return tid, name


def get_taxonomy(tid, log):
    """
    Get taxonomy info from Entrez

    :param tid: Taxonomy ID number
    :type tid: str
    :param log: Log file handle
    :type log: File
    :return: Taxonomy path
    :rtype: list
    """
    tax = None
    try:
        count = 0
        handle = Entrez.efetch(db='taxonomy', id=tid)
        for rec in Entrez.read(handle):
            try:
                tax = rec['Lineage']
                count += 1
            except KeyError:
                log.write(tid + '\tLineage not found\n')
                continue
    except:
        raise

    if count > 1:
        log.write(tid + '\tMore than one Lineage found: ' + str(count))
    if tax is None:
        return tax
    else:
        return tax.split('; ')


###############################################################################
# ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description='')
parser.add_argument('infile', help='Input file containing GI numbers')
parser.add_argument('outdir', help='Output directory')
parser.add_argument('email', help='Email address for Entrez')
parser.add_argument('-s', '--skipfile', help='GI numbers to skip',
                    default=None)
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()

# Check that input file exists
if not os.path.isfile(args.infile):
    print(args.infile, 'does not exist', file=sys.stderr)
    parser.print_usage()
    exit_script()

# Check that output directory exists
if not os.path.isdir(args.outdir):
    print(args.outdir, 'does not exist', file=sys.stderr)
    parser.print_usage()
    exit_script()

# Global variables
gi_file = args.infile
out_dir = args.outdir
vbs = args.verbose
skip_file = args.skipfile
gi_regex = re.compile('gi\|(\d+)\|')
Entrez.email = args.email

###############################################################################
# LOAD INPUT FILE
###############################################################################
# Create log file
log = open(os.path.join(args.outdir, 'log.txt'), 'w', buffering=1)
log.write(timestamp() + ' Starting script\n')
gi_counts = {}  # Hold gi numbers and counts

# Build set of GIs to skip if supplied
gi_to_skip = set()
if skip_file:
    print_status('Loading skip file')
    with open(skip_file, 'r') as f:
        for l in f:
            gi_to_skip.add(l.rstrip('\n'))
    print_status('Finished skip file')

# Read in GI file
print_status('Loading input file\n')
with open(gi_file, 'r') as f:
    for li, l in enumerate(f, start=1):
        if vbs:
            reprint('Reading entry ' + str(li))

        match = gi_regex.search(l)
        if not match:
            continue

        gi = match.group(1)
        if gi in gi_to_skip:
            continue

        if gi in gi_counts:
            gi_counts[gi] += 1
        else:
            gi_counts[gi] = 1

print_status('Finished input file')
print_status('Loaded {} unique GIs'.format(len(gi_counts)))

###############################################################################
# BEGIN ENTREZ QUERIES
###############################################################################
# Open output file
out_file = os.path.join(out_dir, 'tax_info.txt')
print_status('Creating output file ' + out_file)
missing_data = set()  # Hold GI numbers that were unsuccessful

print_status('Begin Entrez queries')
with open(out_file, 'w', buffering=1) as f:
    # Header info
    f.write('gi\tcount\ttaxonomy\n')

    # Iterate through GIs
    numGi = len(gi_counts)
    for i, gi in enumerate(sorted(gi_counts), start=1):
        if vbs:
            reprint('Working on {} ({} out of {})'.format(gi, i, numGi))
        tid = None
        name = None
        tax = None
        # Make Entrez call
        try:
            # Get TaxId first
            # If there is no information retrieved, set as None
            tid, name = get_summary(gi, log)
            if tid is None:
                log.write(gi + '\tNo tax ID found\n')
                sys.stderr.write(gi + '\tNo tax ID found\n')
                sys.stderr.flush()
                missing_data.add(gi)
                continue

            elif tid == '':
                log.write(gi + '\tTax ID was blank\n')
                sys.stderr.write(gi + '\tTax ID was blank\n')
                sys.stderr.flush()
                missing_data.add(gi)
                continue

            # Get taxonomy info next
            if tid is not None:
                tax = get_taxonomy(tid, log)
                if tax is None:
                    log.write(tid + '\tNo tax info found\n')
                    sys.stderr.write(tid + '\tNo tax info found\n')
                    sys.stderr.flush()
                    missing_data.add(gi)
                    continue

        except Exception as e:
            log.write(gi + '\tUnexpected error occurred: ' + str(e) + '\n')
            sys.stderr.write(gi + '\tUnexpected error occurred: '
                             + str(e) + '\n')
            sys.stderr.flush()
            missing_data.add(gi)
            continue

        # Print out tax info
        # Check if we were able to retrieve anything
        if tid is None:
            name = ''
        if tax is None:
            tax = ['']

        # Remove 'cellular organisms' from list
        if tax[0] == 'cellular organisms':
            tax = tax[1:]
        tax_str = '\t'.join(tax)
        count = str(gi_counts[gi])
        f.write('\t'.join([gi, count, name, tax_str]) + '\n')

print_status('Entrez queries complete')

# Write out missing data
if len(missing_data) == 0:
    print_status('All GI numbers were found')
else:
    print_status(
        'A total of {} GI numbers were not found'.format(len(missing_data)))
    with open(os.path.join(out_dir, 'missing_gi.txt'), 'r') as f:
        for gi in sorted(missing_data):
            f.write(gi + '\n')

log.write(timestamp() + ' Script complete\n')
log.close()
print_status('Script complete!')
