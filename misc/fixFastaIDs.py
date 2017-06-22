#!/usr/local/bin/python3
# fixFastaIDs.py
# Find duplicate sequence IDs and rename them
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com
# Created on 21 Jun 2017
# Updated on 21 Jun 2017

from __future__ import print_function, absolute_import, division
import sys
import os
import time
import datetime
import argparse


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


def exit_script(num=1):
    """
    Exit script.
    """
    sys.exit(num)


###############################################################################
# ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description='Find duplicate sequence IDs '
                                 'and rename them')
parser.add_argument('fasta', help='FASTA file')
parser.add_argument('output_name', help='Output filename')
parser.add_argument('output_dir', help='Output directory')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()
ff = args.fasta
on = args.output_name
od = args.output_dir
vbs = args.verbose

# Check FASTA file exists
if not os.path.isfile(ff):
    print_status('FASTA file "' + ff + '" does not exist.')

###############################################################################
# BEGIN PROCESSING FASTA FILE
###############################################################################
seq_ids = set()  # Contains sequence IDs already seen
new_seq_id = 0
new_seq_name = 'SEQ'

# Open new file to write to
fout = open(os.path.join(od, on), 'w')

seq_num = 0
# Iterate through FASTA file
if vbs:
    print_status('Processing FASTA file')
with open(ff, 'r') as f:
    for line in f:
        line = line.rstrip('\n')

        # Check if we're at the sequence ID
        if line.startswith('>'):
            seq_num += 1
            if vbs:
                print_status('Processed ' + str(seq_num) + ' sequences', '\r')

            sid = line[1:]

            # Check if sequence ID was already seen
            if sid in seq_ids:
                new_seq_id += 1
                sid = new_seq_name + str(new_seq_id)
            seq_ids.add(sid)
            fout.write('>' + sid + '\n')
        # At the sequence
        else:
            fout.write(line + '\n')

if vbs:
    print_status('Processed ' + str(seq_num) + ' sequences')
    print_status('Renamed ' + str(new_seq_id) + ' sequences')
fout.close()
print_status('Program complete!')
