#!/usr/local/bin/python3
# fixFastaIDs.py
# Find duplicate sequence IDs and rename them
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com
# Created on 21 Jun 2017
# Updated on 22 Jun 2017

from __future__ import print_function, absolute_import, division
import sys
import os
import time
import datetime
import argparse
import locale


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


def format_num(num):
    """
    Format the number with commas.

    :param num: Number to format
    :type num: int or float
    :return: Formatted number
    :rtype: str
    """
    if isinstance(num, int):
        return locale.format('%d', num, grouping=True)
    elif isinstance(num, float):
        return locale.format('%.1f', num, grouping=True)
    else:
        return None


def process_rate(num_seqs, stime):
    """
    Calculate amount of time to process sequences

    :param num_seqs: Number of sequences processed
    :type num_seqs: int
    :param stime: Start time in seconds (since the epoch)
    :type stime: float
    :return: Number of sequences processed per second
    :rtype: float
    """
    etime = time.time()
    return num_seqs / (etime - stime)


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

# Setup locale for printing out large numbers
locale.setlocale(locale.LC_ALL, 'en_US')

###############################################################################
# BEGIN PROCESSING FASTA FILE
###############################################################################
seq_ids = set()  # Contains sequence IDs already seen
new_seq_id = 0
new_seq_name = 'SEQ'

# Open new file to write to
fout = open(os.path.join(od, on), 'w')

# Verbose output variables
seq_num = 0
start_time = time.time()
nseq_for_timing = 10000
rate = 0

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
                # Calculate rate of process (# sequences / second)
                if seq_num % nseq_for_timing == 0:
                    rate = process_rate(nseq_for_timing, start_time)
                    rate = format_num(rate)
                    start_time = time.time()

                sn = format_num(seq_num)
                ps = 'Processed {} sequences [~ {} seqs/sec]'.format(sn, rate)
                ps = '{:{}}'.format(ps, 15)
                print_status(ps, '\r')

            sid = line[1:]
            # Check if sequence ID was already seen
            if sid in seq_ids:
                # Rename sequence ID
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
