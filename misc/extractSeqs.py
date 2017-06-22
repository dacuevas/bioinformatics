#!/usr/local/bin/python3
# extractSeqs.py
# Pull out FASTA sequences from BLAST output format 6
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com
# Created on 14 Jun 2017
# Updated on 21 Jun 2017

from __future__ import print_function, absolute_import, division
import sys
import os
import time
import datetime
import argparse
import re
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


def parse_blast(line, line_num):
    """
    Parse a line from BLAST output in format 6.

    :param line: input line
    :type line: str
    :param line_num: line number
    :type line_num: int
    :return: a mapping BLAST column to value
    :rtype: dict
    """
    contents = line.rstrip('\n').split('\t')
    if len(contents) < 12:
        print_status('Less than 12 columns in line ' + str(line_num))
        exit_script()

    blast_map = dict()
    blast_map['qseqid'] = contents[0]
    blast_map['sseqid'] = contents[1]
    blast_map['pident'] = contents[2]
    blast_map['length'] = contents[3]
    blast_map['mismatch'] = contents[4]
    blast_map['gapopen'] = contents[5]
    blast_map['qstart'] = contents[6]
    blast_map['qend'] = contents[7]
    blast_map['sstart'] = contents[8]
    blast_map['send'] = contents[9]
    blast_map['evalue'] = contents[10]
    blast_map['bitscore'] = contents[11]
    return blast_map


def get_seq_ids(bf):
    """
    Read in BLAST output file and get query sequence IDs.

    :param bf: BLAST output filepath
    :type bf: str
    :return: query IDs taken from BLAST output file
    :rtype: set
    """
    seq_ids = set()
    if vbs:
        print_status('Parsing BLAST output file')
    with open(bf, 'r') as f:
        for ln, line in enumerate(f, start=1):
            blast = parse_blast(line, ln)
            # Add query ID to set
            seq_ids.add(blast['qseqid'])
    if vbs:
        print_status('Obtained ' + format_num(len(seq_ids)) + ' unique IDs')
    return seq_ids


def print_seqs(ff, seq_ids, on, od):
    """
    Read through FASTA sequences and print out specific sequences.

    :param ff: FASTA filepath
    :type ff: str
    :param seq_ids: sequence IDs to print out
    :type seq_ids: set
    :param on: Output filename
    :type on: str
    :param od: Output directory
    :type od: str
    :return: None
    """
    # Open output file handle
    fout = open(os.path.join(od, on), 'w')

    print_flag = False  # Flag to print sequence
    print_seq = ''  # Holds sequence
    if vbs:
        print_status('Filtering through FASTA file')
    # Verbose output variables
    seq_num = 0
    start_time = time.time()
    nseq_for_timing = 10000
    rate = 0
    num_filtered = 0
    with open(ff, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            # Check if we are looking at a sequence ID
            if re.match(r'>', line):
                seq_num += 1
                if vbs:
                    # Calculate rate of process (# sequences / second)
                    if seq_num % nseq_for_timing == 0:
                        rate = process_rate(nseq_for_timing, start_time)
                        rate = format_num(rate)
                        start_time = time.time()

                    sn = format_num(seq_num)
                    ps = 'Processed {} sequences [~ {} seqs/sec]'.format(sn,
                                                                         rate)
                    ps = '{:{}}'.format(ps, 15)
                    print_status(ps, '\r')

                # When we reach the sequence ID we have to check if there is a
                # sequence that needs to be printed
                if print_flag:
                    fout.write(print_seq + '\n')
                    # Reset variables
                    print_flag = False
                    print_seq = ''

                # Extract sequence ID
                qid = line[1:]

                # Check if sequence ID is in the set
                if qid in seq_ids:
                    fout.write(line + '\n')
                    print_flag = True
                    num_filtered += 1

            # Check if we need to print this sequence
            elif print_flag:
                print_seq += line

    # At the end of the file we have to check if there is a sequence that
    # still needs to be printed
    if print_flag:
        fout.write(print_seq + '\n')
    if vbs:
        print_status('Processed ' + str(seq_num) + ' sequences')
        print_status('Resulted in ' + str(num_filtered) + ' sequences')
    fout.close()

###############################################################################
# ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description='Extract FASTA sequences based '
                                 'on BLAST output')
parser.add_argument('fasta', help='FASTA filepath')
parser.add_argument('blast', help='BLAST output filepath')
parser.add_argument('output_name', help='Output filename')
parser.add_argument('output_dir', help='Output directory')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()
ff = args.fasta
bf = args.blast
on = args.output_name
od = args.output_dir
vbs = args.verbose

# Check FASTA file exists
if not os.path.isfile(ff):
    print_status('FASTA file "' + ff + '" does not exist.')
    exit_script()

# Check blast file exists
if not os.path.isfile(bf):
    print_status('BLAST file "' + bf + '" does not exist.')
    exit_script()

# Setup locale for printing out large numbers
locale.setlocale(locale.LC_ALL, 'en_US')

###############################################################################
# BEGIN PROCESSING FILES
###############################################################################
# Parse BLAST file
seq_ids = get_seq_ids(bf)

# Print out sequences
print_seqs(ff, seq_ids, on, od)

print_status('Program complete!')
