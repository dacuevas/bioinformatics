#!/usr/local/bin/python3
# extractSeqs.py
# Pull out FASTA sequences from BLAST output format 6
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com
# Created on 14 Jun 2017
# Updated on 14 Jun 2017

from __future__ import absolute_import, print_function, division
import sys
import os
import time
import datetime
import argparse


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
def timestamp():
    """Return time stamp"""
    t = time.time()
    fmt = '[%Y-%m-%d %H:%M:%S]'
    return datetime.datetime.fromtimestamp(t).strftime(fmt)


def print_status(msg, end='\n'):
    """Print status message"""
    print('{}    {}'.format(timestamp(), msg), file=sys.stderr, end=end)
    sys.stderr.flush()


def exit_script(num=1):
    """Exit script"""
    sys.exit(num)


###############################################################################
# ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description='Extract FASTA sequences based '
                                 'on BLAST output')
parser.add_argument('fasta', help='Profiling results file')
parser.add_argument('blast', help='Binning results file')
parser.add_argument('output_name', help='Output name')
parser.add_argument('output_dir', help='Output directory')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()
ff = args.fasta
bf = args.blast
on = args.output_name
od = args.output_dir

# Check fasta file exists
if not os.path.isfile(ff):
    print_status('FASTA file "' + ff + '" does not exist.')
    exit_script()

# Check blast file exists
if not os.path.isfile(bf):
    print_status('BLAST file "' + bf + '" does not exist.')
    exit_script()

###############################################################################
# BEGIN PROCESSING FILE
###############################################################################
