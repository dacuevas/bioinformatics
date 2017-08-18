from __future__ import print_function, absolute_import, division
import sys
import os
import time
import datetime
import argparse
from Bio import SeqIO


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


def check_file(file_path, directory_path=None):
    """
    Check if file exists.
    """
    if directory_path:
        file_path = os.path.join(directory_path, file_path)
    return os.path.isfile(file_path)
    
    
def exit_script(num=1):
    """
    Exit script.
    """
    sys.exit(num)

###############################################################################
# ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description='Convert SFF file to '
                                             'FASTA or FASTQ')
parser.add_argument('sff', help='Input SFF file path')
parser.add_argument('out_file', help='Output file path '
                                     '(without format extension) ')
parser.add_argument('format', help='Output format',
                    choices=['fasta', 'fastq'], default='fasta')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()

# Check to see if input file exists
if not check_file(args.sff):
    parser.print_help()
    exit_script()

###############################################################################
# RUN CONVERSION
###############################################################################
if args.verbose:
    print_status('Converting ' + args.sff + ' into a ' + args.format + ' file')
SeqIO.convert(args.sff, 'sff', args.out_file + '.' + args.format, args.format)
print_status('Script complete!')
