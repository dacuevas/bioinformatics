#!/usr/local/bin/python3
# simple_query.py
# Simple example to query NCBI FTP.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 07 Feb 2017
# Updated on 07 Feb 2017

from __future__ import absolute_import, print_function
from ftplib import FTP
import os
import sys

# Set base URL
FTP_URL = 'ftp.ncbi.nlm.nih.gov'

# Declare resource path
resource = 'genomes/refseq/assembly_summary_refseq.txt'

# Conenct to FTP site
try:
    print('Connecting to', FTP_URL)
    ftp = FTP(FTP_URL)
    ftp.login()  # Username and password not required for NCBI
    print('Connected!')

    # Make request for file
    data = []  # data will hold all lines of file
    print('Making request for', resource)
    resp = ftp.retrlines('RETR ' + resource,
                         callback=lambda x: data.append(x))
    print('Request complete!')

    print('Data file contains', len(data), 'rows')
    print('Columns:')
    for idx, c in enumerate(data[1].split('\t'), start=1):
        print(idx, ': ', c, sep='')

    input()

    print('\nFirst row:')
    for idx, d in enumerate(data[2].split('\t'), start=1):
        print('Column', idx, ': ', d, sep='')

    # Close connection
    ftp.close()

except:
    sys.exit('There was an error during the FTP connection')

print('Script complete.', file=sys.stderr)
