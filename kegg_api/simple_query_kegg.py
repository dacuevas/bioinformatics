#!/usr/local/bin/python3
# simple_query_kegg.py
# Simple example to query KEGG API using the requests module.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 07 Mar 2017
# Updated on 07 Mar 2017

from __future__ import absolute_import, print_function
import requests
import os
import sys

# Set base URL
API_BASE_URL = 'http://rest.kegg.jp/'

# Declare resource path and parameters
resource = 'conv/eco/ncbi-geneid'

# Issue request
full_path = os.path.join(API_BASE_URL, resource)
response = requests.get(full_path)

# Check that status code is 200 = good
if response.status_code == 200:
    print('Response is good!')
else:
    print('There was an error with the request: status code = ',
          response.status_code,
          file=sys.stderr)
    sys.exit(1)

print('Request URL =', response.url)
print('Text:')
if response.text.strip('\n') == '':
    sys.exit('No text')

lines = response.text.split('\n')
print('Size =', len(lines))
print('First 10 lines:')
for n in range(10):
    print(lines[n])
    geneid, kegg = lines[n].split('\t')
    print('GeneID:', geneid)
    print('KEGG:', kegg)
