#!/usr/local/bin/python3
# simple_query_mgrast.py
# Simple example to query MG-RAST API using the requests module.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 07 Feb 2017
# Updated on 07 Feb 2017

from __future__ import absolute_import, print_function
import requests
import json
import os
import sys

# Set base URL
API_BASE_URL = 'http://api.metagenomics.anl.gov/'

# Declare resource path and parameters
resource = 'project/mgp128'
my_params = {'verbosity': 'full'}

# Issue request
full_path = os.path.join(API_BASE_URL, resource)
response = requests.get(full_path, params=my_params)

# Check that status code is 200 = good
if response.status_code == 200:
    print('Response is good!')
else:
    print('There was an error with the request: status code = ',
          response.status_code,
          file=sys.stderr)
    sys.exit(1)

print('Request URL =', response.url)

input('Press enter to see response data keys\n')

# Load response as a json dictionary
data = response.json()

print('Keys:')
for idx, k in enumerate(data, start=1):
    print(idx, ': ', k, sep='')
print('\n')

print('Project name:', data['name'])

input('Press enter to see metagenomes\n')
print('There are', len(data['metagenomes']), 'metagenomes in this project.')
for idx, m in enumerate(data['metagenomes'], start=1):
    print('METAGENOME', idx)
    print('Name:', m['name'])
    print('Biome:', m['biome'])
    print('Sequence type:', m['sequence_type'])
    print('Sequence method:', m['sequencing_method'])
    print('Number of sequences:', m['sequences'])
    print('Number of basepairs:', m['basepairs'])
    print('\n==============================================\n')
    input('')

input('Press enter to see all response data from API\n')

# Print out results using json
print(json.dumps(data, indent=4, sort_keys=True))

print('Script complete.', file=sys.stderr)
