#!/usr/local/bin/python3
# reactions_to_kegg.py
# Use PyFBA reaction IDs to obtain KEGG KO IDs and their associated pathways.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 19 Jul 2017
# Updated on 19 Jul 2017

import requests
import argparse
import PyFBA
import os
import sys
import re


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
def parse_response(response):
    """
    Parse the KEGG API response text. Text is all plain text without an
    easily digestible format. Feature headers are listed at the beginning of
    each line. For the GET reaction response, the headers are:
    1) ENTRY
    2) NAME
    3) DEFINITION
    4) EQUATION
    5) COMMENT
    6) RCLASS
    7) ENZYME
    8) PATHWAY
    9) ORTHOLOGY

    :param response: API response text
    :type response: str
    :return: None
    """
    # Response contains newlines
    response = response.split('\n')

    in_orthology = False  # Flag to indicate to obtain all KO values
    for line in response:
        if re.match(r'NAME', line):
            parse_name_line(line)

        elif re.match(r'ENZYME', line):
            parse_enzyme_line(line)

        elif re.match(r'ORTHOLOGY', line) or in_orthology:
            in_orthology = True
            parse_orthology_line(line)


def parse_name_line(line):
    """

    :param line:
    :return:
    """
    pass


def parse_enzyme_line(line):
    """

    :param line:
    :return:
    """
    pass


def parse_orthology_line(line):
    """

    :param line:
    :return:
    """
    pass


###############################################################################
# ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description='Use PyFBA reaction IDs to '
                                 'obtain KEGG KO IDs and their '
                                 'associated pathways')
parser.add_argument('mseed_to_kegg', help='ModelSEED reaction mapper file')
parser.add_argument('model_name', help='Model name')
parser.add_argument('model_dir', help='Model directory')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()

# Check that the mapper file exists
if not os.path.isfile(args.mseed_to_kegg):
    sys.exit('ModelSEED to KEGG mapper file does not exist!')
# Check that model directory exists
if not os.path.isdir(args.model_dir):
    sys.exit('Model directory does not exist!')


###############################################################################
# PROCESS MODEL REACTIONS
###############################################################################
# Load ModelSEED reaction ID to KEGG reacton ID mapper
mseed_to_kegg = {}
with open(args.mseed_to_kegg, 'r') as f:
    # Read header line
    header = f.readline()
    for l in f:
        l = l.rstrip('\n')
        contents = l.split('\t')
        mseed, kegg = contents[:2]

        # Check if reaction ID has already been encountered
        if mseed in mseed_to_kegg:
            print(mseed + ' already encountered in mapper file',
                  file=sys.stderr)
        else:
            mseed_to_kegg[mseed] = kegg

print(str(len(mseed_to_kegg)) + ' reacton IDs found in mapper file',
      file=sys.stderr)

# Load model
model = PyFBA.model.load_model(args.model_dir, args.model_name)

# Set KEGG API URL
API_BASE_URL = 'http://rest.kegg.jp/'

# Output header
print('mseed_id\tkegg_id\tec\tpathway')

for mseed_rxn in model.reactions:
    if mseed_rxn in mseed_to_kegg:
        kegg_rxn = mseed_to_kegg[mseed_rxn]
    else:
        print(mseed_rxn + ' not found in mapper. Skipping.', file=sys.stderr)
        continue

    # Set resource path
    resource = 'get/' + kegg_rxn

    # Issue request
    full_path = os.path.join(API_BASE_URL, resource)
    response = requests.get(full_path)

    # Check that status code is 200 = good
    if response.status_code != 200:
        print('There was an error with the request: status code = ',
              response.status_code,
              file=sys.stderr)
        print('ModelSEED id: ' + mseed_rxn, file=sys.stderr)
        print('KEGG id: ' + kegg_rxn, file=sys.stderr)
        continue

    if response.text.strip('\n') == '':
        print('Response was empty for ' + kegg_rxn, file=sys.stderr)
        continue

    parse_response()
