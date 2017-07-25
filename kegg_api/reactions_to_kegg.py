#!/usr/local/bin/python3
# reactions_to_kegg.py
# Use PyFBA reaction IDs to obtain KEGG KO IDs and their associated pathways.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 19 Jul 2017
# Updated on 23 Jul 2017

from __future__ import print_function, absolute_import, division
import requests
import argparse
import PyFBA
import os
import sys
import re


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
def parse_response(response, mseed, kegg):
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
    :param mseed: Model SEED reaction ID
    :type mseed: str
    :param kegg: KEGG reaction ID
    :type kegg: str
    :return: None
    """
    # Check if response contains all items
    if not re.search(r'NAME', response):
        print('No NAME info for', mseed, ':', kegg, file=sys.stderr)
        no_name = True
    else:
        no_name = False

    if not re.search(r'ENZYME', response):
        print('No ENZYME info for ', mseed, ':', kegg, file=sys.stderr)
        no_enz = True
    else:
        no_enz = False

    if not re.search(r'ORTHOLOGY', response):
        print('No ORTHOLOGY info for ', mseed, ':', kegg, file=sys.stderr)
        no_orth = True
    else:
        no_orth = False

    # Check if there is no information in the response
    if no_name and no_enz and no_orth:
        print('\t\t\t')
        return

    # Response contains newlines
    response = response.split('\n')

    all_lines = []
    in_orthology = False
    for line in response:
        if re.match(r'(///)|(REFERENCES)', line):
            # Reached the end of the orthology lines
            parse_orthology_lines(all_lines)
            break

        elif re.match(r'NAME', line):
            parse_name_line(line)
            if no_enz and no_orth:
                print('\t')
                return

        elif re.match(r'ENZYME', line):
            # Check if there was name info missing
            if no_name:
                print('\t', end='')
            parse_enzyme_line(line)

        elif re.match(r'ORTHOLOGY', line):
            # Check if there was name and/or enzyme info missing
            if no_name and no_enz:
                print('\t\t')
            elif no_enz:
                print('\t')
            in_orthology = True
            ko = line.split()[1]
            all_lines = [ko]

        elif in_orthology:
            ko = line.split()[0]
            all_lines.append(ko)


def parse_name_line(line):
    """
    Parse the NAME line from the KEGG GET response

    :param line: Name line from the response
    :type line: str
    :return: None
    """
    name = re.split(r'\s{2,}', line)[1]
    print(name + '\t', end='')


def parse_enzyme_line(line):
    """
    Parse the ENZYME line from the KEGG GET response

    :param line: ENZYME line from the response
    :type line: str
    :return: None
    """
    enzymes = line.split()
    print(';'.join(enzymes[1:]) + '\t', end='')


def parse_orthology_lines(lines):
    """
    Parse the ORTHOLOGY (KO) lines from the KEGG GET response

    :param lines: ORTHOLOGY lines from the response
    :type lines: list
    :return: None
    """
    print(';'.join(lines))


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

print(str(len(mseed_to_kegg)), 'reacton IDs found in mapper file',
      file=sys.stderr)

# Load model
model = PyFBA.model.load_model(args.model_dir, args.model_name)
print('Model contains', str(model.number_of_reactions()), 'reactions',
      file=sys.stderr)

# Set KEGG API URL
API_BASE_URL = 'http://rest.kegg.jp/'

# Output header
print('mseed_id', 'equation', 'kegg_id', 'name', 'ec','pathway', sep='\t')

# Iterate through reaction list from model
for i, mseed_rxn in enumerate(model.reactions, start=1):
    print('Processing reaction', str(i), 'of',
          str(model.number_of_reactions()),
          end='\r', file=sys.stderr)
    if mseed_rxn in mseed_to_kegg:
        kegg_rxn = mseed_to_kegg[mseed_rxn]
    else:
        print(mseed_rxn, 'not found in mapper. Skipping.', file=sys.stderr)
        continue

    # Set resource path
    resource = 'get/' + kegg_rxn

    # Issue request
    full_path = os.path.join(API_BASE_URL, resource)
    response = requests.get(full_path)

    # Check that status code is 200 = good
    if response.status_code != 200:
        print('There was an error with the request: status code =',
              response.status_code,
              file=sys.stderr)
        print('ModelSEED id:', mseed_rxn, file=sys.stderr)
        print('KEGG id:', kegg_rxn, file=sys.stderr)
        print(mseed_rxn, model.reactions[mseed_rxn].equation,
              kegg_rxn, 'None', 'None', 'None', sep='\t')
        continue

    if response.text.strip('\n') == '':
        print('Response was empty for', kegg_rxn, file=sys.stderr)
        continue

    # Print reaction ID
    print(mseed_rxn, model.reactions[mseed_rxn].equation, kegg_rxn, '',
          end='', sep='\t')
    parse_response(response.text, mseed_rxn, kegg_rxn)

print('\nScript completed!', file=sys.stderr)

