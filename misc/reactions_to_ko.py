#!/usr/local/bin/python3
# reactions_to_ko.py
# Use PyFBA reaction IDs to obtain KEGG KO IDs and their associated pathways.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 02 Aug 2017
# Updated on 02 Aug 2017
from __future__ import print_function, absolute_import, division
import sys
import os
import time
import datetime
import argparse
import PyFBA
import requests
import re


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


def make_kegg_query(ms, ko, log):
    # Set resource path
    resource = 'get/' + ko

    # Issue request
    full_path = os.path.join(KEGG_BASE_URL, resource)
    response = requests.get(full_path)

    # Check that status code is 200 = good
    if response.status_code != 200:
        log.write('There was an error with the request '
                  'regarding ModelSEED reaction ' + ms
                  + ' and KO ' + ko + ', status code = '
                  + response.status_code + '\n')
        continue

    # Check that response is not empty
    elif response.text.strip('\n') == '':
        log.write('Response was empty for KO '
                  + ko + '\n')
        continue

    # Parse response
    parse_kegg_response(response.text, mseed_rxn, ko, log_out)


def parse_kegg_response(res, ms, ko, log):
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

    :param res: API response text
    :type res: str
    :param ms: Model SEED reaction ID
    :type ms: str
    :param ko: KEGG KO ID
    :type ko: str
    :param log: Log output file handle
    :type log: File
    :return: None
    """
    # Check if response contains KO info
    if not re.search(r'ORTHOLOGY', res):
        log.write('No ORTHOLOGY info for ' + ms + ':' + ko + '\n')
        return False

    # Response contains newlines
    res = res.rstrip().split('\n')

    curr_section = ''
    for line in res:
        # First 12 characters contain SECTION
        # Rest of characters are DATA
        kegg_section = line[:12].strip()
        kegg_data = line[12:].strip()
        if kegg_section != '':
            curr_section = kegg_section

        elif curr_section == 'ORTHOLOGY':
            ko_id = re.match(r'K\d+', kegg_data).group(0)
            print(ko_id)


###############################################################################
# ARGUMENT PARSING
###############################################################################
parser = argparse.ArgumentParser(description='Obtain EC and KO values for '
                                 'ModelSEED reactions')
parser.add_argument('model', help='Model name')
parser.add_argument('modeldir', help='Model directory')
parser.add_argument('mseeddir', help='ModelSEED directory')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()

# Check that model directory exists
if not os.path.isdir(args.modeldir):
    print_status('Model directory does not exist')
    exit_script()

# Check that ModelSEED directory exists
if not os.path.isdir(args.mseeddir):
    print_status('ModelSEED directory does not exist')
    exit_script()

# Check that ModelSEED directory is new
alias_file = os.path.join(args.mseeddir,
                          'Biochemistry/Aliases/Reactions_Aliases.tsv')
if not os.path.isfile(alias_file):
    print_status('ModelSEED Reactions_Aliases.tsv file does not exist')
    exit_script()

# Set base URLs
KEGG_BASE_URL = 'http://rest.kegg.jp/'
BIGG_BASE_URL = 'http://bigg.ucsd.edu/api/v2/'

###############################################################################
# BEGIN PROCESSING
###############################################################################
# aliases will contain Reactions_Aliases.tsv info
aliases = {}
# Load ModelSEED alias file
if args.verbose:
    print_status('Parsing alias file')
with open(alias_file, 'r') as f:
    header = f.readline()  # Skipping header
    for l in f:
        ms, oldms, alias, db = l.rstrip('\n').split('\t')

        if ms.strip() == '':
            continue

        # Find which database it belongs to
        if db.startswith('KEGG'):
            db = 'KEGG'

        # ModelSEED reaction column can have multiple
        for ms_curr in ms.split('|'):
            if ms_curr not in aliases:
                aliases[ms_curr] = {db: set()}
            elif db not in aliases[ms_curr]:
                aliases[ms_curr][db] = set()

            aliases[ms_curr][db].add(alias)

if args.verbose:
    print_status('Alias file parsing complete')

# Load model
if args.verbose:
    print_status('Loading model')

model = PyFBA.model.load_model(args.modeldir, args.model)
n_rxns = str(model.number_of_reactions())
if args.verbose:
    print_status('Model contains ' + n_rxns + ' reactions')

###############################################################################
# DATABASE QUERY
###############################################################################
# Create log file
log_out = open('log.txt', 'w')

# Iterate through model reactions
if args.verbose:
    print_status('Processing model reactions')
for i, mseed_rxn in enumerate(model.reactions, start=1):
    print('Processing reaction', str(i), 'of', n_rxns,
          end='\r', file=sys.stderr)

    # Check if reaction exists in alias dictionary
    if mseed_rxn not in aliases:
        log_out.write(mseed_rxn + ' not in alias file\n')
        continue

    # Get databases
    for db in aliases[mseed_rxn]:
        if db == 'MetaCyc' or db == 'PlantCyc':
            log_out.write(mseed_rxn + ' has ' + db + ' info. Skipping\n')

        elif db == 'KEGG':
            # May have multiple KO identifiers
            for ko in aliases[db]:
                make_kegg_query(mseed_rxn, ko, log_out)

log_out.close()