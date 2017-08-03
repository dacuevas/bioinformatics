#!/usr/local/bin/python3
# reactions_to_ko.py
# Use PyFBA reaction IDs to obtain KEGG KO IDs and their associated pathways.
#
# Author: Daniel A Cuevas (dcuevas08.at.gmail.com)
# Created on 02 Aug 2017
# Updated on 03 Aug 2017
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


def make_kegg_query(ms, ko, log, save_kos):
    """
    Setup a query to the KEGG API

    :param ms: ModelSEED reaction ID
    :type ms: str
    :param ko: KEGG KO ID
    :type ko: str
    :param log: Log output file handle
    :type log: File
    :param save_kos: Set to add KO IDs to
    :type save_kos: set
    :return: None
    """
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
                  + str(response.status_code) + '\n')
        return

    # Check that response is not empty
    elif response.text.strip('\n') == '':
        log.write('Response was empty for KO ' + ko + '\n')
        return

    # Parse response
    parse_kegg_response(response.text, ms, ko, log, save_kos)


def parse_kegg_response(res, ms, ko, log, save_kos):
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
    :param ms: ModelSEED reaction ID
    :type ms: str
    :param ko: KEGG KO ID
    :type ko: str
    :param log: Log output file handle
    :type log: File
    :param save_kos: Set to add KO IDs to
    :type save_kos: set
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

        if curr_section == 'ORTHOLOGY':
            ko_id = re.match(r'K\d+', kegg_data).group(0)
            save_kos.add(ko_id)


def make_bigg_query(ms, bi, log):
    """
    Setup a query to the BiGG API

    :param ms: ModelSEED reaction ID
    :type ms: str
    :param bi: BiGG reaction ID
    :type bi: str
    :param log: Log output file handle
    :type log: File
    :return: KO IDs and EC numbers
    :rtype: dict
    """
    # Set resource path
    resource = 'universal/reactions/' + bi

    # Issue request
    full_path = os.path.join(BIGG_BASE_URL, resource)
    response = requests.get(full_path)

    # Check that status code is 200 = good
    if response.status_code != 200:
        log.write('There was an error with the request '
                  'regarding ModelSEED reaction ' + ms
                  + ' and BiGG ' + bi + ', status code = '
                  + str(response.status_code) + '\n')
        return None

    # Check that response is not empty
    elif response.text.strip('\n') == '':
        log.write('Response was empty for BiGG ' + bi + '\n')
        return None

    # Parse response
    return parse_bigg_response(response.json(), ms, bi, log)


def parse_bigg_response(res, ms, bi, log):
    """
    Parse the BiGG API response text. Text is all plain text in JSON format.
    The fields of interest are the KEGG Reaction ID or the EC number.

    :param res: API JSON response
    :type res: dict
    :param ms: ModelSEED reaction ID
    :type ms: str
    :param bi: BiGG reaction ID
    :type bi: str
    :param log: Log output file handle
    :type log: File
    :return: KO IDs and EC numbers
    :rtype: dict
    """
    data = {'KO': set(), 'EC': set()}  # Data to return

    # Check if any database info exists
    db_info = res['database_links']
    if len(db_info) == 0:
        log.write('No database info for BiGG ' + bi + '\n')

    # Check for KEGG
    elif 'KEGG Reaction' in db_info:
        # May have multiple KO identfiers
        for item in db_info['KEGG Reaction']:
            if 'id' not in item:
                log.write('KEGG reaction found but no ID for BiGG '
                          + bi + ' and ModelSEED ' + ms + ' \n')
            data['KO'].add(item['id'])

    # Check for EC number of KEGG does not exist
    elif 'EC Number' in db_info:
        # May have multiple EC numbers
        for item in db_info['EC Number']:
            if 'id' not in item:
                log.write('EC number found but no ID for BiGG '
                          + bi + ' and ModelSEED ' + ms + ' \n')
            data['EC'].add(item['id'])

    # No KEGG or EC
    else:
        log.write('No KEGG Reaction or EC Number for BiGG '
                  + bi + ' and ModelSEED ' + ms + ' \n')

    return data


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
reactions_file = os.path.join(args.mseeddir,
                              'Biochemistry/reactions.tsv')
if not os.path.isfile(alias_file):
    print_status('ModelSEED Reactions_Aliases.tsv file does not exist')
    exit_script()

if not os.path.isfile(reactions_file):
    print_status('ModelSEED reactions.tsv file does not exist')
    exit_script()

# Set base URLs
KEGG_BASE_URL = 'http://rest.kegg.jp/'
BIGG_BASE_URL = 'http://bigg.ucsd.edu/api/v2/'

###############################################################################
# BEGIN PROCESSING
###############################################################################
aliases = {}  # aliases will contain Reactions_Aliases.tsv info
save_kos = set()  # save_kos will contain all KO IDs to print out

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

# Load ModelSEED reactions file
# This file does not have database information but does provide an alias
# The default will be to guess it is a KEGG KO ID if it begins with 'R'. If
# not, then it might be a BiGG ID
if args.verbose:
    print_status('Parsing reactions file')
with open(reactions_file, 'r') as f:
    header = f.readline()
    for l in f:
        l = l.rstrip('\n')
        contents = l.split('\t')
        ms, alias = contents[:2]

        # Check alias name
        if alias.startswith('R'):
            db = 'KEGG'
        else:
            db = 'BiGG'

        # Check if this reaction was not in the aliases file
        if ms not in aliases:
            aliases[ms] = {db: set()}
        elif db not in aliases[ms]:
            aliases[ms][db] = set()

        aliases[ms][db].add(alias)

if args.verbose:
    print_status('Reactions file parsing complete')

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

    # Get databases
    for db in aliases[mseed_rxn]:
        #######################################
        # METACYC & PLANTCYC
        #######################################
        if db == 'MetaCyc' or db == 'PlantCyc':
            log_out.write(mseed_rxn + ' database is ' + db + '. Skipping\n')

        #######################################
        # KEGG
        #######################################
        elif db == 'KEGG':
            # May have multiple KO identifiers
            for ko in aliases[mseed_rxn][db]:
                make_kegg_query(mseed_rxn, ko, log_out, save_kos)

        #######################################
        # BIGG
        #######################################
        elif db == 'BiGG':
            # May have multiple BiGG identifiers
            for bigg in aliases[mseed_rxn][db]:
                # Get dictionary of KO IDs and EC numbers
                bigg_info = make_bigg_query(mseed_rxn, bigg, log_out)
                if bigg_info == None:
                    continue

                # Make queries to KEGG database with KO
                for ko in bigg_info['KO']:
                    make_kegg_query(mseed_rxn, ko, log_out, save_kos)

                # Make queries to KEGG database with EC
                for ec in bigg_info['EC']:
                    make_kegg_query(mseed_rxn, ec, log_out, save_kos)

        else:
            log_out.write(mseed_rxn + ' database is ' + db
                          + ' and not supported\n')

###############################################################################
# PRINT OUT DATA
###############################################################################
print('\n'.join(save_kos))
log_out.close()
print_status('Script complete!')
