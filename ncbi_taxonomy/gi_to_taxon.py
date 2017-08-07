from __future__ import print_function, absolute_import, division
import sys
import os
import time
import datetime
import argparse
import taxon

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


###############################################################################
# ARGUMENT PARSING
###############################################################################
desc = 'Collect taxonomy information for GI numbers'
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('gifile', help='Input file of GI numbers. Gi numbers '
                    'should be second value on each line')
parser.add_argument('outdir', help='Output directory')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose output')

args = parser.parse_args()

# Check that input file exists
if not os.path.isfile(args.gifile):
    print(args.gifile, 'does not exist', file=sys.stderr)
    parser.print_usage()
    exit_script()

# Check that output directory exists
if not os.path.isdir(args.outdir):
    print(args.outdir, 'does not exist', file=sys.stderr)
    parser.print_usage()
    exit_script()

# Global variables
gi_file = args.gifile
out_dir = args.outdir
vbs = args.verbose

###############################################################################
# LOAD INPUT FILE
###############################################################################
# Create log file
log = open(os.path.join(args.outdir, 'log.txt'), 'w', buffering=1)
log.write(timestamp() + ' Starting script\n')
gi_counts = {}  # Hold gi numbers and counts

# Read in GI file
print_status('Loading input file')
with open(args.gifile) as f:
    for li, l in enumerate(f, start=1):
        if vbs:
            print('Reading entry', li, end='\r', file=sys.stderr)
        gi = l.split()[1].split('|')[1]

        if gi in gi_counts:
            gi_counts[gi] += 1

        else:
            gi_counts[gi] = 1

print_status('Finished input file')
print_status('Loaded {} unique GIs'.format(len(gi_counts)))

###############################################################################
# LOAD NCBI TAXONOMY FILES
###############################################################################
print_status('Loading NCBI taxonomy files')
print_status('Loading GI => TAXID database')

gi_to_taxid = taxon.readGiTaxId('nucl')
print_status('GI => TAXID database loaded')
print_status('Converting GIs to TAXIDs')

taxids = {}
for i, gi in enumerate(gi_counts, start=1):
    if vbs:
        print('Converted', i, 'out of ', len(gi_counts), 'GI values',
              end='\r', file=sys.stderr)
    try:
        tid = gi_to_taxid
        taxids[gi] = tid

    except KeyError:
        msg = 'GI ' + gi + ' not found in gi_taxid file'
        log.write(msg + '\n')
        continue

print_status('GI to TAXID conversion complete')

gi_to_taxid.clear()  # Clear to remove from memory

print_status('Loading TAXID => NAME database')
names, blastnames = taxon.readNames()

# Not using blastnames currently, so deleting
blastnames.clear()
print_status('TAXID => NAME database loaded')

print_status('Loading taxonomy info database')
taxa = taxon.readTaxa()

print_status('Taxonomy info database loaded')

###############################################################################
# CONNECT GI TO TAXONOMY INFO
###############################################################################
all_data = {}  # Holds all data for output

num_taxid = len(taxids)
print_status('Gathering taxonomy information')
for ti, tax_id in enumerate(taxids.values(), start=1):
    if vbs:
        print('Working on tax ID ', ti, 'out of ', num_taxid,
              end='\r', file=sys.stderr)
    try:
        curr_node = taxa[tax_id]

    except KeyError:
        msg = 'TAXID ' + tax_id + ' not found in nodes file'
        log.write(msg + '\n')
        continue

    start_tid = tax_id
    all_data[start_tid] = {'species': None,
                           'genus': None,
                           'family': None,
                           'order': None,
                           'class': None}

    # Find all taxonomy hierarchy for this tax id
    # Loop until Phylum is reached. Phylum is right above Class
    # End at Domain in case
    while curr_node.rank == 'phylum' and curr_node.parent != 1:
        curr_name = ''
        try:
            curr_name = names[curr_node.taxid].name

        except KeyError:
            msg = 'TAXID ' + curr_node.taxid + ' not found in names file'
            log.write(msg + '\n')

        # Set name
        all_data[start_tid][curr_node.rank] = curr_name

        # Get parent and repeat
        curr_node = taxa[curr_node.parent]

print_status('Completed taxonomy information')

###############################################################################
# OUTPUT
###############################################################################
out_file = os.path.join(out_dir, 'tax_info.tsv')
print_status('Creating output file ' + out_file)

with open(out_file, 'w') as f:
    # Header info
    f.write('gi\tcount\tspecies\tgenus\tfamily\torder\tclass\n')

    for gi, tax_id in taxids.items():
        count = gi_counts[gi]
        species = all_data[tax_id]['species']
        genus = all_data[tax_id]['genus']
        family = all_data[tax_id]['family']
        order = all_data[tax_id]['order']
        clss = all_data['class']

        # Check if any are None; set to some default value
        default = ''
        if species is None:
            species = default
        if genus is None:
            genus = default
        if family is None:
            family = default
        if order is None:
            order = default
        if clss is None:
            clss = default

        f.write('\t'.join([gi, count, species, genus, family, order, clss]) +
                '\n')

log.write(timestamp() + ' Script complete\n')
log.close()
print_status('Script complete!')
