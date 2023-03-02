#!/usr/bin/env python
import re
import operator
import argparse
from functools import reduce
from Bio import AlignIO
from Bio import SeqIO
#from Bio import Alphabet

#gapped_protein = Alphabet.Gapped(Alphabet.IUPAC.protein, '-')

def get_parts(part_file_name):
    parts = {}
    for line in open(part_file_name):
        match = re.match(r"\W*CHARSET\W+([\w\.\-]+)\W*=\W*(\d+)\W*-\W*(\d+)", line, re.IGNORECASE)
        if match is not None:
            part_name, part_start, part_end = match.groups()
            # number non-unique partition names
            parts[part_name] = {'start':int(part_start)-1, 'end':int(part_end)}
    return parts

def pick_parts(alignment, parts, parts_sample):
    sub_sequences = []
    new_parts = []
    curr_length = 1
    for part in parts_sample:
        part_length = parts[part]['end']-parts[part]['start']
        new_parts.append((part, curr_length, curr_length+part_length-1))
        curr_length+=part_length
        sub_sequences.append( alignment[:, parts[part]['start']:parts[part]['end']] )
    return reduce(operator.add, sub_sequences), new_parts

def print_parts(parts, out):
    with open(out, 'w') as out_handle:
        for part in parts:
            print('DNA, {}={}-{};'.format(*part), file=out_handle)

desc = """part_pick_phylip.py
Extracts the specified partitions from a given alignment, outputs them and a new raxml-style  partition file.
"""

parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 prog='pick_parts.py')
parser.add_argument('-p', '--prefix',
                    required=True,
                    help='phylip & nexus file name prefix.')
parser.add_argument('-o', '--out-dir',
                    required=True,
                    help='Directory to put sub_sampled phylip into')
parser.add_argument('-i', '--parts',
                    required=True,
                    help='file with list of partition names to sample, one per line')
args = parser.parse_args()

all_parts = get_parts(args.prefix + '.nex')
parts_subset = [x.rstrip() for x in open(args.parts)]
orig_super_matrix = AlignIO.read(args.prefix + '.phy', 'phylip-relaxed')

new_super_matrix, new_parts = pick_parts(orig_super_matrix, all_parts, parts_subset)
AlignIO.write(new_super_matrix, args.prefix+'_subsampled.phy', 'phylip-relaxed')
print_parts(new_parts, args.prefix+'_subsampled.parts')
