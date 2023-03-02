#!/usr/bin/env python
import os
import re
import sys
import argparse
from collections import defaultdict
from os.path import  abspath, basename, dirname, join, splitext
from Bio import AlignIO
from Bio import SeqIO
#from Bio import Alphabet

def get_parts(part_file_name):
    parts = []
    partnames = defaultdict(lambda: 0)
    for line in open(part_file_name):
        match = re.match(r"CHARSET\W+([\w\.\-]+)\W*=\W*(\d+)\W*-\W*(\d+)", line, re.IGNORECASE)
        if match is not None:
            part_name, part_start, part_end = match.groups()
            # number non-unique partition names
            if partnames[part_name] ==0:
                partnames[part_name]+=1
            else:
                part_name = "{}_{:02d}".format(part_name, partnames[part_name])
                partnames[part_name]+=1
            # zero-index the intervals 
            parts.append((part_name, int(part_start)-1, int(part_end)))
    return parts

def read_taxon_table(taxon_table_fn):
    names_map = defaultdict(lambda: 'na')
    table = open(taxon_table_fn, 'r')
    header = table.readline().split()
    tax_id = header.index('ncbi_tax_id')
    new_name = header.index('relabelled_name')
    for line in table:
        tmp = line.rstrip('\n').split('\t')
        names_map[tmp[new_name]] = tmp[tax_id]
    return names_map

desc = """phylippart_to_multifasta.py
Simple script to take a phylip and nexus partitions file pair and convert it into a miltifasta where each partition for each taxon is a separate sequence.
"""

parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 prog='prepare_phylip.py')
parser.add_argument('-p', '--prefix',
                    required=True,
                    help='phylip & nexus file name prefix.')
parser.add_argument('-o', '--out-dir',
                    required=True,
                    help='Directory to put converted fasta')
parser.add_argument('-t', '--taxon-table',
                    required=True,
                    help='path to taxon_table.tsv')
args = parser.parse_args()

prefix = args.prefix
out_dir = args.out_dir
part_file_name = prefix + '.nex'
phylip_file_name = prefix + '.phy'
matrix  = basename(prefix)

name2taxid = read_taxon_table(args.taxon_table)
parts = get_parts(part_file_name)
#gapped_protein = Alphabet.Gapped(Alphabet.IUPAC.protein, '-')
#gapped_protein

alignment = AlignIO.read(phylip_file_name, 'phylip-relaxed')
fasta_out = open(join(out_dir, '{}.fa'.format(matrix)),  'w')
for part_name, start, stop in parts:
    for i, taxon in enumerate(alignment):
        tmp_seq = taxon[start:stop]
        #tmp_seq.seq.alphabet = gapped_protein
        old_seq = tmp_seq.seq
        tmp_seq.seq = tmp_seq.seq.ungap()
        if len(tmp_seq) > 1:
            tmp_seq.id = ':'.join([matrix, taxon.id, name2taxid[taxon.id], part_name])
            SeqIO.write(tmp_seq, fasta_out, 'fasta')
fasta_out.close()