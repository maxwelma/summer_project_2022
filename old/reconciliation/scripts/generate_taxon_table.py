#!/usr/bin/env python3
# adapted from :
# https://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython
import re
import os
import sys
import glob
import urllib
import pickle
from collections import defaultdict
from nexus import NexusReader # https://pypi.org/project/python-nexus/
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'casey.dunn@yale.edu'

clades = { # clades to bin samples into
    'Cnidaria': 'Cnidaria',
    'Bilateria': 'Bilateria',
    'Placozoa': 'Placozoa',
    'Porifera': 'Porifera',
    'Ctenophora': 'Ctenophora',
    'Choanoflagellida': 'Choanoflagellida',
    'Ministeria': 'Filasterea',
    'Capsaspora': 'Filasterea',
    'Ichthyosporea': 'Ichthyosporea',
    'Fungi': 'Fungi',
}

def get_pickles():
    """be kind to API, don't ask for what we already know"""
    if os.path.isfile('../taxonomy_info/known_ids.pickle'):
        k_ids = pickle.load(open('../taxonomy_info/known_ids.pickle', 'rb'))
    else:
        k_ids = defaultdict(lambda: {})

    if os.path.isfile('../taxonomy_info/known_tax.pickle'):
        k_tax = pickle.load(open('../taxonomy_info/known_tax.pickle', 'rb'))
    else:
        k_tax = defaultdict(lambda: {})
    return k_ids, k_tax

def get_names(file_name):
    """read in a phy or nex file, return the taxa names found"""
    names = []
    if file_name.endswith('.phy'):
        phy = SeqIO.parse(file_name, 'phylip-relaxed')
        names = [t.name for t in phy]
    elif file_name.endswith('.nex'):
        nex = NexusReader(file_name)
        try:
            names = [t for t in nex.taxa]
        except AttributeError as e:
            # when there isn't a separate taxa block, get names from sequences
            names = [block[0] for block in nex.blocks['data']]
    return names

def sanitize_name(name):
    # trim trailing characters
    new_name = name.rstrip('_. ')
    # trim numeric only suffixes 
    number_trimmed = re.match(r'(.*)_\d+$', new_name)
    if number_trimmed is not None:
        new_name = number_trimmed.groups()[0]
    # deal with sp
    tmp_split = new_name.split('_')
    if len(tmp_split) > 1 and (tmp_split[1] == 'sp' or tmp_split[1] == 'species'):
        new_name = name.split('_')[0]
    # please forgive my hack :-/
    if new_name.startswith('Mertensiid'):
        new_name = "Undescribed mertensiid sp. 3"
    return re.sub('[ _-]', '+', new_name)

def get_tax_id(query_name):
    """sanitize query name, search for and return NCBI ID"""
    search_term = sanitize_name(query_name)
    search = Entrez.esearch(term = search_term, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    try:
        tax_id = record['IdList'][-1]
    except IndexError as e:
        print("Query {} ({}) not found.".format(query_name, search_term))
        tax_id = 'nan'
    return tax_id

def get_tax_data(tax_id):
    """Fetch taxonomy info with NCBI ID"""
    if tax_id == 'nan':
        return [{'Lineage':'na'}]
    try:
        search = Entrez.efetch(id = tax_id, db = "taxonomy", retmode = "xml")
        return Entrez.read(search)
    except Exception as e:
        print("taxonomy lookup error: ", e)
    return [{'Lineage':'na'}]

def get_clade(lineage):
    for clade in clades:
        if clade in lineage:
            return clades[clade]
    return 'na'

def get_manual_entries():
    fh = open('../taxonomy_info/manual_taxonomy_map.tsv')
    fh.readline() # discard header
    manual = defaultdict(lambda: {})
    for line in fh:
        tmp = line.rstrip().split('\t')
        manual[tmp[0]][tmp[1]] = tmp[2]
    return manual

if __name__ == '__main__':

    known_ids, known_tax = get_pickles()
    manual_entries = get_manual_entries()
    table_columns = ["original_matrix", "relabelled_name", "clade_assignment",
                     "ncbi_tax_id", "ncbi_taxonomy", "matrix_name"]
    sep = '\t'
    table_out = open("../taxonomy_info/taxon_table.tsv", 'w')
    print(sep.join(table_columns), file=table_out)

    # glob hack to match nex and phy
    for file_name in sorted(glob.glob('../considered_data/**/*.[pn][eh][xy]', recursive=True)):
        tmp_meta = file_name.split('/')
        manuscript = tmp_meta[2]
        names_list = []
        print("Reading file", file_name)
        try:
            names_list = get_names(file_name)
        except Exception as e:
            print("Error loading", file_name, ":\n", e)
        for original_name in names_list:
            new_name = 'na'
            manual = False
            if original_name in manual_entries[manuscript]:
                manual = True
                new_name = manual_entries[manuscript][original_name]
                print('Renaming "{}" to "{}".'.format(original_name, new_name))
                known_ids[manuscript][original_name] = get_tax_id(new_name)
            
            if original_name not in known_ids[manuscript]:
                known_ids[manuscript][original_name] = get_tax_id(original_name)

            if known_ids[manuscript][original_name] not in known_tax[manuscript]:
                known_tax[manuscript][known_ids[manuscript][original_name]] = get_tax_data(known_ids[manuscript][original_name])[0]

            # if there is an NCBI Scientific name
            if not manual:
                new_name = original_name.replace('.', '')
                new_name = original_name.rstrip('_.')

            this_clade = get_clade(known_tax[manuscript][known_ids[manuscript][original_name]]['Lineage'])

            print(sep.join([ file_name, new_name, this_clade, known_ids[manuscript][original_name],
                                known_tax[manuscript][known_ids[manuscript][original_name]]['Lineage'],
                                original_name]), file=table_out)

    table_out.close()

    # can't pickle lambda functions, convert defaultdicts to dicts
    pickle.dump(dict(known_ids), open('../taxonomy_info/known_ids.pickle', 'wb'))
    pickle.dump(dict(known_tax), open('../taxonomy_info/known_tax.pickle', 'wb'))
