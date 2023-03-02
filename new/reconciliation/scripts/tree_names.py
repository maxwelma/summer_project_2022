#!/usr/bin/env python
from glob import glob
from os import path
import pandas as pd
from dendropy import Tree

taxon_table_file = 'reconciliation/taxonomy_info/taxon_table.tsv'
taxon_table = pd.read_csv(taxon_table_file, sep='\t')

for treefile in glob('trees_new/*.contree'):
    # e.g. 'Borowiec2015-Total1080.WAG.contree' or 'Dunn2008.poisson_C60.contree'
    manuscript_matrix = path.split(treefile)[1].split('.', 1)[0]
    tmp = manuscript_matrix.split('-')
    if len(tmp) == 1:
        manuscript, matrix = tmp[0], tmp[0]
    else:
        manuscript, matrix = tmp
    this_taxon_table = taxon_table[taxon_table['original_matrix'].str.contains(manuscript)]
    tree = Tree.get_from_path(treefile, schema='newick')
    taxa = [x.replace(' ', '_') for x in tree.taxon_namespace.labels()]
    for taxon in taxa:
        if this_taxon_table[this_taxon_table['relabelled_name']==taxon].shape[0] == 0:
            relab = this_taxon_table[this_taxon_table['matrix_name']==taxon]
            if relab.shape[0] == 0:
                print(manuscript, matrix, taxon)