#!/usr/bin/env python
from glob import glob
from collections import defaultdict as dd
# import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from goatools import obo_parser
import networkx as nx

go = obo_parser.GODag('../blast/GO/go-basic.obo')
ribo_go_id = 'GO:0003735'

# read in blast results
def read_blast_results(dirpattern, colnames, pident_cutoff, evalue_cutoff):
    df = None
    for i, infile in enumerate(sorted(glob(dirpattern))):
        tmp_df = pd.read_csv(infile, sep='\t', names=colnames)

        tmp_df = tmp_df[((tmp_df['pident']>pident_cutoff) & 
                         (tmp_df['evalue']<evalue_cutoff))]
        if i ==0:
            df = tmp_df
        else:
            df = pd.concat((df, tmp_df))
    return df

# get all parents from a list of GO ids
def find_parents(term_id, go, go_ids_set, _return=False):
    for term_next in go[term_id].parents:
        go_ids_set.update({term_next.id})
        # find all parents
        find_parents(term_next.id, go, go_ids_set)
    if(_return):
        return go_ids_set

# check for presence of a specifc GO id in this list and all parents
def check_for_go_id(go_ids_list, **kwargs):
    term_id_set = set(go_ids_list)
    for _go_id in go_ids_list:
        find_parents(_go_id, go, term_id_set)
    return kwargs['go_id'] in term_id_set

def annotate_subgraph(i, subgraph, busco_hits, swiss_hits):
    density = nx.density(subgraph)
    nnodes = subgraph.number_of_nodes()
    nodes_list = list(subgraph.nodes)
    annotate_nodes = []
    for node in nodes_list:
        busco_row = busco_hits[((busco_hits['matrix'] ==node[0]) & 
                            (busco_hits['part']   ==node[1]))][['Busco_Id', 'Description']]
        swiss_row = swiss_hits[((swiss_hits['matrix'] ==node[0]) & 
                                (swiss_hits['part']   ==node[1]))][['sseqid', 'swissprot_desc']]
        b_id = ''
        b_desc = ''
        s_id = ''
        s_desc = ''
        s_go = ''
        ribo_found = ''

        if busco_row.shape[0] != 0:
            b_id = busco_row['Busco_Id'].iloc[0]
            b_desc = busco_row['Description'].iloc[0]
        if swiss_row.shape[0] != 0:
            s_id = swiss_row['sseqid'].iloc[0]
            s_desc = swiss_row['swissprot_desc'].iloc[0]
            go_lookup = swiss2go[swiss2go['UniProtKB-AC'] == s_id.split('.')[0]]
            if go_lookup.shape[0] != 0:
                s_go = go_lookup['GO'].iloc[0]
                if not type(s_go) is type(np.nan):
                    go_list = s_go.split('; ')
                    ribo_found = check_for_go_id(go_list, go_id=ribo_go_id)
        annotate_nodes.append([i] + list(node) + [len(subgraph[node]), nnodes, density,
                                                b_id, b_desc, s_id, s_desc, s_go, ribo_found])
    return annotate_nodes

# Read in All vs All results (ava) from diaomond blast runs
# cutoffs: chose to retain enough results to keep at least 1 hit for ~90% of partitions
ava_pident_cutoff = 50
ava_evalue_cutoff = 1e-5
ava_colnames = ['qseqid', 'sseqid', 'pident', 'length', 
                'mismatch', 'gapopen', 'qstart', 'qend', 
                'sstart', 'send', 'evalue', 'bitscore']
ava_df = read_blast_results('../blast/diamond_results/*_permissive.txt', 
                            ava_colnames, 
                            ava_pident_cutoff, 
                            ava_evalue_cutoff)

# filter out self-matches
ava_df = ava_df[ava_df['qseqid'] != ava_df['sseqid']]

# explode the ava blast qseq and sseq columns
ava_qseq_colnames = ['matrix1', 'taxon1', 'taxid1', 'part1']
ava_sseq_colnames = ['matrix2', 'taxon2', 'taxid2', 'part2']
ava_qseq = ava_df['qseqid'].str.split(':', expand=True)
ava_sseq = ava_df['sseqid'].str.split(':', expand=True)
ava_qseq.columns = ava_qseq_colnames
ava_sseq.columns = ava_sseq_colnames
ava_df = pd.concat((ava_qseq, ava_sseq, ava_df), axis=1)

# count matrix:partition <-> mscript:matrix:partition 
parts_groupby_columns = ['matrix1', 'part1', 'matrix2', 'part2']
parts_ava_df = ava_df[ava_df['part1'] != ava_df['part2']]
parts = parts_ava_df.groupby(parts_groupby_columns).size().reset_index(name='count')
parts.to_csv('../blast/graphs/partitions_graph.tsv', sep='\t', index=False)
# count matrix:partition <-> mscript:matrix:partition 
taxon_groupby_columns = ['matrix1', 'taxon1', 'matrix2', 'taxon2']
taxa = ava_df.groupby(taxon_groupby_columns).size().reset_index(name='count')
taxa.to_csv('../blast/graphs/taxon_graph.tsv', sep='\t', index=False)

# read in BUSCO results
# columns:
# Busco_Id, Status, Sequence, Score, Length, Description, #Species, #Single-Copy, #Multi-Copy, ProtMedianLength, EvoRate, BiologicalProcesses, InterProDomains
busco_df = pd.read_csv('../blast/graphs/busco_metazoa_results.tsv', sep='\t')
bseq_colnames = ['matrix', 'taxon', 'taxid', 'part']
bseq = busco_df['Sequence'].str.split(':', expand=True)
bseq.columns = bseq_colnames
busco_df = pd.concat((bseq, busco_df), axis=1)
busco_hits = busco_df.groupby(['matrix','part', 'Busco_Id','Description']).size().reset_index(name='count')
busco_hits = busco_hits.sort_values(['matrix', 'part', 'count'], ascending=[True, True, False]).drop_duplicates(['matrix', 'part'])

# read in all vs SwissProt results
swiss_colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'swissprot_desc']
# cutoffs: chose to retain enough results to keep at least 1 hit for ~98% of partitions
swiss_pident_cutoff = 50
swiss_evalue_cutoff = 1e-15
swiss_df = read_blast_results('../blast/swissprot_results/*_swissprot.txt', swiss_colnames, swiss_pident_cutoff, swiss_evalue_cutoff)
swiss_qseq_colnames = ['matrix', 'taxon', 'taxid', 'part']
swiss_qseq = swiss_df['qseqid'].str.split(':', expand=True)
swiss_qseq.columns = swiss_qseq_colnames
swiss_df = pd.concat((swiss_qseq, swiss_df), axis=1)
swiss_hits = swiss_df.groupby(['matrix', 'part', 'sseqid', 'swissprot_desc']).size().reset_index(name='count')
# drop dupes keeps first, which is highest count or best hit if order hasn't changed from sort
swiss_hits = swiss_hits.sort_values(['matrix', 'part', 'count'], ascending=[True, True, False]).drop_duplicates(['matrix', 'part'])
swiss_hits.to_csv('../blast/graphs/swiss_hits.tsv', sep='\t', index=False)

# get swissprot -> GO mappings
swiss2go = pd.read_csv('../blast/graphs/swissprot_to_GO.tsv', sep='\t')

# construct graph from ava partition-level hits
parts_graph = nx.Graph()
parts_graph.add_edges_from([(tuple(x[:2]), tuple(x[2:4]), {'weight':int(x[4])}) for x in parts.values])
# divide graph into separate components, sorted by size
# https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.connected_components.html
components = nx.connected_components(parts_graph)

component_rows = []
components_to_plot = []
tossed_nodes = []
component_number = 0

for component in sorted(components, key=len, reverse=True):
    subgraph = nx.Graph(parts_graph.subgraph(component))
    node_num_edges = {x:len(subgraph[x]) for x in component}
    # bridges = nx.bridges(subgraph)
    vitality = nx.closeness_vitality(subgraph)
    node_num_edges = pd.DataFrame({'nodes':[x for x in component], 'num_edges':[len(subgraph[x]) for x in component]}).sort_values(by='num_edges', ascending=False)
    stat_summary = node_num_edges.describe()
    max_edges_node = node_num_edges.max()['nodes']
    # if the node with most edges has > 2stddev and its removal will split the component into two
    if (node_num_edges['num_edges'].max() > node_num_edges['num_edges'].std()*2+node_num_edges['num_edges'].mean()) and (np.isneginf(vitality[max_edges_node])):
        tossed_nodes.append(max_edges_node)
        subgraph.remove_node(max_edges_node)
        sub_components = nx.connected_components(subgraph)
        for sub_component in sub_components:
            component_rows += annotate_subgraph(component_number, parts_graph.subgraph(sub_component), busco_hits, swiss_hits)
            component_number += 1
    else:
        component_rows += annotate_subgraph(component_number, subgraph, busco_hits, swiss_hits)
        component_number += 1

    if component_number < 9:
        components_to_plot.append(subgraph)

annotated_components_columns = ['component_number', 'matrix', 'partition_name', 'edges',
                                'nodes_in_component', 'component_density', 'BUSCO_ID',
                                'BUSCO_description', 'SwissProt_accession',
                                'SwissProt_description', 'GO_annotations', 'ribo_found']

out = open('../blast/graphs/partition_components_split_annotated.tsv', 'w')
# write header
print('\t'.join(annotated_components_columns), file=out)
#write out data
for row in component_rows:
    print('\t'.join([str(x) for x in row]), file=out) 
out.close()

out = open('../blast/graphs/discarded_nodes.tsv', 'w')
# write header
print('\t'.join(['matrix', 'partition_name']), file=out)
#write out data
for node in tossed_nodes:
    print('\t'.join(node), file=out) 
out.close()

'''
fig = plt.figure(figsize=(12,12))
fig.subplots(3, 3)
for i, (component, ax) in enumerate(zip(components_to_plot, fig.axes)):
    positions = nx.spring_layout(component, k=0.5, iterations=100)
    plt.sca(ax)
    nx.draw_networkx(component, with_labels=False, pos=positions, alpha=0.7,
                     node_size=12, node_color='grey', linewidths=0.4)
    ax.axis('off')
fig.tight_layout()
plt.savefig('../blast/graphs/nine_largest_component_networks.pdf')
'''