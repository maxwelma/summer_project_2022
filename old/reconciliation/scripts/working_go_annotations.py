#!/usr/bin/env python
from goatools import obo_parser
import pandas as pd

def find_parents(term_id, go, go_ids_set, _return=False):
    for term_next in go[term_id].parents:
        go_ids_set.update({term_next.id})
        # find all parents
        find_parents(term_next.id, go, go_ids_set)
    if(_return):
        return go_ids_set

def check_for_go_id(go_ids_list, **kwargs):
    term_id_set = set(go_ids_list)
    for _go_id in go_ids_list:
        find_parents(_go_id, go, term_id_set)
    return kwargs['go_id'] in term_id_set


go = obo_parser.GODag('../blast/GO/go-basic.obo')
components = pd.read_csv('../blast/graphs/partition_components_annotated.tsv', sep='\t')
components_go = components['GO_annotations'].str.split('; ').fillna('')
melted_gos = pd.melt(components['GO_annotations'].str.split('; ', expand=True))['value'].dropna()
print(melted_gos.apply(lambda x: go[x].namespace).value_counts())
print(melted_gos.apply(lambda x: go[x].name).value_counts()[0:20])