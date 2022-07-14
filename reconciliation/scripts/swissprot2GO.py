#!/usr/bin/env python
import pandas as pd

swiss_hit_ids = [x.split('.')[0] for x in open('../blast/swissprot_results/unique_gene_ids.txt')]
# get from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
uniprotkb = pd.read_csv('idmapping_selected.tab', names=["UniProtKB-AC", "GO"], usecols=[0,6], index_col=0, sep='\t', dtype={"UniProtKB-ID": 'object', "GO": 'object'})

uniprotkb.loc[swiss_hit_ids].to_csv('../blast/graphs/swissprot2GO.tsv', sep='\t')
