#!/usr/bin/env python
import pandas as pd

part_components_df = pd.read_table('../blast/graphs/partition_components_annotated.tsv', sep='\t')
part_components_df['SwissProt_description_short'] = part_components_df['SwissProt_description'].str.split(r'e: Full=|;', expand=True).iloc[:,1]
grouped_descriptions = part_components_df.fillna('na').groupby(['component_number', 'nodes_in_component', 'BUSCO_description', 'SwissProt_description_short'], as_index=False).size().to_frame(name='node_count')
grouped_descriptions.to_excel('../blast/graphs/component_annotations.xlsx')