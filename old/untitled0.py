# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:47:17 2022

@author: matil
"""

#subset moody et al table with gene information
import pandas as pd

in_df = pd.read_excel(r"C:\Users\matil\Downloads\summer_proj\Moody_et_al\elife-66695-supp1-v2.xlsx", header=1,sep="\t",sheet_name='Marker_gene_screening_sum_full')

in_df = in_df.loc[in_df['MarkerGene'].isin(['p0131', 'p0151', 'p0159', 'p0174', 'p0181', 'p0287', 'p0306', 'p0364','p0000', 'p0011', 'p0020', 'p0091', 'p0094', 'p0202','p0073','p0093'])]

in_df.to_excel(r"C:\Users\matil\Downloads\summer_proj\Moody_et_al\elife-66695-supp1-subset.xlsx", index=False,header=True)

in_df = in_df.drop_duplicates()
in_df.to_csv(r"C:\Users\matil\Downloads\summer_proj\Moody_et_al\elife-66695-supp1-subset.txt", index=False,header=True,sep="\t")

