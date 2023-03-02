# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 13:00:38 2022

@author: matil
"""
import pandas as pd
import numpy as np
import seaborn
partitions = pd.read_csv(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\add_new_data_test\part_map_glob_new.csv', sep="\t")
partitions_all_info = pd.read_csv(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\files_mine\part_map_glob.csv', sep="\t")
partitions_mc=partitions_all_info[["matrix","component_number"]]
matrices = partitions_mc[['matrix']].drop_duplicates()
partition_group = partitions_mc.groupby("component_number")['matrix'].apply(list)

co_ma_df = pd.DataFrame(0, index=np.arange(len(partition_group)),columns = list(matrices.matrix))


for co_nr,df in zip(partition_group.index.values.astype(int), partition_group):
    for matrix in list(matrices.matrix):
        if matrix in df:
            print(co_nr)
            co_ma_df.at[co_nr, matrix] = 1
            
co_ma_df.to_csv('component_matrix_occurence.cvs',sep='\t')
co_ma_df.to_excel('component_matrix_occurence.xlsx')

seaborn.heatmap(co_ma_df)
small_df = co_ma_df.drop(["Hejnol2009","Borowiec2015_Best108","Borowiec2015_Total1080","Simion2017"],axis=1)
small_df = small_df.loc[(small_df!=0).any(1)]
seaborn.heatmap(small_df,cmap=seaborn.color_palette("hls", 2))