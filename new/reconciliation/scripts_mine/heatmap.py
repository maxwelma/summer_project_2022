# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 13:00:38 2022

@author: matil
"""
import pandas as pd
import numpy as np
import seaborn

#Take partition map and output table of yes or no per component per matrix
#make plot showing occurance of component in each matrix

#make dataframe from reshaped partition map (csv) 
#(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\add_new_data_test\part_map_glob_new.csv')
#
def reshape_parts(partition_file):
    import pandas as pd
    from compare_comp_nrs import process_df
    parts = pd.read_csv(partition_file,sep='\t')
    parts,group = process_df(parts)
    #parts = parts[["matrix","component_number"]]
    matrices = parts[['matrix']].drop_duplicates()
    parts = parts.groupby("component_number")['matrix'].apply(list)
    co_ma_df = pd.DataFrame(0, index=np.arange(len(partition_group)),columns = list(matrices.matrix))
    #make true or false df
    for co_nr,df in zip(partition_group.index.values.astype(int), partition_group):
        for matrix in list(matrices.matrix):
            if matrix in df:
                co_ma_df.at[co_nr, matrix] = 1
    return co_ma_df


    
def colour_porifera(df,matrices_por):
    for matrix in matrices_por:
        if matrix in df.columns:
            df[[matrix]] = df[[matrix]].mask(df > 0, 2) #green    
    return df

matrices_por = ["Philippe2009","Nosenko2013_nonribo_9187", "Nosenko2013_ribo_11057", "Nosenko2013_ribo_14615", "Simion2017"]
partitions = colour_porifera(reshape_parts(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\files_mine\part_map_glob.csv'),matrices_por)

#Save as tables
#partitions.to_csv('component_matrix_occurence.cvs',sep='\t')
#partitions.to_excel('component_matrix_occurence.xlsx')

#load table directly

partitions = pd.read_csv('component_matrix_occurence.cvs',sep='\t')

#Set colours for heatmap
colors = ["white","green","orange"]
seaborn.heatmap(partitions, cmap=seaborn.color_palette(colors), cbar=False)

#remove large datasets
small_df = partitions.drop(["Hejnol2009","Borowiec2015_Best108","Borowiec2015_Total1080","Simion2017"],axis=1)
small_df = small_df.loc[(small_df!=0).any(1)]
seaborn.heatmap(small_df,cmap=seaborn.color_palette(colors),cbar=False)

#keep only large datasets
large_df = partitions[["Hejnol2009","Borowiec2015_Best108","Borowiec2015_Total1080","Simion2017"]]
large_df = large_df.loc[(large_df!=0).any(1)]
seaborn.heatmap(large_df,cmap=seaborn.color_palette(colors),cbar=False)

#keep only porifera
colors_pori = ["white","green"]
pori_df = partitions[matrices_por]
pori_df = pori_df.loc[(pori_df!=0).any(1)]
seaborn.heatmap(pori_df,cmap=seaborn.color_palette(colors_pori),cbar=False)



#partitions = pd.read_csv(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\add_new_data_test\part_map_glob_new.csv', sep="\t")
#partitions_all_info = pd.read_csv(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\files_mine\part_map_glob.csv', sep="\t")
#partitions_mc=partitions_all_info[["matrix","component_number"]]
#matrices = partitions_mc[['matrix']].drop_duplicates()
#partition_group = partitions_mc.groupby("component_number")['matrix'].apply(list)
#
##make df of 0s
#co_ma_df = pd.DataFrame(0, index=np.arange(len(partition_group)),columns = list(matrices.matrix))
#
##make true or false df
#for co_nr,df in zip(partition_group.index.values.astype(int), partition_group):
#    for matrix in list(matrices.matrix):
#        if matrix in df:
#            print(co_nr)
#            co_ma_df.at[co_nr, matrix] = 1

#porifera = green

##ctenofora = orange
#matrices_cten = list(matrices.drop(matrices[matrices['matrix'].isin(matrices_por)].index).matrix)
##make true or false df colour based
#for co_nr,df in zip(partition_group.index.values.astype(int), partition_group):
#    for matrix in matrices_por:
#        if matrix in df:
#            print(co_nr)
#            co_ma_df.at[co_nr, matrix] = 1 #green
#    for matrix in matrices_cten:
#        if matrix in df:
#            print(co_nr)
#            co_ma_df.at[co_nr, matrix] = 2 #orange

