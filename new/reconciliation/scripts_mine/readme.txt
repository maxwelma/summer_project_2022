add_comp_nr_1_col.py
#find component number for each sequence in fasta file based on the fasta header name 

add_comp_nrs.py
#find component number for each sequence in blast db file (fasta) based on the fasta header name and output per component fasta files

compare_comp_nrs.py
#from diamond results (*_permissive.txt) use component map (Table S9 from supplementary material OR partition_map_global from R) to find what component maps to each sequence. Make a dataframe with a subset of sequences where at least one sequence doesnt map to any component and  (if any) a subset om sequences that dont match component number.

supp_mat2df.py
#Turn supplementart table s9 into file map for mapping of component numers

colour_trees.sh        
#Bash script for tree colouring

heatmap.py       
#Take partition map and output table of yes or no per component per matrix, make plot showing occurance of component in each matrix (colour coded by paper conclusion)

overlap_rectangles.R
#make figure similar to supplementary figure showing overlaps

shape_partition.py
#For dosReis data, reshape the input partition file to be in desired format
