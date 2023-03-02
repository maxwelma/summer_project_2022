# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:26:14 2022

@author: matil
"""


# to do: compare ids in ava to matrix-partition in partitions and add component numbers to ava as column

#for sseqid and for qseqid 
#if name (xseqid) in ava is in list of names in partitions (matrix_partition) add component nr as new column (xcomp_nr)

#for i,name in zip(ava_df.index.values.astype(int), ava_df.sseqid):
#    for co_nr , name_list in zip(partitions_df.component_number, partitions_df.matrix_partition):
#        if name in name_list:
#            ava_df.at[i,'s_comp_nr'] = co_nr
#
#for i,name in zip(ava_df.index.values.astype(int), ava_df.qseqid):
#    for co_nr , name_list in zip(partitions_df.component_number, partitions_df.matrix_partition):
#        if name in name_list:
#            ava_df.at[i,'q_comp_nr'] = co_nr
#            
#ava_df[['s_comp_nr','q_comp_nr']]
#
#non_match_comp_df = ava_df[ava_df['s_comp_nr'] != ava_df['q_comp_nr']]

#partitions_df["matrix_partition"]= partitions_df.matrix_partition.str.replace(" ","")
#partitions_df["matrix_partition"] = partitions_df.matrix_partition.str.split(", ")
#
#partitions_df_expanded =  partitions_df.explode('matrix_partition')



    j = (i + 1) / len(ava_df)
    sys.stdout.write('\r')
    # the exact output you're looking for:
    sys.stdout.write("[%-20s] %d%%" % ('='*int(20*j), 100*j))
    #sys.stdout.flush()
    sleep(0.25)
    
    #ava_df.at[i,'qcomp_nr'] = co_nr
    
    # to do: compare ids in ava to matrix-partition in partitions and add component numbers to ava as column

#for sseqid and for qseqid 
#if name (xseqid) in ava is in list of names in partitions (matrix_partition) add component nr as new column (xcomp_nr)

#for i,name in zip(ava_df.index.values.astype(int), ava_df.sseqid):
#    for co_nr , name_list in zip(partitions_df.component_number, partitions_df.matrix_partition):
#        if name in name_list:
#            ava_df.at[i,'s_comp_nr'] = co_nr
#
#for i,name in zip(ava_df.index.values.astype(int), ava_df.qseqid):
#    for co_nr , name_list in zip(partitions_df.component_number, partitions_df.matrix_partition):
#        if name in name_list:
#            ava_df.at[i,'q_comp_nr'] = co_nr
#            
#ava_df[['s_comp_nr','q_comp_nr']]
#
#non_match_comp_df = ava_df[ava_df['s_comp_nr'] != ava_df['q_comp_nr']]
    
#partitions_df["matrix_partition"]= partitions_df.matrix_partition.str.replace(" ","")
#partitions_df["matrix_partition"] = partitions_df.matrix_partition.str.split(", ")
#
#partitions_df_expanded =  partitions_df.explode('matrix_partition')

#for i,qname,sname in zip(ava_df.index.values.astype(int), ava_df.qseqid,ava_df.sseqid):
#    ava_df = add_co_nr(i,ava_df,'scomp_nr','sname',sname)
#    ava_df = add_co_nr(i,ava_df,'qcomp_nr','qname',qname)
#    #print(i)
#    #ava_df.at[i,'qcomp_nr'] = co_nr
    
    #rename ids as matrix and partition, remove taxa and positions
#supp_mat_df = supp_mat_df[["component_number","matrix","partition"]]
#supp_mat_df_sep = supp_mat_df.copy()
#supp_mat_df_sep["matrix_partition"] = supp_mat_df["matrix"] + supp_mat_df["partition"].astype(str)+","
#supp_mat_df["matrix_partition"] = supp_mat_df["matrix"] + supp_mat_df["partition"].astype(str)
#supp_mat_df=supp_mat_df[["component_number","matrix_partition"]] #expanded (original) file
#supp_mat_df_sep=supp_mat_df_sep[["component_number","matrix_partition"]] #non expanded file, to be grouped
#group_df = supp_mat_df_sep.groupby('component_number').agg({'matrix_partition':'sum'}) #group, non expanded file

#headers = []
#with open(infile, "r") as f:
#    for record in SeqIO.parse(f, "fasta"):
#        headers.append(record.description)
#headers_df = pd.DataFrame(headers, columns={"fasta_header"})


