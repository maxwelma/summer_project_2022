# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:44:27 2022

@author: matil
"""

#Idea: add component number for each seq based on the blast query name and export to fastas
import pandas as pd
from Bio import SeqIO
from compare_comp_nrs import add_co_nr,process_df
from glob import glob
import argparse
# read in blast results, function copied from blast_partitions_graph.py
def read_blast_results(dirpattern, colnames, pident_cutoff, evalue_cutoff):
    df = None
    for i, infile in enumerate(sorted(glob(dirpattern))):
        tmp_df = pd.read_csv(infile, sep='\t', names=colnames)

        tmp_df = tmp_df[((tmp_df['pident']>pident_cutoff) & 
                         (tmp_df['evalue']<evalue_cutoff))]
        if i ==0:
            df = tmp_df
        else:
            df = pd.concat((df, tmp_df),ignore_index=True)
    return df
#add component number
def edit_df(df,partition_df):
    for i,name in zip(df.index.values.astype(int), df.name):
        df = add_co_nr(i,df,'comp_nr','name',name,partition_df)
        print('\rProcessing row:'+str(i)+'/'+str(len(df)))
    return df

#read fasta headers into dataframe
def read_blast_db(dirpattern,colnames):
    l1 = []
    l2 = []
    with open(dirpattern, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            l1.append(record.name)
            l2.append(record.seq)
                #print(record)
#        if i ==0:
#            l = l
#        else:
#            df = pd.concat((l, tmp_l),ignore_index=True)
    df = pd.DataFrame({colnames[0]:l1,colnames[1]:l2})
    return df
  
def main():
    desc = """Script that
    1)Reads in fasta format database 
    2) adds component number for each seq based on the blast query 
    3) exports sequences associated with each component to separate fasta files"""
    
    parser = argparse.ArgumentParser(description=desc,
                                 usage='%(prog)s --in-file file.fasta --out-dir out_dir/ --part-map path/to/partition_map.csv', 
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 prog='make_fasta.py.py')
    parser.add_argument('-i', '--in-file',
                        required=True,
                        help='fasta format blast database (see reconciliation/db/animal_root.fa)')
    parser.add_argument('-o', '--out-dir',
                        required=True,
                        help='Directory to put resulting fasta files')
    parser.add_argument('-p', '--part-map',
                        required=True,
                        help='Csv format table of matrixes and partitions and corresponding component numbers. (see files_mine/part_map_glob.csv')
    args = parser.parse_args()
    
    infile = args.in_file
    outpath = args.out_dir
    partitions_df = pd.read_csv(args.out_dir, sep="\t")
    
#    #import fasta headers from blast queries to use as name
#    infile = r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\reconciliation\reconciliation\blast\db\animal_root.fa'
#    #set outpath for fastas
#    outpath = r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\files_mine\fasta'
#    #import partition map 
#    partitions_df = pd.read_csv(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\files_mine\part_map_glob.csv', sep="\t")
#    
    #make df from infile
    colnames = ["fasta_header","seq"]
    headers_df = read_blast_db(infile,colnames)
    
    #add name column
    headers_df = headers_df.join(headers_df.fasta_header.str.split(':', expand = True))[[colnames[0],colnames[1],0,3]]
    headers_df["name"] = headers_df[0] + headers_df[3]
    headers_df = headers_df.drop([0,3],axis=1)
    
    #process partition map
    partitions_df,group_df_r = process_df(partitions_df)
    
    #add component numbers
    headers_df_processed = edit_df(headers_df.copy(),partitions_df)
    
    #split dataframe on component nr
    dfs = dict(tuple(headers_df_processed.groupby('comp_nr')))
    
    #write outfile
    #list of lists of turned into biopython record
    records = []
    for i in dfs:
        sep_rec = []
        for l,j in zip(dfs[i].fasta_header,dfs[i].seq):
            sep_rec.append(SeqIO.SeqRecord(j, id = l,description="Component_nr:"+str(int(i))))
        records.append([sep_rec,i])
        
    for i in range(0,len(records)):
        SeqIO.write(records[i][0], outpath+'/component'+str(int(records[i][1]))+".fasta", "fasta")


if __name__ == "__main__":
    main()