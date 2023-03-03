
#### Steps

working directory is `reconciliation/scripts`

- run `./generate_taxon_table.py` to generate taxon table.
- run `./gen_phylips.sh` to generate phylips.
- run `./gen_fastas.sh` to generate fastas.
- generate single fasta with `cat ../blast/queries/*.fa > ../blast/db/animal_root.fa`
- on your own:
    - run blast+ against swissprot, example commands in `../blast/swissprot_results/blast_jobs.txt`
    - run diamondBLAST all vs all, examples in `../blast/diamond_results/blast_jobs.txt`
        - make db with `cd ../blast/db && diamond makedb --in animal_root.fa -d animal_root_diamond`
    - run BUSCO:
        ```
        run_BUSCO.py -c 20 -o animal_root_metazoa_busco -m prot -l $EBROOTBUSCO/datasets/metazoa_odb9 -i ../db/animal_root.fa
        run_BUSCO.py -c 20 -o animal_root_eukaryota_busco -m prot -l $EBROOTBUSCO/datasets/eukaryota_odb9 -i ../db/animal_root.fa
        ```
- run `./blast_partitions_graph.py` to generate summary graphs of the diamond blast results in `../blast/graphs`
- run `./merge_busco_results.py` to generate a table of busco results, `../blast/graphs/busco_metazoa_results.tsv`