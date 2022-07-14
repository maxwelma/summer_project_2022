 echo ../../data_processed/matrices/Borowiec2015_Best108; ./phylippart2multifasta.py -p ../../data_processed/matrices/Borowiec2015_Best108 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Borowiec2015_Total1080; ./phylippart2multifasta.py -p ../../data_processed/matrices/Borowiec2015_Total1080 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Chang2015; ./phylippart2multifasta.py -p ../../data_processed/matrices/Chang2015 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Dunn2008; ./phylippart2multifasta.py -p ../../data_processed/matrices/Dunn2008 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Hejnol2009; ./phylippart2multifasta.py -p ../../data_processed/matrices/Hejnol2009 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Moroz2014_3d; ./phylippart2multifasta.py -p ../../data_processed/matrices/Moroz2014_3d -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Nosenko2013_nonribo_9187; ./phylippart2multifasta.py -p ../../data_processed/matrices/Nosenko2013_nonribo_9187 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Nosenko2013_ribo_11057; ./phylippart2multifasta.py -p ../../data_processed/matrices/Nosenko2013_ribo_11057 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Nosenko2013_ribo_14615; ./phylippart2multifasta.py -p ../../data_processed/matrices/Nosenko2013_ribo_14615 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Philippe2009; ./phylippart2multifasta.py -p ../../data_processed/matrices/Philippe2009 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Ryan2013_est; ./phylippart2multifasta.py -p ../../data_processed/matrices/Ryan2013_est -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Simion2017; ./phylippart2multifasta.py -p ../../data_processed/matrices/Simion2017 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Whelan2015_D10; ./phylippart2multifasta.py -p ../../data_processed/matrices/Whelan2015_D10 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Whelan2015_D1; ./phylippart2multifasta.py -p ../../data_processed/matrices/Whelan2015_D1 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Whelan2015_D20; ./phylippart2multifasta.py -p ../../data_processed/matrices/Whelan2015_D20 -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Whelan2017_full; ./phylippart2multifasta.py -p ../../data_processed/matrices/Whelan2017_full -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 echo ../../data_processed/matrices/Whelan2017_strict; ./phylippart2multifasta.py -p ../../data_processed/matrices/Whelan2017_strict -o ../blast/queries -t ../taxonomy_info/taxon_table.tsv
 cat ../blast/queries/*.fa > ../blast/db/animal_root.fa
