
library(ggplot2)
require(gridExtra)
source( "functions.R" )

########################################################################################## load matrices

taxonomy_reference = read_tsv("../reconciliation/taxonomy_info/taxon_table.tsv")
clades = c( "Fungi", "Ichthyosporea", "Filasterea", "Choanoflagellida", "Ctenophora", 
						"Porifera", "Placozoa", "Bilateria", "Cnidaria" )
# raw columns:
# component_number, matrix, partition_name, edges, nodes_in_component, component_density, BUSCO_ID, BUSCO_description,
# SwissProt_accession, SwissProt_description, GO_annotations, ribo_found
partition_map_global =
	read_tsv("../reconciliation/blast/graphs/partition_components_split_annotated.tsv") %>% #component nr first occurance
	dplyr::rename(partition = partition_name) %>%
	mutate( component_number = as.character(component_number) )

taxa =
	taxonomy_reference %>%
	distinct( relabelled_name, clade_assignment, ncbi_tax_id ) %>%
	dplyr::rename( taxon = relabelled_name, clade = clade_assignment ) %>%
	mutate( clade = factor( clade, levels = clades )  )


matrix_path = "../data_processed/phy_files"
phylip_file_names = list.files(path = matrix_path, pattern = ".+\\.phy$", full.names = TRUE)

#sequence_matrices = foreach( phylip_file = phylip_file_names) %dopar% parse_phylip( phylip_file )   #Error messages?

sequence_matrices = lapply(phylip_file_names , parse_phylip)

# Matrix gene composition

busco_results =
	read_tsv("../reconciliation/blast/graphs/busco_metazoa_results.tsv") %>%
	filter( Status != "Missing" )
# Multiple fields are in one colon delimited string. Need to parse them out.
# Busco result example
# "Moroz2014:ED3a:Capitella:51293:0241"
# manuscript:matrix:species:NCBI_taxon_id:partition

Bs =
	str_split_fixed( busco_results$Sequence, ":", 4 ) %>%
	as_tibble()

names( Bs ) = c( "matrix", "species", "ncbi_taxon_id", "partition" )

Bs %<>% mutate( ncbi_taxon_id = as.integer(ncbi_taxon_id) )

busco_results %<>% bind_cols( Bs )

# Combine and summarize results
busco_distinct =  busco_results %>% select( matrix, partition, Description ) %>% distinct(  )

partition_to_busco_map =
	busco_distinct %>%
	group_by( matrix, partition ) %>%
	summarise( BUSCO = names(which(table(Description) == max(table(Description)))[1]) )

partition_map_global %<>% left_join( partition_to_busco_map, by = c("matrix", "partition") )

# Matrix overlap
matrix_overlap =
	lapply(sequence_matrices, function(x) lapply(sequence_matrices, function(y) compute_matrix_overlap(x, y))) %>%
	unlist(recursive = FALSE) %>%
	bind_rows()

overlap_rects =  
	lapply(sequence_matrices, function(x) lapply(sequence_matrices, function(y) overlap_rects(x, y))) %>%
	unlist(recursive = FALSE) %>%
	bind_rows()

########################################################################################## save data

#Genes missing component numbers
comp_nr_dunn = sequence_matrices[[1]]@partitions$component_number
missing_nr_dunn = comp_nr_dunn[str_detect(comp_nr_dunn,"Dunn")==TRUE]
missing_nr_dunn = substring(missing_nr_dunn, 10,17)

fileConn<-file("../files_mine/missing_nr_dunn2008.txt")
writeLines(missing_nr_dunn, fileConn)
close(fileConn)

comp_nr_phil = sequence_matrices[[2]]@partitions$component_number
missing_nr_phil = comp_nr_phil[str_detect(comp_nr_phil,"Philippe")==TRUE]
missing_nr_phil = substring(missing_nr_phil, 14,22)

fileConn<-file("../files_mine/missing_nr_philippe2009.txt")
writeLines(missing_nr_phil, fileConn)
close(fileConn)

########################################################################################## plot

df = overlap_rects


p_list <- list()
i=0

for (val in seq(1,nrow(df),2))
{
	df_int <- df[c(val,val+1),]
	title=paste(df_int$matrix_1[1],df_int$matrix_2[1])
	i = i+1
	#nam <- paste("plot", val, sep = "")
	
	#assign(nam, ggplot(df_int, aes_string(xmin = df_int$xmin, xmax =df_int$xmax, ymin = df_int$ymin, ymax =df_int$ymax )) +
				 	#geom_rect(fill = "pink", alpha = 0.4) )#+ theme_void() + theme(plot.background = element_rect(fill = "black"))
	p_list[[i]] <- ggplot(df_int, aes_string(xmin = df_int$xmin, xmax =df_int$xmax, ymin = df_int$ymin, ymax =df_int$ymax )) +
		geom_rect(fill = c(1,2), alpha = 0.4) +
		  ggtitle(title) + theme(plot.title = element_text(size = 4))
}

do.call(grid.arrange,p_list)


########################################################################################## skräp

#xmin_l <- df$xmin
#xmax_l <- df$xmax
#ymin_l <- df$ymin
#ymax_l <- df$ymax


#plot_all <- ggplot(df, aes_string(xmin = df$xmin, xmax =df$xmax, ymin = df$ymin, ymax =df$ymax )) +
#	geom_rect(color = NA, fill = "pink", alpha = 0.4) #+ theme_void() + theme(plot.background = element_rect(fill = "black"))

#df_int <- df[c(4,5),]

#plot_intersecting <- ggplot(df_int, aes_string(xmin = df_int$xmin, xmax =df_int$xmax, ymin = df_int$ymin, ymax =df_int$ymax )) +
#	geom_rect(fill = "pink", alpha = 0.4) #+ theme_void() + theme(plot.background = element_rect(fill = "black"))

#plot_em = ggplot(df)+geom_blank()

#grid.arrange(plot2, plot_em, plot4, plot7, ncol=floor(sqrt(nrow(df))), nrow=ceiling(sqrt(nrow(df))))

#lapply(dflist, function(df) {
	# Do some complex operations on each data frame, df
	# More steps
	
	# Make sure the last thing is NULL. The last statement within the function will be
	# returned to lapply, which will try to combine these as a list across all data frames.
	# You don't actually care about this, you just want to run the function.
	
	
#	NULL
#})

#grid.arrange(plot2, plot_em, plot4, plot7, ncol=2, nrow=2)

#plot + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))