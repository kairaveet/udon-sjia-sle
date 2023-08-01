## see cell type frequencies in cna results 

## real disease (mas, lung disease, and active considered as "disease")

setwd("/Users/tha8tf/cna_results_paper_revisions/association_test_results_batch_effects/")

clusters = read.table("/Users/tha8tf/cluster_numer_to_type_translation_sjia.txt",sep = "\t",stringsAsFactors = F,header = T)

file_list = c("neighborhood_coefficents_Disease_Status_num",
              "neighborhood_coefficents_il1_inhibition",
              "neighborhood_coefficents_il18",
              "neighborhood_coefficents_lung_disease",
              "neighborhood_coefficents_real_disease",
              "neighborhood_coefficents_steroids",
              "neighborhood_coefficents_systemic_features",
              "neighborhood_coefficents_udon_pheno_complement_u6",
              "neighborhood_coefficents_udon_pheno_apoptotic_u11",
              "neighborhood_coefficents_udon_pheno_ifn_u12")

for (i in file_list){
  cna_ncorrs_i = read.table(paste0(i,".txt"), sep = "\t",stringsAsFactors = F, header = T)
  
  cna_ncorrs_i$cell_type = NA
  
  for (clus in unique(cna_ncorrs_i$seurat_clusters)){
    
    cna_ncorrs_i[cna_ncorrs_i$seurat_clusters == clus,9] = clusters[clusters$Cell.type.cluster==clus,2]
  }
  
  cluster_table = as.data.frame(table(cna_ncorrs_i$cell_type))
  rownames(cluster_table) = cluster_table$Var1
  
  expanded_population = cna_ncorrs_i[cna_ncorrs_i$poscells == "True",]
  depleted_population = cna_ncorrs_i[cna_ncorrs_i$negcells == "True",]
  
  
  composition_expanded_population = as.data.frame(table(expanded_population$cell_type))
  composition_depleted_population = as.data.frame(table(depleted_population$cell_type))
  
  rownames(composition_expanded_population) = composition_expanded_population$Var1
  composition_expanded_population$n_cells = cluster_table[rownames(composition_expanded_population),2]
  composition_expanded_population$perc_per_population = (composition_expanded_population$Freq/composition_expanded_population$n_cells)*100
  composition_expanded_population$Var1 = as.character(composition_expanded_population$Var1)
  composition_expanded_population = composition_expanded_population[order(composition_expanded_population$perc_per_population), ]
  composition_expanded_population$Var1 = factor(composition_expanded_population$Var1, levels = composition_expanded_population$Var1)
  
  rownames(composition_depleted_population) = composition_depleted_population$Var1
  composition_depleted_population$n_cells = cluster_table[rownames(composition_depleted_population),2]
  composition_depleted_population$perc_per_population = (composition_depleted_population$Freq/composition_depleted_population$n_cells)*100
  composition_depleted_population$Var1 = as.character(composition_depleted_population$Var1)
  composition_depleted_population = composition_depleted_population[order(composition_depleted_population$perc_per_population), ]
  composition_depleted_population$Var1 = factor(composition_depleted_population$Var1, levels = composition_depleted_population$Var1)
  
  
  write.table(composition_expanded_population, paste0(i,"_composition_expanded_population.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  write.table(composition_depleted_population, paste0(i,"_composition_depleted_population.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  
  pdf(paste0(i,"_bar_plot_perc_per_celltype_expanded_populations.pdf"), width = 8, height = 10)
    p = ggplot(composition_expanded_population, aes(x = Var1, y = perc_per_population)) + 
      geom_bar(stat="identity", width=0.7, fill = "red") + theme_classic() + coord_flip()
    show(p)
    dev.off()

  
  pdf(paste0(i,"_bar_plot_perc_per_celltype_depleted_populations.pdf"), width = 8, height = 10)
    p = ggplot(composition_depleted_population, aes(x = Var1, y = perc_per_population)) + 
      geom_bar(stat="identity", width=0.7, fill = "blue") + theme_classic() + coord_flip()
    show(p)
    dev.off()
  
}




