library(Seurat)

setwd("/data/salomonis2/LabFiles/Kairavee/cna_evaluations_paper_revisions/")

so = readRDS("/data/salomonis2/LabFiles/Kairavee/Emely/PrepFilesforCNA/SeuratIntegratedObject_2NEW__subset_for_cna_withmetadata.rds")

so_counts = GetAssayData(object = so, slot = "counts")

cnaf = read.table("/data/salomonis2/LabFiles/Kairavee/cna_evaluations_paper_revisions/NAM_nbhdXpc.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F)


cor_counts = cor(t(as.matrix(so_counts)),cnaf,use= "pairwise.complete.obs",method = "pearson")

write.table(cor_counts,"cor_counts_cna_nam-pcs_genes.txt",sep = "\t",row.names = T,col.names = T,quote = F)

so_norm = GetAssayData(object = so, slot = "data")

cor_norm = cor(t(as.matrix(so_norm)),cnaf,use= "pairwise.complete.obs",method = "pearson")

write.table(cor_counts,"cor_norm_cna_nam-pcs_genes.txt",sep = "\t",row.names = T,col.names = T,quote = F)