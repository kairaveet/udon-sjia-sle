library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(xlsx)
source("~/FishersTest_ForClinicalFeatures.R")

## read the UDON groups file in 

library(readxl)
UDON_SLE_k25_21c <- read_excel("/Users/tha8tf/Downloads/UDON-SLE_k25-21c.xlsx")

groups_udon_sle = as.data.frame(UDON_SLE_k25_21c)

## read clinical data in 

ClinicalMetadata <- read_excel("/Users/tha8tf/Downloads/ClinicalMetadata.xlsx")

clinical_data = as.data.frame(ClinicalMetadata)

## some of the data here is not binary so we need to make binary 

idx = clinical_data == "NA"
clinical_data[idx] = NA

numeric_clinical_data = clinical_data[,8:54]
rownames(numeric_clinical_data) = clinical_data$SampleID

## clean up the numeric clinical data to make it binary based on average values of clinical parameters

numeric_clinical_data[ , ] <- apply(numeric_clinical_data[,], 2, function(x) as.numeric(as.character(x)))
colmeans_clinical_data = colMeans(numeric_clinical_data,na.rm = T)
numeric_clinical_data_og = numeric_clinical_data

for (c in 1:ncol(numeric_clinical_data)){
  
  numeric_clinical_data[,c] = numeric_clinical_data[,c] >= colmeans_clinical_data[c]
  numeric_clinical_data[,c] = 1*numeric_clinical_data[,c]
}

groups_udon_sle[,7:53] = numeric_clinical_data[groups_udon_sle$SampleID,]
groups_udon_sle$`UDON cluster (ICGS2 name)` = paste0("ICGS_Cluster_",groups_udon_sle$`UDON cluster (ICGS2 name)`)
groups_udon_sle[,5] = groups_udon_sle$`Cluster Name (SJIA PBMC aligned)`
colnames(groups_udon_sle)[44] = "DSDNA_Caps"
groups_udon_sle$age = gsub("SLE.*","",groups_udon_sle$SampleID)
print(unique(groups_udon_sle$age))


numeric_clinical_data_allcovs = numeric_clinical_data

cols_to_remove = names(which(colSums(numeric_clinical_data_allcovs) < 4))

groups_udon_sle_og = groups_udon_sle

groups_udon_sle = groups_udon_sle[,!names(groups_udon_sle) %in% cols_to_remove]

## let's deal with the categorical data here like race ethnicity etc
rownames(clinical_data) = clinical_data$SampleID
groups_udon_sle[,55:57] = clinical_data[groups_udon_sle$SampleID,5:7]

groups_udon_sle$Female = (groups_udon_sle$Gender == "F")*1
groups_udon_sle$Male = (groups_udon_sle$Gender == "M")*1

groups_udon_sle$AfricanAmerican = (groups_udon_sle$Race == "AA")*1
groups_udon_sle$Caucasian = (groups_udon_sle$Race == "C")*1
groups_udon_sle$Hispanic_R = (groups_udon_sle$Race == "H")*1
groups_udon_sle$Asian = (groups_udon_sle$Race == "As")*1

groups_udon_sle$Hispanic_E = (groups_udon_sle$Ethnicity == "H")*1
groups_udon_sle$NonHispanic = (groups_udon_sle$Ethnicity == "NH")*1

groups_udon_sle = groups_udon_sle[,-(55:57)]

vec = 7:62

p_vec = c()
p_list <- vector(mode = "list", length = length(vec))
names(p_list) = colnames(groups_udon_sle)[vec]

rowname_vec = c("B_Intermediate_c27","CD4_Cytotoxic_T_cell","Naive_CD8_Cytotoxic_T_cell","CD16_Mono",                 
                "Mono_c15", "CD14_Mono", "Macrophages", "preDCs", "Eryth", "Platelet", "CD4_T_Cells",
                "Plasmablast", "B_Intermediate", "Lymphoid_Cell_Cycle", "Platelet_Meg", "NK_CD56_Bright",          
                "HSPC", "MAIT", "CD8_TEM", "NK", "pDCs", "B_Naive", "B_Memory", "dNT", "CD4_TCM", "T_Reg",
                "CD8_T_Cells", "CD8_Naive","CD4_Naive")

for (clinc in vec){
  
  print(colnames(groups_udon_sle)[clinc])
  
  final_list = fisherstest_forclinicalfeats_cmh(groups_df = groups_udon_sle,clinicaltype = "1",col_number = clinc,p_val = 0.1,number_of_samples = 2,batch_col = "age")
  
  pm = final_list$pval
  pm = pm[rowname_vec,]
  #if (ncol(as.data.frame(pm)) != 21){next}
  
  #if (ncol(as.data.frame(pm)) != 21){next}
  
  p_lgm = pm < 0.1
  p_lgm = p_lgm*1
  
  orm = final_list$OR
  orm[orm=="Inf"] = NA
  orm = (as.matrix(orm))*1
  
  
  p_vec = c(p_vec,pm)
  p_list[[colnames(groups_udon_sle)[clinc]]] = as.data.frame(pm)
  
  

  if (sd(as.matrix(p_lgm),na.rm = T) == 0 || is.na(sd(as.matrix(p_lgm),na.rm = T))) {
    plot1=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = paste0(colnames(groups_udon_sle)[clinc],"_PVal_0.05"),breaks = c(0,1))
  } else {
    plot1=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = paste0(colnames(groups_udon_sle)[clinc],"_PVal_0.05"))
  }


  if (sd(as.matrix(orm),na.rm = T) == 0 || is.na(sd(as.matrix(orm),na.rm = T))){
    plot2=pheatmap::pheatmap(as.matrix(orm),color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = paste0(colnames(groups_udon_sle)[clinc],"_ORVal"),breaks = c(0,1))
  }
  else{
    plot2=pheatmap::pheatmap(as.matrix(orm),color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = paste0(colnames(groups_udon_sle)[clinc],"_ORVal"))
  }

  plot_list=list()
  plot_list[['plot1']]=plot1[[4]]
  plot_list[['plot2']]=plot2[[4]]

  pdf(paste0("SLE_CMH_PVal_0.1__post_filter_for_fdr_lessthan2", as.character(colnames(groups_udon_sle)[clinc]), ".pdf"),width = 15,height = 8)
  x = grid.arrange(grobs=plot_list, ncol=2,nrow=1)
  show(x)

  dev.off()
  
}

for (i in names(p_list)){
  
  print(i)
  
  ym = reshape2::melt(as.matrix(p_list[[i]]))
  ym$celltype = rownames(p_list[[i]])
  
  df_i_p_adj = p.adjust(ym$value,method="BH")
  print("BH:")
  idx2 = which(df_i_p_adj < 0.1)
  print(ym[idx2,])
}


vec = 7:62

p_vec = c()
p_list <- vector(mode = "list", length = length(vec))
names(p_list) = colnames(groups_udon_sle)[vec]

rowname_vec = c("B_Intermediate_c27","CD4_Cytotoxic_T_cell","Naive_CD8_Cytotoxic_T_cell","CD16_Mono",                 
                "Mono_c15", "CD14_Mono", "Macrophages", "preDCs", "Eryth", "Platelet", "CD4_T_Cells",
                "Plasmablast", "B_Intermediate", "Lymphoid_Cell_Cycle", "Platelet_Meg", "NK_CD56_Bright",          
                "HSPC", "MAIT", "CD8_TEM", "NK", "pDCs", "B_Naive", "B_Memory", "dNT", "CD4_TCM", "T_Reg",
                "CD8_T_Cells", "CD8_Naive","CD4_Naive")

for (clinc in vec){
  
  print(colnames(groups_udon_sle)[clinc])
  
  final_list = fisherstest_forclinicalfeats_adult_specific_all(groups_df = groups_udon_sle,clinicaltype = "1",col_number = clinc,p_val = 0.1,number_of_samples = 2,batch_col = "age")
  
  pm = final_list$pval
  pm = pm[rowname_vec,]
  #if (ncol(as.data.frame(pm)) != 21){next}
  
  #if (ncol(as.data.frame(pm)) != 21){next}
  
  p_lgm = pm < 0.1
  p_lgm = p_lgm*1
  
  orm = final_list$OR
  orm[orm=="Inf"] = NA
  orm = (as.matrix(orm))*1
  
  
  p_vec = c(p_vec,pm)
  p_list[[colnames(groups_udon_sle)[clinc]]] = as.data.frame(pm)
  
  
  
  if (sd(as.matrix(p_lgm),na.rm = T) == 0 || is.na(sd(as.matrix(p_lgm),na.rm = T))) {
    plot1=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = paste0(colnames(groups_udon_sle)[clinc],"_PVal_0.1"),breaks = c(0,1))
  } else {
    plot1=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = paste0(colnames(groups_udon_sle)[clinc],"_PVal_0.1"))
  }
  
  
  if (sd(as.matrix(orm),na.rm = T) == 0 || is.na(sd(as.matrix(orm),na.rm = T))){
    plot2=pheatmap::pheatmap(as.matrix(orm),color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = paste0(colnames(groups_udon_sle)[clinc],"_ORVal"),breaks = c(0,1))
  }
  else{
    plot2=pheatmap::pheatmap(as.matrix(orm),color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = paste0(colnames(groups_udon_sle)[clinc],"_ORVal"))
  }
  
  plot_list=list()
  plot_list[['plot1']]=plot1[[4]]
  plot_list[['plot2']]=plot2[[4]]
  
  pdf(paste0("SLE_ADULT_Specific_PVal_0.1_all_associations_lessthan2_samples_", as.character(colnames(groups_udon_sle)[clinc]), ".pdf"),width = 15,height = 8)
  x = grid.arrange(grobs=plot_list, ncol=2,nrow=1)
  show(x)
  
  dev.off()
  
}

p_list_adult = p_list
for (i in names(p_list_adult)){
  
  print(i)
  
  ym = reshape2::melt(as.matrix(p_list_adult[[i]]))
  ym$celltype = rownames(p_list_adult[[i]])
  
  df_i_p_adj = p.adjust(ym$value,method="BH")
  print("BH:")
  idx2 = which(df_i_p_adj < 0.1)
  print(ym[idx2,])
}



vec = 7:62

p_vec = c()
p_list <- vector(mode = "list", length = length(vec))
names(p_list) = colnames(groups_udon_sle)[vec]

rowname_vec = c("B_Intermediate_c27","CD4_Cytotoxic_T_cell","Naive_CD8_Cytotoxic_T_cell","CD16_Mono",                 
                "Mono_c15", "CD14_Mono", "Macrophages", "preDCs", "Eryth", "Platelet", "CD4_T_Cells",
                "Plasmablast", "B_Intermediate", "Lymphoid_Cell_Cycle", "Platelet_Meg", "NK_CD56_Bright",          
                "HSPC", "MAIT", "CD8_TEM", "NK", "pDCs", "B_Naive", "B_Memory", "dNT", "CD4_TCM", "T_Reg",
                "CD8_T_Cells", "CD8_Naive","CD4_Naive")

for (clinc in vec){
  
  print(colnames(groups_udon_sle)[clinc])
  
  final_list = fisherstest_forclinicalfeats_child_specific_all(groups_df = groups_udon_sle,clinicaltype = "1",col_number = clinc,p_val = 0.1,number_of_samples = 2,batch_col = "age")
  
  pm = final_list$pval
  pm = pm[rowname_vec,]
  #if (ncol(as.data.frame(pm)) != 21){next}
  
  #if (ncol(as.data.frame(pm)) != 21){next}
  
  p_lgm = pm < 0.1
  p_lgm = p_lgm*1
  
  orm = final_list$OR
  orm[orm=="Inf"] = NA
  orm = (as.matrix(orm))*1
  
  
  p_vec = c(p_vec,pm)
  p_list[[colnames(groups_udon_sle)[clinc]]] = as.data.frame(pm)
  
  
  
  if (sd(as.matrix(p_lgm),na.rm = T) == 0 || is.na(sd(as.matrix(p_lgm),na.rm = T))) {
    plot1=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = paste0(colnames(groups_udon_sle)[clinc],"_PVal_0.1"),breaks = c(0,1))
  } else {
    plot1=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = paste0(colnames(groups_udon_sle)[clinc],"_PVal_0.1"))
  }
  
  
  if (sd(as.matrix(orm),na.rm = T) == 0 || is.na(sd(as.matrix(orm),na.rm = T))){
    plot2=pheatmap::pheatmap(as.matrix(orm),color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = paste0(colnames(groups_udon_sle)[clinc],"_ORVal"),breaks = c(0,1))
  }
  else{
    plot2=pheatmap::pheatmap(as.matrix(orm),color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = paste0(colnames(groups_udon_sle)[clinc],"_ORVal"))
  }
  
  plot_list=list()
  plot_list[['plot1']]=plot1[[4]]
  plot_list[['plot2']]=plot2[[4]]
  
  pdf(paste0("SLE_CHILD_Specific_PVal_0.1_all_associations_lessthan2_samples_", as.character(colnames(groups_udon_sle)[clinc]), ".pdf"),width = 15,height = 8)
  x = grid.arrange(grobs=plot_list, ncol=2,nrow=1)
  show(x)
  
  dev.off()
  
}

p_list_child = p_list
for (i in names(p_list_child)){
  
  print(i)
  
  ym = reshape2::melt(as.matrix(p_list_child[[i]]))
  ym$celltype = rownames(p_list_child[[i]])
  
  df_i_p_adj = p.adjust(ym$value,method="BH")
  print("BH:")
  idx2 = which(df_i_p_adj < 0.1)
  print(ym[idx2,])
}

