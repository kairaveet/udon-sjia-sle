fisherstest_forclinicalfeats <- function(groups_df,clinicaltype,col_number){
  
  celltypes_unique = unique(groups_df$V5)
  diseasesubtypes_unique = unique(groups_df$V4)
  cluster_unique = unique(groups_df$V3)
  
  # loop through each cell type and each subtype 
  disease_subtype = clinicaltype
  
  pval_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(pval_matrix) = celltypes_unique
  colnames(pval_matrix) = cluster_unique
  
  oddsratio_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(oddsratio_matrix) = celltypes_unique
  colnames(oddsratio_matrix) = cluster_unique
  
  for (cluster in cluster_unique){
    
    for (celltype in celltypes_unique){
      
      ## pseudobulks from that cluster
      pseudobulks_cluster = groups_df[groups_df$V3 == cluster,]
      pseudobulks_other_cluster = groups_df[!groups_df$V3 == cluster,]
      
      # celltype_disease_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 == disease_subtype & pseudobulks_cluster$V5 == celltype,])
      # notcelltype_disease_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 == disease_subtype & pseudobulks_cluster$V5 != celltype,])
      # 
      # notdisease_celltype_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 != disease_subtype & pseudobulks_cluster$V5 == celltype,])
      # notdisease_notcelltype_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 != disease_subtype & pseudobulks_cluster$V5 != celltype,])
      
      diseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 == disease_subtype & pseudobulks_cluster$V5 == celltype,])
      notdiseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 != disease_subtype & pseudobulks_cluster$V5 == celltype,])
      
      diseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster$V4 == disease_subtype & pseudobulks_other_cluster$V5 == celltype,])
      notdiseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster$V4 != disease_subtype & pseudobulks_other_cluster$V5 == celltype,])
      
      fisher_df = data.frame("in ClusterX" = c(diseasetype_celltype_centroids_in_cluster,notdiseasetype_celltype_centroids_in_cluster), "Not_in_ClusterX" = c(diseasetype_celltype_centroids_in_othercluster,notdiseasetype_celltype_centroids_in_othercluster), row.names = c("Disease_CellType", "NotDisease_CellType"))
      
      test_results_df = fisher.test(as.matrix(fisher_df))
      
      pval_matrix[celltype,cluster] = test_results_df$p.value
      oddsratio_matrix[celltype,cluster] = test_results_df$estimate[["odds ratio"]]
      
    }
  }
  
  final_list = list(pval = pval_matrix, OR = oddsratio_matrix)
  return(final_list)
  
}


fisherstest_forclinicalfeats <- function(groups_df,clinicaltype,col_number,p_val,number_of_samples){
  na_idx = which(is.na(groups_df[,col_number]))
  if (length(na_idx) > 0 ){
  groups_df = groups_df[-na_idx,]
  }
  celltypes_unique = unique(groups_df[,5])
  diseasesubtypes_unique = unique(groups_df[,col_number])
  cluster_unique = unique(groups_df[,3])
  # loop through each cell type and each subtype
  disease_subtype = clinicaltype
  pval_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(pval_matrix) = celltypes_unique
  colnames(pval_matrix) = cluster_unique
  oddsratio_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(oddsratio_matrix) = celltypes_unique
  colnames(oddsratio_matrix) = cluster_unique
  for (cluster in cluster_unique){
    for (celltype in celltypes_unique){
      ## pseudobulks from that cluster
      pseudobulks_cluster = groups_df[groups_df[,3]==cluster,]
      pseudobulks_other_cluster = groups_df[!groups_df[,3] == cluster,]
      # celltype_disease_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 == disease_subtype & pseudobulks_cluster$V5 == celltype,])
      # notcelltype_disease_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 == disease_subtype & pseudobulks_cluster$V5 != celltype,])
      #
      # notdisease_celltype_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 != disease_subtype & pseudobulks_cluster$V5 == celltype,])
      # notdisease_notcelltype_pseudobulks = nrow(pseudobulks_cluster[pseudobulks_cluster$V4 != disease_subtype & pseudobulks_cluster$V5 != celltype,])
      diseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number] == disease_subtype & pseudobulks_cluster[,5] == celltype,])
      notdiseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number]  != disease_subtype & pseudobulks_cluster[,5] == celltype,])
      diseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  == disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
      notdiseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  != disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
      fisher_df = data.frame("in ClusterX" = c(diseasetype_celltype_centroids_in_cluster,notdiseasetype_celltype_centroids_in_cluster), "Not_in_ClusterX" = c(diseasetype_celltype_centroids_in_othercluster,notdiseasetype_celltype_centroids_in_othercluster), row.names = c("Disease_CellType", "NotDisease_CellType"))
      test_results_df = fisher.test(as.matrix(fisher_df),alternative = "greater")
      pval_matrix[celltype,cluster] = test_results_df$p.value
      oddsratio_matrix[celltype,cluster] = test_results_df$estimate[["odds ratio"]]
      # if (cluster == "ICGS_Cluster_20"){
      #   if (celltype == "NK_CD56_Bright"){
      #     print(fisher_df)
      #   }
      # }
      
      if (test_results_df$p.value < p_val){
        if (fisher_df[1,1] < number_of_samples){
          pval_matrix[celltype,cluster] = 1 ## put NA here rather than call it 1
          #pval_matrix[celltype,cluster] = NA
          #print("Less than 3 samples situation occured")
        }
      }
      
    }
  }
  final_list = list(pval = pval_matrix, OR = oddsratio_matrix)
  return(final_list)
}


#if (test_results_df$p.value < p_val){
if (fisher_df[1,1] < number_of_samples){
  #pval_matrix[celltype,cluster] = 1 ## put NA here rather than call it 1
  pval_matrix[celltype,cluster] = NA
  print("Less than 3 samples situation occured")
}
#}


fisherstest_forclinicalfeats_cmh <- function(groups_df,clinicaltype,col_number,p_val,number_of_samples,batch_col){
  na_idx = which(is.na(groups_df[,col_number]))
  if (length(na_idx) > 0 ){
    groups_df = groups_df[-na_idx,]
  }
  
  celltypes_unique = unique(groups_df[,5])
  diseasesubtypes_unique = unique(groups_df[,col_number])
  cluster_unique = unique(groups_df[,3])
  batch_unique = unique(groups_df[,batch_col])
  # loop through each cell type and each subtype
  disease_subtype = clinicaltype
  pval_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(pval_matrix) = celltypes_unique
  colnames(pval_matrix) = cluster_unique
  oddsratio_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(oddsratio_matrix) = celltypes_unique
  colnames(oddsratio_matrix) = cluster_unique
  
  groups_df_og = groups_df

  
  for (cluster in cluster_unique){
    for (celltype in celltypes_unique){
      ## initiate an empty 3-d array to be filled out
      print(cluster)
      print(celltype)
      contigency_list <- vector(mode = "list", length = length(batch_unique))
      names(contigency_list) = batch_unique

      for (batch_component_idx in 1:length(batch_unique)){
        groups_df = groups_df_og
        batch_component = batch_unique[batch_component_idx]
        print(batch_component)
        groups_df = groups_df[groups_df[,batch_col] == batch_component,]
        ## pseudobulks from that cluster
        pseudobulks_cluster = groups_df[groups_df[,3]==cluster,]
        pseudobulks_other_cluster = groups_df[!groups_df[,3] == cluster,]
        diseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number] == disease_subtype & pseudobulks_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number]  != disease_subtype & pseudobulks_cluster[,5] == celltype,])
        diseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  == disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  != disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        contigency_df = data.frame("in ClusterX" = c(diseasetype_celltype_centroids_in_cluster,notdiseasetype_celltype_centroids_in_cluster), "Not_in_ClusterX" = c(diseasetype_celltype_centroids_in_othercluster,notdiseasetype_celltype_centroids_in_othercluster), row.names = c("Disease_CellType", "NotDisease_CellType"))
        
        ## put this df as one of the "stacking pieces" in the 3d array 
        contigency_list[[batch_component]] = as.matrix(contigency_df)
      }
        #contigency_array <- array(c(contigency_list$a,contigency_list$c), dim=c(2, 2, length(batch_unique)))
        contigency_array <- array(unlist(contigency_list), dim=c(2, 2, length(batch_unique)))
        
        if (length(contigency_list) == 1){
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        } else if (any(contigency_array[1,1 ,] < number_of_samples)) {
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        } else {
          test_results_df <- tryCatch(
            {test_results_df = mantelhaen.test(contigency_array,alternative = "greater")
            }, error = function(e){
              print(e)
              # low_sample_size_idx = which(apply(contigency_array, 3L, sum) < 2) ## the array which has sample size/ total # of samples in contigency df less than 1 
              # print(paste("doing fisher's test for",cluster, "and",celltype, "within",batch_unique[low_sample_size_idx]))
              # contigency_array = contigency_array[, , -low_sample_size_idx]
              # test_results_df = fisher.test(contigency_array,alternative = "greater")
              test_results_df = vector(mode = "list", length = 2)
              names(test_results_df) = c("p.value","estimate")
              test_results_df$p.value = NA
              test_results_df$estimate = NA
              return(test_results_df)
            }
          ) 
        }
        
        
        pval_matrix[celltype,cluster] = test_results_df$p.value
        
        odds_ratio <- tryCatch(
          {odds_ratio = test_results_df$estimate[["common odds ratio"]]
          }, error = function(e){
            print(e)
            odds_ratio = test_results_df$estimate
            return(odds_ratio)
          }
        ) 
          
        oddsratio_matrix[celltype,cluster] = odds_ratio
    }
  }
  
  final_list = list(pval = pval_matrix, OR = oddsratio_matrix)
  return(final_list)
}


#if (test_results_df$p.value < p_val){
if (fisher_df[1,1] < number_of_samples){
  #pval_matrix[celltype,cluster] = 1 ## put NA here rather than call it 1
  pval_matrix[celltype,cluster] = NA
  print("Less than 3 samples situation occured")
}
#}



fisherstest_forclinicalfeats_adult_specific <- function(groups_df,clinicaltype,col_number,p_val,number_of_samples,batch_col){
  na_idx = which(is.na(groups_df[,col_number]))
  if (length(na_idx) > 0 ){
    groups_df = groups_df[-na_idx,]
  }
  
  celltypes_unique = unique(groups_df[,5])
  diseasesubtypes_unique = unique(groups_df[,col_number])
  cluster_unique = unique(groups_df[,3])
  batch_unique = unique(groups_df[,batch_col])
  # loop through each cell type and each subtype
  disease_subtype = clinicaltype
  pval_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(pval_matrix) = celltypes_unique
  colnames(pval_matrix) = cluster_unique
  oddsratio_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(oddsratio_matrix) = celltypes_unique
  colnames(oddsratio_matrix) = cluster_unique
  
  groups_df_og = groups_df
  
  
  for (cluster in cluster_unique){
    for (celltype in celltypes_unique){
      ## initiate an empty 3-d array to be filled out
      print(cluster)
      print(celltype)
      contigency_list <- vector(mode = "list", length = length(batch_unique))
      names(contigency_list) = batch_unique
      
      for (batch_component_idx in 1:length(batch_unique)){
        groups_df = groups_df_og
        batch_component = batch_unique[batch_component_idx]
        print(batch_component)
        groups_df = groups_df[groups_df[,batch_col] == batch_component,]
        ## pseudobulks from that cluster
        pseudobulks_cluster = groups_df[groups_df[,3]==cluster,]
        pseudobulks_other_cluster = groups_df[!groups_df[,3] == cluster,]
        diseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number] == disease_subtype & pseudobulks_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number]  != disease_subtype & pseudobulks_cluster[,5] == celltype,])
        diseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  == disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  != disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        contigency_df = data.frame("in ClusterX" = c(diseasetype_celltype_centroids_in_cluster,notdiseasetype_celltype_centroids_in_cluster), "Not_in_ClusterX" = c(diseasetype_celltype_centroids_in_othercluster,notdiseasetype_celltype_centroids_in_othercluster), row.names = c("Disease_CellType", "NotDisease_CellType"))
        
        ## put this df as one of the "stacking pieces" in the 3d array 
        contigency_list[[batch_component]] = as.matrix(contigency_df)
      }
      #contigency_array <- array(c(contigency_list$a,contigency_list$c), dim=c(2, 2, length(batch_unique)))
      contigency_array <- array(unlist(contigency_list), dim=c(2, 2, length(batch_unique)))
      
      if (length(contigency_list) == 1){
        ## check if the only available matrix is adult or child 
        if ("a" %in% names(contigency_list)){
          print(paste("doing fisher's test for",cluster, "and",celltype))
          contigency_array = contigency_array[, , 1]
          test_results_df = fisher.test(contigency_array,alternative = "greater")
        } else {
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        }
        
      } else if (contigency_array[1,1,2] < number_of_samples) { ## if the number of samples in kids much less then check if there are enough for adults
        if (contigency_array[1,1,1] < number_of_samples){
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        } else {
          print(paste("doing fisher's test for",cluster, "and",celltype))
          contigency_array = contigency_array[, , 1]
          test_results_df = fisher.test(contigency_array,alternative = "greater")
        }
      } else {
            # low_sample_size_idx = which(apply(contigency_array, 3L, sum) < 2) ## the array which has sample size/ total # of samples in contigency df less than 1 
            # print(paste("doing fisher's test for",cluster, "and",celltype, "within",batch_unique[low_sample_size_idx]))
            # contigency_array = contigency_array[, , -low_sample_size_idx]
            # test_results_df = fisher.test(contigency_array,alternative = "greater")
            test_results_df = vector(mode = "list", length = 2)
            names(test_results_df) = c("p.value","estimate")
            test_results_df$p.value = NA
            test_results_df$estimate = NA
          }
      
      
      pval_matrix[celltype,cluster] = test_results_df$p.value
      
      odds_ratio <- tryCatch(
        {odds_ratio = test_results_df$estimate[["odds ratio"]]
        }, error = function(e){
          print(e)
          odds_ratio = test_results_df$estimate
          return(odds_ratio)
        }
      ) 
      
      oddsratio_matrix[celltype,cluster] = odds_ratio
    }
  }
  
  final_list = list(pval = pval_matrix, OR = oddsratio_matrix)
  return(final_list)
}





fisherstest_forclinicalfeats_child_specific <- function(groups_df,clinicaltype,col_number,p_val,number_of_samples,batch_col){
  na_idx = which(is.na(groups_df[,col_number]))
  if (length(na_idx) > 0 ){
    groups_df = groups_df[-na_idx,]
  }
  
  celltypes_unique = unique(groups_df[,5])
  diseasesubtypes_unique = unique(groups_df[,col_number])
  cluster_unique = unique(groups_df[,3])
  batch_unique = unique(groups_df[,batch_col])
  # loop through each cell type and each subtype
  disease_subtype = clinicaltype
  pval_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(pval_matrix) = celltypes_unique
  colnames(pval_matrix) = cluster_unique
  oddsratio_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(oddsratio_matrix) = celltypes_unique
  colnames(oddsratio_matrix) = cluster_unique
  
  groups_df_og = groups_df
  
  
  for (cluster in cluster_unique){
    for (celltype in celltypes_unique){
      ## initiate an empty 3-d array to be filled out
      print(cluster)
      print(celltype)
      contigency_list <- vector(mode = "list", length = length(batch_unique))
      names(contigency_list) = batch_unique
      
      for (batch_component_idx in 1:length(batch_unique)){
        groups_df = groups_df_og
        batch_component = batch_unique[batch_component_idx]
        print(batch_component)
        groups_df = groups_df[groups_df[,batch_col] == batch_component,]
        ## pseudobulks from that cluster
        pseudobulks_cluster = groups_df[groups_df[,3]==cluster,]
        pseudobulks_other_cluster = groups_df[!groups_df[,3] == cluster,]
        diseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number] == disease_subtype & pseudobulks_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number]  != disease_subtype & pseudobulks_cluster[,5] == celltype,])
        diseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  == disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  != disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        contigency_df = data.frame("in ClusterX" = c(diseasetype_celltype_centroids_in_cluster,notdiseasetype_celltype_centroids_in_cluster), "Not_in_ClusterX" = c(diseasetype_celltype_centroids_in_othercluster,notdiseasetype_celltype_centroids_in_othercluster), row.names = c("Disease_CellType", "NotDisease_CellType"))
        
        ## put this df as one of the "stacking pieces" in the 3d array 
        contigency_list[[batch_component]] = as.matrix(contigency_df)
      }
      #contigency_array <- array(c(contigency_list$a,contigency_list$c), dim=c(2, 2, length(batch_unique)))
      contigency_array <- array(unlist(contigency_list), dim=c(2, 2, length(batch_unique)))
      
      if (length(contigency_list) == 1){
        ## check if the only available matrix is adult or child 
        if ("c" %in% names(contigency_list)){
          print(paste("doing fisher's test for",cluster, "and",celltype))
          contigency_array = contigency_array[, , 1]
          test_results_df = fisher.test(contigency_array,alternative = "greater")
        } else {
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        }
        
      } else if (contigency_array[1,1,1] < number_of_samples) { ## if the number of samples in kids much less then check if there are enough for adults
        if (contigency_array[1,1,2] < number_of_samples){
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        } else {
          print(paste("doing fisher's test for",cluster, "and",celltype))
          contigency_array = contigency_array[, , 1]
          test_results_df = fisher.test(contigency_array,alternative = "greater")
        }
      } else {
        # low_sample_size_idx = which(apply(contigency_array, 3L, sum) < 2) ## the array which has sample size/ total # of samples in contigency df less than 1 
        # print(paste("doing fisher's test for",cluster, "and",celltype, "within",batch_unique[low_sample_size_idx]))
        # contigency_array = contigency_array[, , -low_sample_size_idx]
        # test_results_df = fisher.test(contigency_array,alternative = "greater")
        test_results_df = vector(mode = "list", length = 2)
        names(test_results_df) = c("p.value","estimate")
        test_results_df$p.value = NA
        test_results_df$estimate = NA
      }
      
      
      pval_matrix[celltype,cluster] = test_results_df$p.value
      
      odds_ratio <- tryCatch(
        {odds_ratio = test_results_df$estimate[["odds ratio"]]
        }, error = function(e){
          print(e)
          odds_ratio = test_results_df$estimate
          return(odds_ratio)
        }
      ) 
      
      oddsratio_matrix[celltype,cluster] = odds_ratio
    }
  }
  
  final_list = list(pval = pval_matrix, OR = oddsratio_matrix)
  return(final_list)
}




fisherstest_forclinicalfeats_adult_specific_all <- function(groups_df,clinicaltype,col_number,p_val,number_of_samples,batch_col){
  na_idx = which(is.na(groups_df[,col_number]))
  if (length(na_idx) > 0 ){
    groups_df = groups_df[-na_idx,]
  }
  
  celltypes_unique = unique(groups_df[,5])
  diseasesubtypes_unique = unique(groups_df[,col_number])
  cluster_unique = unique(groups_df[,3])
  batch_unique = sort(unique(groups_df[,batch_col]))
  # loop through each cell type and each subtype
  disease_subtype = clinicaltype
  pval_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(pval_matrix) = celltypes_unique
  colnames(pval_matrix) = cluster_unique
  oddsratio_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(oddsratio_matrix) = celltypes_unique
  colnames(oddsratio_matrix) = cluster_unique
  
  groups_df_og = groups_df
  
  
  for (cluster in cluster_unique){
    for (celltype in celltypes_unique){
      ## initiate an empty 3-d array to be filled out
      print(cluster)
      print(celltype)
      contigency_list <- vector(mode = "list", length = length(batch_unique))
      names(contigency_list) = batch_unique
      
      for (batch_component_idx in 1:length(batch_unique)){
        groups_df = groups_df_og
        batch_component = batch_unique[batch_component_idx]
        print(batch_component)
        groups_df = groups_df[groups_df[,batch_col] == batch_component,]
        ## pseudobulks from that cluster
        pseudobulks_cluster = groups_df[groups_df[,3]==cluster,]
        pseudobulks_other_cluster = groups_df[!groups_df[,3] == cluster,]
        diseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number] == disease_subtype & pseudobulks_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number]  != disease_subtype & pseudobulks_cluster[,5] == celltype,])
        diseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  == disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  != disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        contigency_df = data.frame("in ClusterX" = c(diseasetype_celltype_centroids_in_cluster,notdiseasetype_celltype_centroids_in_cluster), "Not_in_ClusterX" = c(diseasetype_celltype_centroids_in_othercluster,notdiseasetype_celltype_centroids_in_othercluster), row.names = c("Disease_CellType", "NotDisease_CellType"))
        
        ## put this df as one of the "stacking pieces" in the 3d array 
        contigency_list[[batch_component]] = as.matrix(contigency_df)
      }
      #contigency_array <- array(c(contigency_list$a,contigency_list$c), dim=c(2, 2, length(batch_unique)))
      contigency_array <- array(unlist(contigency_list), dim=c(2, 2, length(batch_unique)))
      
      if ("a" %in% names(contigency_list)){
        if (contigency_array[1,1,1] < number_of_samples){
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        } else {
          print(paste("doing fisher's test for",cluster, "and",celltype))
          contigency_array = contigency_array[, , 1]
          test_results_df = fisher.test(contigency_array,alternative = "greater")
        }
      } else {
          test_results_df = vector(mode = "list", length = 2)
          names(test_results_df) = c("p.value","estimate")
          test_results_df$p.value = NA
          test_results_df$estimate = NA
        }
      
      pval_matrix[celltype,cluster] = test_results_df$p.value
      
      odds_ratio <- tryCatch(
        {odds_ratio = test_results_df$estimate[["odds ratio"]]
        }, error = function(e){
          print(e)
          odds_ratio = test_results_df$estimate
          return(odds_ratio)
        }
      ) 
      
      oddsratio_matrix[celltype,cluster] = odds_ratio
    }
  }
  
  final_list = list(pval = pval_matrix, OR = oddsratio_matrix)
  return(final_list)
}




fisherstest_forclinicalfeats_child_specific_all <- function(groups_df,clinicaltype,col_number,p_val,number_of_samples,batch_col){
  na_idx = which(is.na(groups_df[,col_number]))
  if (length(na_idx) > 0 ){
    groups_df = groups_df[-na_idx,]
  }
  
  celltypes_unique = unique(groups_df[,5])
  diseasesubtypes_unique = unique(groups_df[,col_number])
  cluster_unique = unique(groups_df[,3])
  batch_unique = sort(unique(groups_df[,batch_col]))
  # loop through each cell type and each subtype
  disease_subtype = clinicaltype
  pval_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(pval_matrix) = celltypes_unique
  colnames(pval_matrix) = cluster_unique
  oddsratio_matrix = data.frame(matrix(ncol = length(cluster_unique), nrow = length(celltypes_unique)))
  rownames(oddsratio_matrix) = celltypes_unique
  colnames(oddsratio_matrix) = cluster_unique
  
  groups_df_og = groups_df
  
  
  for (cluster in cluster_unique){
    for (celltype in celltypes_unique){
      ## initiate an empty 3-d array to be filled out
      print(cluster)
      print(celltype)
      contigency_list <- vector(mode = "list", length = length(batch_unique))
      names(contigency_list) = batch_unique
      
      for (batch_component_idx in 1:length(batch_unique)){
        groups_df = groups_df_og
        batch_component = batch_unique[batch_component_idx]
        print(batch_component)
        groups_df = groups_df[groups_df[,batch_col] == batch_component,]
        ## pseudobulks from that cluster
        pseudobulks_cluster = groups_df[groups_df[,3]==cluster,]
        pseudobulks_other_cluster = groups_df[!groups_df[,3] == cluster,]
        diseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number] == disease_subtype & pseudobulks_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_cluster = nrow(pseudobulks_cluster[pseudobulks_cluster[,col_number]  != disease_subtype & pseudobulks_cluster[,5] == celltype,])
        diseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  == disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        notdiseasetype_celltype_centroids_in_othercluster = nrow(pseudobulks_other_cluster[pseudobulks_other_cluster[,col_number]  != disease_subtype & pseudobulks_other_cluster[,5] == celltype,])
        contigency_df = data.frame("in ClusterX" = c(diseasetype_celltype_centroids_in_cluster,notdiseasetype_celltype_centroids_in_cluster), "Not_in_ClusterX" = c(diseasetype_celltype_centroids_in_othercluster,notdiseasetype_celltype_centroids_in_othercluster), row.names = c("Disease_CellType", "NotDisease_CellType"))
        
        ## put this df as one of the "stacking pieces" in the 3d array 
        contigency_list[[batch_component]] = as.matrix(contigency_df)
      }
      #contigency_array <- array(c(contigency_list$a,contigency_list$c), dim=c(2, 2, length(batch_unique)))
      contigency_array <- array(unlist(contigency_list), dim=c(2, 2, length(batch_unique)))
      
      if ("c" %in% names(contigency_list)){
        if(length(contigency_list) == 1){
          if (contigency_array[1,1,1] < number_of_samples){
            test_results_df = vector(mode = "list", length = 2)
            names(test_results_df) = c("p.value","estimate")
            test_results_df$p.value = NA
            test_results_df$estimate = NA
          } else {
            print(paste("doing fisher's test for",cluster, "and",celltype))
            contigency_array = contigency_array[, , 1]
            test_results_df = fisher.test(contigency_array,alternative = "greater")
          }
        } else {
          if (contigency_array[1,1,2] < number_of_samples){
            test_results_df = vector(mode = "list", length = 2)
            names(test_results_df) = c("p.value","estimate")
            test_results_df$p.value = NA
            test_results_df$estimate = NA
          } else {
            print(paste("doing fisher's test for",cluster, "and",celltype))
            contigency_array = contigency_array[, , 2]
            test_results_df = fisher.test(contigency_array,alternative = "greater")
          }
        }
        
        
      } else {
        test_results_df = vector(mode = "list", length = 2)
        names(test_results_df) = c("p.value","estimate")
        test_results_df$p.value = NA
        test_results_df$estimate = NA
      }
      
      pval_matrix[celltype,cluster] = test_results_df$p.value
      
      odds_ratio <- tryCatch(
        {odds_ratio = test_results_df$estimate[["odds ratio"]]
        }, error = function(e){
          print(e)
          odds_ratio = test_results_df$estimate
          return(odds_ratio)
        }
      ) 
      
      oddsratio_matrix[celltype,cluster] = odds_ratio
    }
  }
  
  final_list = list(pval = pval_matrix, OR = oddsratio_matrix)
  return(final_list)
}











