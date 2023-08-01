
udon_groups = read.table("/Users/tha8tf/udon_groups_k15_sjia.txt",sep = "\t",stringsAsFactors = F,header = T)

u4 = udon_groups[udon_groups$Pseudocluster == "10",]
u6 = udon_groups[udon_groups$Pseudocluster == "6",]

u4$PID = paste0(u4$PID, "-fold")
u6$PID = paste0(u6$PID, "-fold")

cna_cells_metadata$u4 = 0

for (i in unique(u4$PID)){
  
  u4_pid = u4[u4$PID == i, ]
  
  for (j in unique(u4_pid$Celltype)){
    
    cna_cells_metadata[(cna_cells_metadata$Patient_ID == i & cna_cells_metadata$seurat_clus_names == j), 11] = 1
  }
  
}


cna_cells_metadata$u6 = 0

for (i in unique(u6$PID)){
  
  u6_pid = u6[u6$PID == i, ]
  
  for (j in unique(u6_pid$Celltype)){
    
    cna_cells_metadata[(cna_cells_metadata$Patient_ID == i & cna_cells_metadata$seurat_clus_names == j), 12] = 1
  }
  
}

u7 = udon_groups[udon_groups$Pseudocluster == "8",]
u11 = udon_groups[udon_groups$Pseudocluster == "9",]
u12 = udon_groups[udon_groups$Pseudocluster == "1",]

u7$PID = paste0(u7$PID, "-fold")
u11$PID = paste0(u11$PID, "-fold")
u12$PID = paste0(u12$PID, "-fold")

cna_cells_metadata$u7 = 0

for (i in unique(u7$PID)){
  
  u7_pid = u7[u7$PID == i, ]
  
  for (j in unique(u7_pid$Celltype)){
    
    cna_cells_metadata[(cna_cells_metadata$Patient_ID == i & cna_cells_metadata$seurat_clus_names == j), 13] = 1
  }
  
}


cna_cells_metadata$u11 = 0

for (i in unique(u11$PID)){
  
  u11_pid = u11[u11$PID == i, ]
  
  for (j in unique(u11_pid$Celltype)){
    
    cna_cells_metadata[(cna_cells_metadata$Patient_ID == i & cna_cells_metadata$seurat_clus_names == j), 14] = 1
  }
  
}


cna_cells_metadata$u12 = 0

for (i in unique(u12$PID)){
  
  u12_pid = u12[u12$PID == i, ]
  
  for (j in unique(u12_pid$Celltype)){
    
    cna_cells_metadata[(cna_cells_metadata$Patient_ID == i & cna_cells_metadata$seurat_clus_names == j), 15] = 1
  }
  
}


