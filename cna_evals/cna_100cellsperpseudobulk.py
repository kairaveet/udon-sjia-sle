
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import cna
np.random.seed(0)

import anndata as ad
d = ad.read_h5ad('/data/salomonis2/LabFiles/Kairavee/cna_results/SeuratIntegratedObject_2NEW__subset_for_cna_withmetadata.h5ad')

#subsetcells_df = pd.read_csv("/data/salomonis2/LabFiles/Kairavee/Emely/PrepFilesforCNA/metaData_AllCells__forSubsetCells",sep="\t")


print(d)
print(d.obs.id.head())
print(d.obs.Disease_Status_num.head())
print(d.obs.Batch.head())

from multianndata import MultiAnnData
d = MultiAnnData(X=d.X, obs =d.obs, sampleid = "Patient_ID")
d.obs_to_sample(['Disease_Status_num','Batch'])
#print(d)
print(d.samplem.head())

seurat_clus = d.obs['seurat_clusters'].to_numpy()

d.obs['seurat_clusters'] = pd.Categorical(seurat_clus)

# compute the UMAP cell-cell similarity graph
sc.pp.neighbors(d)
print("Ran UMAP Cell by Cell")

# compute UMAP coordinates for plotting
sc.tl.umap(d)
print("Found UMAP Coords")

 
plt.figure(figsize=(20, 20))
sc.pl.umap(d, color='seurat_clusters', title='Seurat Clusters')
plt.savefig("umaplot1_seuratclusters.pdf",bbox_inches='tight')
plt.show()
plt.close()

plt.figure(figsize=(20, 20))
sc.pl.umap(d, color='Patient_ID', title='Sample IDs')
plt.savefig("umaplot1_samples.pdf",bbox_inches='tight')
plt.show()
plt.close()

plt.figure(figsize=(20, 20))
sc.pl.umap(d, color='Disease_Status', title='Subtype Information')
plt.savefig("umaplot1_subtypeinfo.pdf",bbox_inches='tight')
plt.show()
plt.close()

## ********** INDEPENDENT of association testing, calculate NAM matrix and get the NAM-PCs **********
cna.tl.nam(d, batches=d.samplem.Batch)

d.uns['NAM_nbhdXpc'].to_csv('NAM_nbhdXpc.txt',sep='\t')

d_kept_cells = d[d.uns['keptcells'],:]

d_kept_cells.obs.to_csv('multianndata_OBS_metadata_sjia_after_cna.txt',sep='\t')

corr_mat = np.corrcoef(d.uns['NAM_nbhdXpc'],d.X)

np.savetxt("corr_matrix_nam_pcs_with_genes.txt",corr_mat, delimiter = "\t")


# plot variance explained by each NAM PC
plt.figure(figsize=(20, 20))
plt.plot(d.uns['NAM_svs']/d.uns['NAM_svs'].sum(), marker='o', linestyle='--')
plt.xlabel('NAM PC')
plt.ylabel('fraction of variance explained')
plt.savefig("variance_explained_sjia_nam-pc_no_covs.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 1
plt.figure(figsize=(20, 20))
d.obs['NAMPC1'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC1'] = d.uns['NAM_nbhdXpc'].PC1
sc.pl.umap(d, color='NAMPC1', cmap='seismic', title='NAM PC1 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC1.pdf",bbox_inches='tight')
plt.show()
plt.close()

#plot neighborhood loadings of NAM PC 2
plt.figure(figsize=(20, 20))
d.obs['NAMPC2'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC2'] = d.uns['NAM_nbhdXpc'].PC2
sc.pl.umap(d, color='NAMPC2', cmap='seismic', title='NAM PC2 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC2.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 3
plt.figure(figsize=(20, 20))
d.obs['NAMPC3'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC3'] = d.uns['NAM_nbhdXpc'].PC3
sc.pl.umap(d, color='NAMPC3', cmap='seismic', title='NAM PC3 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC3.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 4
plt.figure(figsize=(20, 20))
d.obs['NAMPC4'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC4'] = d.uns['NAM_nbhdXpc'].PC4
sc.pl.umap(d, color='NAMPC4', cmap='seismic', title='NAM PC4 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC4.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 5
plt.figure(figsize=(20, 20))
d.obs['NAMPC5'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC5'] = d.uns['NAM_nbhdXpc'].PC5
sc.pl.umap(d, color='NAMPC5', cmap='seismic', title='NAM PC5 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC5.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 6
plt.figure(figsize=(20, 20))
d.obs['NAMPC6'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC6'] = d.uns['NAM_nbhdXpc'].PC6
sc.pl.umap(d, color='NAMPC6', cmap='seismic', title='NAM PC6 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC6.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 7
plt.figure(figsize=(20, 20))
d.obs['NAMPC7'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC7'] = d.uns['NAM_nbhdXpc'].PC7
sc.pl.umap(d, color='NAMPC7', cmap='seismic', title='NAM PC7 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC7.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 8
plt.figure(figsize=(20, 20))
d.obs['NAMPC8'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC8'] = d.uns['NAM_nbhdXpc'].PC8
sc.pl.umap(d, color='NAMPC8', cmap='seismic', title='NAM PC8 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC8.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 9
plt.figure(figsize=(20, 20))
d.obs['NAMPC9'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC9'] = d.uns['NAM_nbhdXpc'].PC9
sc.pl.umap(d, color='NAMPC9', cmap='seismic', title='NAM PC9 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC9.pdf",bbox_inches='tight')
plt.show()
plt.close()


#plot neighborhood loadings of NAM PC 10
plt.figure(figsize=(20, 20))
d.obs['NAMPC10'] = np.nan
d.obs.loc[d.uns['keptcells'], 'NAMPC10'] = d.uns['NAM_nbhdXpc'].PC10
sc.pl.umap(d, color='NAMPC10', cmap='seismic', title='NAM PC10 neighborhood loadings')
plt.savefig("umaplot1_NAM-PC10.pdf",bbox_inches='tight')
plt.show()
plt.close()











