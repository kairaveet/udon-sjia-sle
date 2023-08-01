
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import cna
np.random.seed(0)

import anndata as ad
d = ad.read_h5ad('/data/salomonis2/LabFiles/Kairavee/cna_results/SeuratIntegratedObject_2NEW__subset_for_cna_withmetadata.h5ad')

print(d)
print(d.obs.id.head())
print(d.obs.Batch.head())

from multianndata import MultiAnnData
d = MultiAnnData(X=d.X, obs =d.obs, sampleid = "Patient_ID")
d.obs_to_sample(['Disease_Status_num','Batch'])
print(d.samplem.head())

##### add other metadata to the sample metadata

sample_metadata = pd.read_csv("samplem_metadata_4_all_metadata.txt",sep='\t')

sample_metadata.set_index("Patient_ID",inplace=True)

d.samplem['active'] = sample_metadata['Active']
d.samplem['inactive'] = sample_metadata['Inactive']
d.samplem['lung_disease'] = sample_metadata['Lung_Disease']
d.samplem['mas'] = sample_metadata['MAS']
d.samplem['real_disease'] = sample_metadata['Real_Disease']
d.samplem['il1_inhibition']	= sample_metadata['IL_1_inhibition']
d.samplem['il6r_inhibition'] = sample_metadata['IL6R_inhibition']
d.samplem['steroids'] = sample_metadata['Steriods']
d.samplem['ferritin'] = sample_metadata['Ferritin']
d.samplem['cxcl9'] = sample_metadata['CXCL9']
d.samplem['il18'] = sample_metadata['IL_18']
d.samplem['crp'] = sample_metadata['CRP']
d.samplem['s100a8_a9'] = sample_metadata['S100A8_A9']
d.samplem['s100a12'] = sample_metadata['S100A12']
d.samplem['anc_1000x'] = sample_metadata['ANC_1000x']
d.samplem['systemic_features'] = sample_metadata['Systemic_Features']
d.samplem['arthritis'] = sample_metadata['Active_arthritis']
d.samplem['fever'] = sample_metadata['Fever']
			 												
# d.samplem['active'] = sample_metadata['Active']
# d.samplem['inactive'] = sample_metadata['Inactive']
# d.samplem['lung_disease'] = sample_metadata['Lung_Disease']
# d.samplem['mas'] = sample_metadata['MAS']
# d.samplem['real_disease'] = sample_metadata['Real_Disease']


d.samplem.to_csv("samplem_metadata_5.txt",sep='\t')

seurat_clus = d.obs['seurat_clusters'].to_numpy()

d.obs['seurat_clusters'] = pd.Categorical(seurat_clus)

# compute the UMAP cell-cell similarity graph
sc.pp.neighbors(d)
print("Ran UMAP Cell by Cell")

# compute UMAP coordinates for plotting
sc.tl.umap(d)
print("Found UMAP Coords")


# ################################################################### Test for clinical covariates #######################################################################################
# 
sample_metadata_df = d.samplem
clinical_vars = sample_metadata_df.columns
clinical_vars = clinical_vars.drop('Batch')

for col in clinical_vars:
	print("testing print")
	print(col)
	
	# perform association test for case/ctrl status, controlling for sex as a covariate and accounting for potential batch effect
	res = cna.tl.association(d,                   #dataset
            d.samplem[col],                   #sample-level attribute of intest (case/control status,       #covariates to control for (in this case just one)
            batches=None)        #batch assignments for each sample so that cna can account for batch effects

	print('\nglobal association p-value:', res.p)

	# visualize which populations are expanded or depleted among case samples relative to cntrls
	plt.figure(figsize=(20, 20))
	cna.pl.umap_ncorr(d,                           #dataset
            res,                               #cna results object
            scatter0={'alpha':0.5, 's':20},    #plotting parameters for neighborhoods that pass FDR
            scatter1={'alpha':0.05, 's':20})   #plotting parameters for neighborhoods that don't pass FDR
	plt.title('p = {:.2e}'.format(res.p))
	plt.savefig(("neighborhood_coefficient_no_batch_effects_subset_" + col + ".pdf"))
	plt.show()

	# Global association test p-value
	print('Global association test p-value: ', res.p, ',', res.k, 'PCs used')

	# Correlation threshold for 5% FDR
	print('Correlation threshold for 5% FDR', res.fdr_5p_t)

	# Correlation threshold for 10% FDR
	print('Correlation threshold for 10% FDR', res.fdr_10p_t)

	try:
		# Number of neighborhoods with local associations at 5% FDR
		n = np.sum(abs(res.ncorrs) > res.fdr_5p_t)
		print('Number of neighborhoods with local associations', n)
	
	except:
		print("no associations found for this clinical covariate")
		continue
	
	FDR_thresh = res.fdr_5p_t
	
	d.obs['ncorrs'] = res.ncorrs
	
	# Positively-associated cells
	d.obs['poscells'] = np.repeat(False, d.obs.shape[0])
	d.obs['poscells'].loc[d.obs['ncorrs']>FDR_thresh] = True
	
	# Negatively-associated cells
	d.obs['negcells'] = np.repeat(False, d.obs.shape[0])
	d.obs['negcells'].loc[d.obs['ncorrs']<-FDR_thresh] = True
	
	to_export = pd.DataFrame(d.obs[['id', 'Disease_Status', 'seurat_clusters', 'Batch', 'ncorrs','poscells', 'negcells']])
	to_export.to_csv(("neighborhood_coefficents_no_batch_" + col + ".txt"), sep='\t')
	
	# d.obs = d.obs.drop(labels = ['poscells','negcells','ncorrs'])
	
###############################################################################################################################################################

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


## write files 

d.uns['NAM_nbhdXpc'].to_csv('NAM_nbhdXpc.txt',sep='\t')

d_kept_cells = d[d.uns['keptcells'],:]



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











