import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import anndata

def do_subsample(adata, subsample, seed=None):
	if not subsample or adata.shape[0]<subsample:
		pass
	elif subsample<1:
		sc.pp.subsample(adata, fraction=subsample, n_obs=None, random_state=seed, copy=False)
	else:
		sc.pp.subsample(adata, fraction=None, n_obs=subsample, random_state=seed, copy=False)
	return adata

def loadSTAR(path, subsample=None, seed=None, make_unique=False, labels=None):
	### Note: Directory must have a 
	''' 
	matrix.mtx
	genes.tsv
	barcodes.tsv
	'''
	adata = sc.read_10x_mtx(
		path,
		cache=True,
        make_unique=make_unique)
	if labels is True:
		cClust = pd.read_csv(path+'coarseClusters.csv', index_col=0)
		adata.obs = pd.merge(adata.obs, cClust,left_index=True, right_index=True,how='left')
		fClust = pd.read_csv(path+'fineClusters.csv', index_col=0)
		adata.obs = pd.merge(adata.obs, fClust, left_index=True,  right_index=True, how='left')
		adata = adata[adata.obs['coarseClusters'].notna() & (adata.obs['fineClusters'].notna())]
	adata = do_subsample(adata, subsample=subsample, seed=seed)
	return adata

def loadAllFriedman(path):
	treatments = ['NonfatVehicle', 'FriedVehicle', 'FriedCRV431', 'FriedLANI', 'FriedCRVLAN']
	Data = {}
	for treatment in treatments:
		Data[treatment] = sc.read_h5ad(path + treatment + '.h5ad')
	return Data