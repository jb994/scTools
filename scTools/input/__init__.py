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

def loadSTAR(path, subsample=None, seed=None, make_unique=False):
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
	adata = do_subsample(adata, subsample=subsample, seed=seed)
	return adata
