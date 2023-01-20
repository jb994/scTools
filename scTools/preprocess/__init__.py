import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse.linalg
from scanpy.pp import filter_cells, filter_genes, normalize_total, log1p, highly_variable_genes

def fullFilterShebang(adata, 
				plot=True, 
				doFilter=True, 
				doNorm=True, 
				flavor='seurat', 
				doMito=False, 
				doBlood=False, 
				selectGenes=800):
	if doFilter: filter(adata)
	if doMito: adata = mitoFilter(adata, plot=plot)
	if doBlood: adata = bloodFilter(adata, plot=plot)
	if doNorm: logNorm(adata)
	adata = hvg(adata, flavor=flavor, plot=plot, selectGenes=selectGenes)
	adata.raw = adata
	print(f"After cell and gene filtering we have {adata.shape}")
	sc.pp.scale(adata, max_value=10)
	return adata


def fullEmbedShebang(adata, n_components=50, plot=True):
	doPCA(adata, n_components=n_components, plot=plot)
	umapify(adata, plot=plot)


def filter(adata,
		min_genes=200,
		min_cells=3):
	filter_cells(adata, min_genes=min_genes)
	filter_genes(adata, min_cells=min_cells)


def mitoFilter(adata,plot=True):
	adata.var['mt'] = adata.var_names.str.startswith('MT-')
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	if plot: sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
	adata = adata[adata.obs.pct_counts_mt < 5, :]
	adata = adata[adata.obs.n_genes_by_counts < 2500, :]
	return adata

def riboFilter(adata,plot=True):
	adata.var['ribo'] = adata.var_names.str.startswith(("Rps","Rpl"))
	sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
	if plot: sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)
	adata = adata[adata.obs.pct_counts_ribo > 5, :]
	adata = adata[adata.obs.n_genes_by_counts < 2500, :]
	return adata


def bloodFilter(adata, plot=True):
	adata.var['blood'] = adata.var_names.str.startswith('Hbb-' or 'Hba-')
	sc.pp.calculate_qc_metrics(adata, qc_vars=['blood'], percent_top=None, log1p=False, inplace=True)
	if plot: sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_blood'],
             jitter=0.4, multi_panel=True)
	adata = adata[adata.obs.pct_counts_blood <10 , :]
	return adata

def logNorm(adata):
	normalize_total(adata)
	log1p(adata)

def hvg(adata, flavor='seurat', selectGenes=1000, plot=True):
	if flavor!='seurat_v3':
		highly_variable_genes(adata, flavor=flavor, min_mean=0.0125, max_mean=3, min_disp=0.5)
	else:
		highly_variable_genes(adata, flavor=flavor, n_top_genes=selectGenes)
	if plot: sc.pl.highly_variable_genes(adata)
	adata = adata[:, adata.var.highly_variable]
	return adata

def doPCA(adata, n_components=50, plot=True, inplace=True):
	if 'highly_variable' not in adata.var: 
		highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=4000)
	sc.tl.pca(adata, n_comps=n_components, use_highly_variable=True)
	if plot: sc.pl.pca_variance_ratio(adata, n_pcs=n_components, log=True)
	if not inplace: return adata

def umapify(adata, plot=True, color=None, title='batch', legend_loc='right margin', redoPCA=False):
	if redoPCA: adata=doPCA(adata, inplace=False)
	sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
	sc.tl.umap(adata)
	if plot: sc.pl.umap(adata, color=color, title=title, legend_loc=legend_loc)


