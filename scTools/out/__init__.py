import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from scTools import process
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

import warnings
from anndata import ImplicitModificationWarning
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

def makeSpreadsheet(df, filename=None, sheetKey = 'cellType', sortKey=None ):
	### Given a df and a sheetKey, will create an excel spreadsheet where each sheet\
	### is for each unique value in the sheetKey column
	if filename is None:
		print("Must give valid filename location")
		quit()
	wb = Workbook()
	ws1 = wb.active
	ws1.title = "All"
	for r in dataframe_to_rows(df, index=False, header=True):
		ws1.append(r)

	for key in np.unique(df[sheetKey]):
		print(f"On {key}")
		keyDF = df[df[sheetKey]==key]
		if sortKey: keyDF = keyDF.sort_values(sortKey)
		ws = wb.create_sheet(title=key)
		for r in dataframe_to_rows(keyDF, index=False, header=True):
			ws.append(r)

	wb.save(filename)
	print(f"File {filename} created")


def rankCellTypes(Data, treatments, silent=False, heatmap=True, key='fineClusters'):
	### Put in a dictionary of AnnDatas and their keys
	### For each pair of AnnDatas a correlation score is determined for\
	### how much each of the celltypes have changed between groups.
	### Returns a DF of treatment comparisons and cellTypes
	cols=[]
	dfData=[]
	allTypes = np.unique(process.catAdata(Data, treatments).obs[key])

	for i,t1 in enumerate(treatments):
		for j,t2 in enumerate(treatments):
			if i<j:
				print(f"{t1} x {t2}")
				cols.append(f"{t1} x {t2}")
				bigData = process.catAdata(Data, [t1,t2])
				bigData = process.dgeNorm(bigData, useRaw=False, scale=False, pca=True)
				#sc.pp.scale(bigData)
				#sc.pp.normalize_per_cell(bigData, counts_per_cell_after=1e4)
				#sc.pp.log1p(bigData)
				#sc.tl.pca(bigData, n_comps=50)
				row=[]
				for cellType in allTypes:
					if not silent: print(cellType)
					if (cellType in np.unique(bigData[bigData.obs['batch']==t1].obs[key])) and (cellType in np.unique(bigData[bigData.obs['batch']==t2].obs[key])):
						adata = bigData[bigData.obs[key]==cellType]
						sc.pp.scale(adata)
						sc.tl.rank_genes_groups(adata, 'batch', method='wilcoxon', use_raw=False, max_iter=2000, pts=False,cor_method='spearman')
						sc.tl.dendrogram(adata, groupby='batch', )
						if not silent: print(adata.uns['dendrogram_batch']['correlation_matrix'][0,1])
						if heatmap: sc.pl.rank_genes_groups_heatmap(adata, vmin=-3, vmax=3, cmap='bwr')
						# mpl.colors.SymLogNorm or mpl.colors.CenteredNorm #norm=mpl.colors.CenteredNorm()
						row.append(adata.uns['dendrogram_batch']['correlation_matrix'][0,1])
					else:
						row.append(None)

				dfData.append(row)
	return pd.DataFrame(np.asarray(dfData).T, columns=cols, index=allTypes)





