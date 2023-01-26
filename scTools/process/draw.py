import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import scanpy as sc
from scanpy.plotting._dotplot import DotPlot
#from numpy import RuntimeWarning

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def heatmapper(adata, groupby='batch'):
	#sc.pl.heatmap(adata, marker_genes, groupby='louvain', figsize=(5, 8),
        #var_group_positions=[(0,1), (11, 12)], use_raw=False, vmin=-3, vmax=3, cmap='bwr',
        #var_group_labels=['B cells', 'dendritic'], var_group_rotation=0, dendrogram='dendrogram_louvain')

	sc.pl.rank_genes_groups_heatmap(adata, groupby=groupby, 
		use_raw=False, vmin=-3, vmax=3, cmap='bwr')


def DotPlotter(adata, var_names, groupby, 
				depth='fineClusters', 
				vmin=0, 
				vmax=3, 
				cmap='bwr',
				categories_order=None):
	### Created because when using the typical sc.pl.dotplot function, it always scales the size of the\
	### dots in relation to %presence of a gene across ALL cell types and not the cell type in question.
	### The size of dots and color now relate to that gene's presence WITHIN the given celltype.

	batches = np.unique(adata.obs['batch'])
	dot_color_df = pd.DataFrame(index=batches)
	dot_size_df = pd.DataFrame(index=batches)

	for cellType, genes in var_names.items():
		print(genes)
		tempAdata = adata[adata.obs[depth]==cellType][:,genes]
		sc.tl.rank_genes_groups(tempAdata, 'batch', method='wilcoxon', key_added = "wilcoxon")
		
		pctCellsWGene = [np.squeeze(np.asarray(tempAdata[tempAdata.obs['batch']==t].X.todense().astype(bool).mean(axis=0))) for t in batches]
		readsPerCell = [list(sc.get.rank_genes_groups_df(tempAdata, group=t, key='wilcoxon').logfoldchanges) if t in np.unique(tempAdata.obs['batch']) else [0]*len(genes) for t in batches]

		sizeTempDF = pd.DataFrame(pctCellsWGene, index=batches, columns=genes)
		colorTempDF = pd.DataFrame(readsPerCell, index=batches, columns=genes)
		colorTempDF=colorTempDF.apply(lambda x: x-x.mean(), axis=0)

		dot_size_df = pd.concat([dot_size_df, sizeTempDF], axis=1, join='inner').fillna(0)
		dot_color_df = pd.concat([dot_color_df, colorTempDF], axis=1, join='inner').fillna(0)

	dp = DotPlot(adata,
		var_names,
		groupby,
		cmap=cmap,
		categories_order=categories_order,
		vmin=vmin, vmax=vmax,
		dot_color_df=dot_color_df,
		dot_size_df=dot_size_df,
		).show()

