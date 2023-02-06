import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import scanpy as sc
import scipy
from scanpy.plotting._dotplot import DotPlot
#from numpy import RuntimeWarning

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def heatmapper(adata, baseTreatment=None, contrastTreatment=None, 
			var_names=None,
			cellTypes=None, 
			cellTypesKey='fineClusters', 
			groupby='batch'):
	### Make a heatmap which plots LogFoldChanges across multiple cell types
	### plt's ax.imshow() is useful for this

	if (baseTreatment is None) or (contrastTreatment is None):
		print("Error: Must define baseTreatment & contrastTreatmtn")
		quit()
	if cellTypes is None:
		cellTypes = np.unique(adata.obs[cellTypesKey])


	heatmapMatrix = []
	row_names = []
	baseData = adata[adata.obs['batch'] == baseTreatment]
	contrastData = adata[adata.obs['batch'] == contrastTreatment]
	genes = [item for v in var_names.values() for item in v]

	for cellType in cellTypes:
		baseTemp = baseData[baseData.obs[cellTypesKey]==cellType, genes]
		contrastTemp = contrastData[contrastData.obs[cellTypesKey]==cellType, genes]

		if scipy.sparse.issparse(baseTemp.X):
			meanOfBase = baseTemp.X.todense().mean(axis=0) #mean of every gene's expression
		else:
			meanOfBase = baseTemp.X.mean(axis=0)

		if scipy.sparse.issparse(contrastTemp.X):
			cellFoldChanges = np.divide(contrastTemp.X.todense(), meanOfBase, axis=1)
		else:
			#print(contrastTemp.X.shape)
			#print(len(meanOfBase))
			cellFoldChanges = np.divide(contrastTemp.X, meanOfBase)
			#print(cellFoldChanges.shape)
		heatmapMatrix.append(np.log1p(cellFoldChanges))

	heatmapMatrix = np.nan_to_num(np.vstack(heatmapMatrix).T, 0)
	print(heatmapMatrix)
	print(heatmapMatrix.shape)

	width=10
	height=len(var_names) * 0.18
	height_ratios = [0, height, 0]
	width_ratios = [width, 0, 0]
	fig = plt.figure(figsize=(width, height))
	norm = plt.Normalize(vmin=-3, vmax=3)
	cmap='bwr'
	axs = gridspec.GridSpec(
		nrows=3,
		ncols=3,
		width_ratios=width_ratios,
		wspace=0.25 / width,
		hspace=0.3 / height,
		height_ratios=height_ratios,
	)
	heatmap_ax = fig.add_subplot(axs[1, 0])
	#kwds.setdefault('interpolation', 'nearest')
	im = heatmap_ax.imshow(heatmapMatrix, aspect='auto', norm=norm, cmap=cmap)
	#heatmap_ax.set_xlim(0 - 0.5, heatmapMatrix.shape[0] - 0.5)
	#heatmap_ax.set_ylim(heatmapMatrix.shape[1] - 0.5, -0.5)
	'''
	heatmap_ax.tick_params(axis='x', bottom=False, labelbottom=False)
	heatmap_ax.set_xlabel('')
	heatmap_ax.grid(False)
	heatmap_ax.tick_params(axis='y', labelsize='small', length=1)
	heatmap_ax.set_yticks(np.arange(len(var_names)))
	heatmap_ax.set_yticklabels(var_names, rotation=0)
	'''
	#fig, ax = plt.subplots()
	#im = ax.imshow(heatmapMatrix)
	heatmap_ax.set_title("Log Fold Changes")
	#fig.tight_layout()
	plt.show()


	#sc.pl.rank_genes_groups_heatmap(adata, groupby=groupby, 
	#	use_raw=False, vmin=-3, vmax=3, cmap='bwr')


def snsHeatmap(df, var_names=None, vmax=3, vmin=-3, cmap='bwr',):
	### Makes a heatmap of log2foldchanges given a df
	if var_names is not None:
		df = df[:,var_names].copy()

	plt.figure(figsize = (5,16))
	ax=sns.heatmap(df,
		vmax=vmax,
		vmin=vmin,
		cmap=cmap,
		)
	#fig.set_size_inches(20, 20)
	#fig.ylabel('Genes')
	#fig.yticks([0.5,1.5])


def DotPlotter(adata, var_names, groupby, 
				depth='fineClusters', 
				ref=None,
				refKey=None,
				vmin=0, 
				vmax=3, 
				rotation=None,
				n_genes=10,
				cmap='bwr',
				categories_order=None,
				title=None,
				save=False):
	### Created because when using the typical sc.pl.dotplot function, it always scales the size of the\
	### dots in relation to %presence of a gene across ALL cell types and not the cell type in question.
	### The size of dots and color now relate to that gene's presence WITHIN the given celltype.

	batches = np.unique(adata.obs['batch'])
	dot_color_df = pd.DataFrame(index=batches)
	dot_size_df = pd.DataFrame(index=batches)

	for cellType, genes in var_names.items():
		#print(genes)
		tempAdata = adata[adata.obs[depth]==cellType][:,genes]
		sc.tl.rank_genes_groups(tempAdata, 'batch', reference=ref, method='wilcoxon', key_added = "wilcoxon", n_genes=n_genes)
		
		pctCellsWGene = [np.squeeze(np.asarray(tempAdata[tempAdata.obs['batch']==t].X.todense().astype(bool).mean(axis=0))) for t in batches]
		sizeTempDF = pd.DataFrame(pctCellsWGene, index=batches, columns=genes)

		readsPerCell = [list(sc.get.rank_genes_groups_df(tempAdata, group=t, key='wilcoxon').logfoldchanges) if (t in np.unique(tempAdata.obs['batch'])) and (t!=ref) else [0]*len(genes) for t in batches]
		colorTempDF = pd.DataFrame(readsPerCell, index=batches, columns=genes)
		#print("\nColorTempDF")
		#print(colorTempDF)
		'''
		if ref is None:
			colorTempDF=colorTempDF.apply(lambda x: x-x.mean(), axis=0)
		else:
			if refKey is not None:
				pass
			else:
				colorTempDF=colorTempDF.apply(lambda x: x-x[ref], axis=0)
		'''

		dot_size_df = pd.concat([dot_size_df, sizeTempDF], axis=1, join='inner').fillna(0)
		dot_color_df = pd.concat([dot_color_df, colorTempDF], axis=1, join='inner').fillna(0)

	dp = DotPlot(adata,
		var_names,
		groupby,
		cmap=cmap,
		categories_order=categories_order,
		vmin=vmin, vmax=vmax,
		var_group_rotation=rotation,
		dot_color_df=dot_color_df,
		dot_size_df=dot_size_df,
		title=title,
		)
	plt.tight_layout()
	dp.show()
	if save is not False:
		dp.fig.savefig(save, bbox_inches='tight')
	#dp.show()



