import numpy as np
import pandas as pd
from collections import Counter
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import scanpy as sc
import scipy
import plotly.graph_objects as go
import networkx as nx
from scanpy.pl import MatrixPlot
from scanpy.plotting._dotplot import DotPlot
from scTools import process

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
		print("Error: Must define baseTreatment & contrastTreatment")
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
			cellFoldChanges = np.divide(contrastTemp.X, meanOfBase)
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
	im = heatmap_ax.imshow(heatmapMatrix, aspect='auto', norm=norm, cmap=cmap)
	heatmap_ax.set_title("Log Fold Changes")
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

def DotPlotter(adata, var_names, groupby, 
				depth='fineClusters', 
				ref=None,
				refKey=None,
				vmin=0, 
				vmax=3, 
				rotation=None,
				n_genes=10,
				cmap='bwr',
				size_df = None,
				color_df=None,
				categories_order=None,
				title=None,
				save=False):
	### Created because when using the typical sc.pl.dotplot function, it always scales the size of the\
	### dots in relation to %presence of a gene across ALL cell types and not the cell type in question.
	### The size of dots and color now relate to that gene's presence WITHIN the given celltype.

	if size_df is None or color_df is None:
		size_df, color_df = sizeNColor(adata, var_names, groupby, ref=ref, n_genes=n_genes, obsKey=depth)

	dp = DotPlot(adata,
		var_names,
		groupby,
		cmap=cmap,
		categories_order=categories_order,
		vmin=vmin, vmax=vmax,
		var_group_rotation=rotation,
		dot_color_df=color_df,
		dot_size_df=size_df,
		title=title,
		)
	plt.tight_layout()
	dp.show()
	if save is not False:
		dp.fig.savefig(save, bbox_inches='tight')


def matrixPlotter(adata, var_names, groupby, 
				depth='fineClusters', 
				ref=None,
				refKey=None,
				vmin=0, 
				vmax=3, 
				rotation=None,
				n_genes=10,
				cmap='bwr',
				color_df = None,
				categories_order=None,
				title=None,
				save=False):
	### Created because when using the typical sc.pl.dotplot function, it always scales the size of the\
	### dots in relation to %presence of a gene across ALL cell types and not the cell type in question.
	### The size of dots and color now relate to that gene's presence WITHIN the given celltype.

	if color_df==None:
		size_df, color_df = sizeNColor(adata, var_names, groupby, ref=ref, n_genes=n_genes, obsKey=depth)

	dp = MatrixPlot(adata,
		var_names,
		groupby,
		cmap=cmap,
		categories_order=categories_order,
		vmin=vmin, vmax=vmax,
		var_group_rotation=rotation,
		values_df=color_df,
		title=title,
		).style(cmap=cmap, edge_color='none')
	plt.tight_layout()
	dp.show()
	if save is not False:
		dp.fig.savefig(save, bbox_inches='tight')

def sizeNColor(adata, var_names, groupby='batch', ref=None, n_genes=5, obsKey='fineClusters'):
	groups = np.unique(adata.obs[groupby])
	sizeDF = pd.DataFrame(index=groups)
	colorDF = pd.DataFrame(index=groups)
	for cellType, genes in var_names.items():
		if len(genes)>0:
			tempAdata = adata[adata.obs[obsKey]==cellType][:,genes]
			sc.tl.rank_genes_groups(tempAdata, groupby, reference=ref, method='wilcoxon', key_added = "wilcoxon", n_genes=n_genes)

			#pctCellsWGene = [np.squeeze(np.asarray(tempAdata[tempAdata.obs[groupby]==t].X.todense().astype(bool).mean(axis=0))) for t in groups]
			#sizeTempDF = pd.DataFrame(pctCellsWGene, index=groups, columns=genes)
			sizeTempDF = getSizeDF(tempAdata, groups, groupby, genes)
			sizeDF = pd.concat([sizeDF, sizeTempDF], axis=1, join='inner').fillna(0)

			#readsPerCell = [list(sc.get.rank_genes_groups_df(tempAdata, group=t, key='wilcoxon').logfoldchanges) if (t in np.unique(tempAdata.obs[groupby])) and (t!=ref) else [0]*len(genes) for t in groups]
			#colorTempDF = pd.DataFrame(readsPerCell, index=groups, columns=genes)
			colorTempDF = getColorDF(tempAdata, groups, groupby, genes, ref=ref)
			colorDF = pd.concat([colorDF, colorTempDF], axis=1, join='inner').fillna(0)

	return sizeDF, colorDF


def getSizeDF(adata, groups, groupby='batch', cols=None):
	pctCellsWGene = [np.squeeze(np.asarray(adata[adata.obs[groupby]==t].X.todense().astype(bool).mean(axis=0))) for t in groups]
	return pd.DataFrame(pctCellsWGene, index=groups, columns=cols)

def getColorDF(adata, groups, groupby='batch', cols=None, ref=False):
	readsPerCell = [list(sc.get.rank_genes_groups_df(adata, group=t, key='wilcoxon').logfoldchanges) if (t in np.unique(adata.obs[groupby])) and (t!=ref) else [0]*len(cols) for t in groups]
	return pd.DataFrame(readsPerCell, index=groups, columns=cols)

def vennDiagram(df, treats, cellType, doFilter=False, title=None, save=None, returnCircles=False):
	from matplotlib_venn import venn2, venn3
	circles = []
	labels = []

	for t1, t2 in treats:
		tempdf = df[(df.treatment==f'{t1} x {t2}') & (df.cellType==cellType)]
		if doFilter: tempdf = process.filterDEGs(tempdf, t1, t2) #Filter p-val<=0.05 |lfc|>0.75
		if tempdf.shape[0]!=0:
			circles.append(set(tempdf['genes']))
			labels.append(f'{t1} vs {t2}')

	if len(circles)==3 or len(circles)==2:
		if len(circles)==3:
			venn3(circles, set_labels = labels)
		elif len(circles)==2:
			venn2(circles, set_labels = labels)
		plt.title(title)
		if save:
			plt.savefig(save)
		plt.show()
		
	else:
		print(f"No DEGs found for {cellType}")

	if returnCircles:
		return circles



def geneNetwork(df, nodeVariable='genes', coexpressionVariable='cellType'):
	### Make a network showing the connection of DEGs across multiple cellTypes
	### df requires columns to match 'nodes' and 'coexpression'

	graph = nx.Graph()
	nodes, edges, attr = getNodesAndEdges(df, nodeVariable, coexpressionVariable)
	addNodes(graph, nodes)
	addEdges(graph, edges)
	edge_trace, node_trace = getTraces(graph, attr)

	return createNetworkFig(edge_trace, node_trace)
	
def getNodesAndEdges(df, nodeVariable='genes', coexpressionVariable='cellType'):
	# Collect Nodes and Edges
	nodes = {}
	edges = {}
	attr = {}
	for gene in np.unique(df[nodeVariable]):
		uniqueCellTypes = np.unique(df[df[nodeVariable]==gene][coexpressionVariable])
		counts = Counter(df[df[coexpressionVariable].isin(uniqueCellTypes)][nodeVariable])
		nodes[gene]=len(uniqueCellTypes)
		edges[gene]=counts
		attr[gene]={'Cell Types':uniqueCellTypes, 'Connections':counts}
		#print(uniqueCellTypes)
		#attr[gene] = ('Cell Types', uniqueCellTypes)
		del edges[gene][gene]
	return nodes, edges, attr

def addNodes(graph, nodes):
	print("Adding Nodes")
	for node in nodes.keys():
		if nodes[node] > 0:
			graph.add_node(node, size = nodes[node])

def addEdges(graph, edges):
	print('Adding Edges')
	for node in edges.keys():
		for co_node in edges[node].keys():
			if edges[node][co_node] > 0:
				graph.add_edge(node, co_node, weight = edges[node][co_node])

def make_edge(x, y, text, width):
	return  go.Scatter(x         = x,
						y         = y,
						line      = dict(width = width,
							color = 'cornflowerblue'),
						hoverinfo = 'text',
						text      = ([text]),
						mode      = 'lines')

def getTraces(graph, attr=None):
	pos_ = nx.spring_layout(graph)
	edge_trace = makeEdgeTraces(graph, pos_)
	node_trace = makeNodeTraces(graph, pos_, attr)
	return edge_trace, node_trace

def makeEdgeTraces(graph, pos_):
	edge_trace = []
	for edge in graph.edges():
		if graph.edges()[edge]['weight'] > 0:
			char_1 = edge[0]
			char_2 = edge[1]
			x0, y0 = pos_[char_1]
			x1, y1 = pos_[char_2]
			text   = char_1 + '--' + char_2 + ': ' + str(graph.edges()[edge]['weight'])
	        
			trace  = make_edge([x0, x1, None], [y0, y1, None], text, 
								width = 0.3*graph.edges()[edge]['weight']**1.75)
			edge_trace.append(trace)
	return edge_trace

def makeNodeTraces(graph, pos_, attr=None):
	# Make a node trace
	print("Creating Node Traces")
	node_trace = go.Scatter(x         = [],
							y         = [],
							text      = [],
							textposition = "top center",
							textfont_size = 10,
							mode      = 'markers+text',
							hoverinfo = 'text',
							hovertext = [],
							marker    = dict(color = [],
							size  = [],
							line  = None))
	# For each node in gene connection Graph, get the position and size and add to the node_trace
	print("Adding Nodes Sizes")
	for node in graph.nodes():
		x, y = pos_[node]
		node_trace['x'] += tuple([x])
		node_trace['y'] += tuple([y])
		node_trace['marker']['color'] += tuple(['cornflowerblue'])
		node_trace['marker']['size'] += tuple([5*graph.nodes()[node]['size']])
		node_trace['text'] += tuple(['<b>' + node + '</b>'])
		if attr is not None: 
			celltypes = ', '.join(attr[node]['Cell Types'])
			connections=[]
			for i,(k,v) in enumerate(attr[node]['Connections'].most_common()):
				toappend = f'{k}: {v}'
				if i!=len(attr[node]['Connections'])-1:
					toappend += ', '
				if i%5==4:
					toappend += '<br>'
				connections.append(toappend)
			connections=''.join(connections)

			node_trace['hovertext'] += tuple( 
			[ f"<b>Gene</b>: {node}<br><b>Cell Types</b>: {celltypes}<br><b>Connections</b>: {connections}" ] 
			)
	return node_trace

def createNetworkFig(edge_trace, node_trace):
	### Layout
	layout = go.Layout(
		paper_bgcolor='rgba(0,0,0,0)', # transparent background
		plot_bgcolor='rgba(0,0,0,0)', # transparent 2nd background
		xaxis =  {'showgrid': False, 'zeroline': False}, # no gridlines
		yaxis = {'showgrid': False, 'zeroline': False}, # no gridlines
	)
	fig = go.Figure(layout=layout)
	print("Add all the Edge Traces to Figure")
	for trace in edge_trace:
	    fig.add_trace(trace)
	print("Adding Node Trace")
	fig.add_trace(node_trace)
	fig.update_layout(showlegend = False) # Remove legend
	fig.update_xaxes(showticklabels = False) # Remove tick labels
	fig.update_yaxes(showticklabels = False)
	return fig


def makeVolcanoPlot(df, 
	effect_size='logfoldchanges', 
	p='pvals', 
	gene='genes',
	snp=None,
	xlabel=None,
	ylabel=None,
	annotation='cellType',
	effect_size_line=[-1,1],
	genomewideline_value=-np.log10(10e-4),
	title = 'Volcano Plot'
	):

	if xlabel is None and effect_size=='logfoldchanges':
		xlabel='Log2 Fold Changes'

	import dash_bio as dashbio
	fig=dashbio.VolcanoPlot(dataframe=df, 
                    effect_size=effect_size, 
                    p=p,
                    gene=gene, 
                    snp=None,
                    xlabel=xlabel,
                    ylabel=None,
                    annotation=annotation,
                    effect_size_line=effect_size_line,
                    genomewideline_value=genomewideline_value,
                   )
	fig.update_layout(title_text=title, title_x=0.5)
	return fig


