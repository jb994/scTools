import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

def heatmapper(adata, groupby='batch'):
	#sc.pl.heatmap(adata, marker_genes, groupby='louvain', figsize=(5, 8),
        #var_group_positions=[(0,1), (11, 12)], use_raw=False, vmin=-3, vmax=3, cmap='bwr',
        #var_group_labels=['B cells', 'dendritic'], var_group_rotation=0, dendrogram='dendrogram_louvain')

	sc.pl.rank_genes_groups_heatmap(adata, groupby=groupby, 
		use_raw=False, vmin=-3, vmax=3, cmap='bwr')