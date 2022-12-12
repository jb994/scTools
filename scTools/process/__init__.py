import numpy as np
import pandas as pd


def replaceClusters(new, old, adata, key='leiden'):
	if type(old)==str:
		leiden = adata.obs['leiden'].map({x:(new if x==old else x) for x in adata.obs['leiden']})
	elif type(old)==list:
		leiden = adata.obs['leiden'].map({x:(new if x in old else x) for x in adata.obs['leiden']})
	adata.obs['leiden'] = leiden