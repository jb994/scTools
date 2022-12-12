import numpy as np
import pandas as pd


def replaceClusters(new, old, adata, key='leiden'):
	if type(old)==str:
		newClusters = adata.obs[key].map({x:(new if x==old else x) for x in adata.obs[key]})
	elif type(old)==list:
		newClusters = adata.obs[key].map({x:(new if x in old else x) for x in adata.obs[key]})
	adata.obs[key] = newClusters

def largeScaleReplace(adata, key='coarseClusters'):
	## Just in case I need to replace a whole bunch again
	monoCells=['Monocytes', 'Monocytes1', 'Monocytes2', 'Monocytes4', 'Monocytes3', 'Macrophages', 
           'Macrophage/Monocyte', 'Kupffer Cells', 'Mono/Macrophages', 'Macrophages2']
	lymphCells = ['T Cells', 'T Cells1', 'T Cell2', 'T Cells2', 'y-delta T Cells', 
              'Y-Delta T Cells','Pre-Dysfunctional T Cell?', 'Naive T Cells']
	endoCells = ['Endothelial Cells1', 'Endothelial Cells2', 'Endothelial Cells3']
	bCells = ['B Cells', 'B Cells2']

    adata.obs[key] = adata.obs['leiden']
    process.replaceClusters('Monocytes', monoCells, adata, key=key)
    process.replaceClusters('T Cells', lymphCells, adata, key=key)
    process.replaceClusters('Endothelial Cells', endoCells, adata, key=key)
    process.replaceClusters('B Cells', bCells, adata, key=key)
    print(np.unique(adata.obs[key]))