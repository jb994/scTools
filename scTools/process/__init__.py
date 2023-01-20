
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import gseapy
import copy

def replaceClusters(new, old, adata, key='leiden'):
    if type(old) == str:
        newClusters = adata.obs[key].map({x: new if x == old else x for x in adata.obs[key]})
    elif type(old) == list:
        newClusters = adata.obs[key].map({x: new if x in old else x for x in adata.obs[key]})
    adata.obs[key] = newClusters


def subClusters(adata, coarseLabel, fineLabel, fineData=None, coarseKey='coarseClusters', fineKey='fineClusters'):
    if fineData is None:
        fineData = adata
    newClusters = []
    for cell in adata.obs.index:
        if cell in fineData.obs.index and fineData.obs[fineKey][cell] == fineLabel:
            print(f"Changing {adata.obs[coarseKey][cell]} for {coarseLabel}")
            newClusters.append(coarseLabel)
        else:
            newClusters.append(adata.obs[coarseKey][cell])
    return newClusters.copy()


def largeScaleReplace(adata, key='coarseClusters'):
    monoCells = [
     "'Monocytes'", "'Monocytes1'", "'Monocytes2'", "'Monocytes4'", 
     "'Monocytes3'", "'Macrophages'", 
     "'Macrophage/Monocyte'",  "'Kupffer Cells'",  "'Mono/Macrophages'",  "'Macrophages2'"]
    lymphCells = ["'T Cells'", "'T Cells1'", "'T Cell2'", "'T Cells2'", "'y-delta T Cells'", 
     "'Y-Delta T Cells'",  "'Pre-Dysfunctional T Cell?'",  "'Naive T Cells'"]
    endoCells = ['Endothelial Cells1', 'Endothelial Cells2', 'Endothelial Cells3']
    bCells = ['B Cells', 'B Cells2']
    adata.obs[key] = adata.obs['leiden']
    process.replaceClusters('Monocytes', monoCells, adata, key=key)
    process.replaceClusters('T Cells', lymphCells, adata, key=key)
    process.replaceClusters('Endothelial Cells', endoCells, adata, key=key)
    process.replaceClusters('B Cells', bCells, adata, key=key)
    print(np.unique(adata.obs[key]))


def batchIt(Data, keys=None, inplace=True):
    if not keys:
        keys = Data.keys()
    if not inplace:
        Data = copy.deepcopy(Data)
    for i, key in enumerate(keys):
        Data[key].obs['batch'] = i
    if not inplace:
        return Data


def batchCorrect(Data):
    batchIt(Data, inplace=True)
    adata = ad.concat([Data[key] for key in Data.keys()])
    sc.pp.combat(adata)


def getSizeFactors(adata):
    adata.X = np.expm1(adata.X)
    counts = adata[np.all(adata, axis=1)]
    logcounts = np.log1p(counts)
    logGeoMeans = np.mean(logcounts, axis=1).reshape(len(logcounts), 1)
    sf = np.expm1(np.median((logcounts - loggeommeans), axis=0))
    return sf


def equalizeReadMeans(adata):
    print('Equalizing Batch Mean Read Counts')
    readMeans = []
    adata.X = np.expm1(adata.X)
    for batch in np.unique(adata.obs['batch']):
        readMeans.append(np.sum(adata[adata.obs['batch'] == batch].obs['total_counts']) / adata[adata.obs['batch'] == batch].shape[0])

    avgReadMean = np.mean(readMeans)
    for i, batch in enumerate(np.unique(adata.obs['batch'])):
        adata[adata.obs['batch'] == batch].X = adata[adata.obs['batch'] == batch].X * (avgReadMean / readMeans[i])

    adata.X = np.log1p(adata.X)
    return adata


def catAdata(adataDict, keys, obsKey=None, obsVal=None, normReads=False, combat=False):
    for k in keys:
        adataDict[k].obs['batch'] = k
    adata = ad.concat([adataDict[k] for k in keys])
    adata.obs_names_make_unique()
    if normReads:
        adata = equalizeReadMeans(adata)
    if combat:
        sc.pp.combat(adata)
    if obsVal:
        adata = adata[adata.obs[obsKey] == obsVal]
    return adata


def tablize(adata, index=None, treatment=None, cellType=None, counts=False):
    arr = np.empty((0, len(adata.uns['rank_genes_groups']['names'][adata.obs['batch'][-1]])))
    cols = []
    if treatment is not None:
        if type(treatment) == list:
            cols.append('treatment')
            arr = np.vstack((arr, [f"{treatment[0]} x {treatment[1]}"] * arr.shape[1]))
            cols.append('baseTreatment')
            arr = np.vstack((arr, [treatment[0]] * arr.shape[1]))
            cols.append('contrastTreatment')
            arr = np.vstack((arr, [treatment[1]] * arr.shape[1]))
        else:
            if type(treatment) == str:
                cols.append('treatment')
                treatment = [treatment] * arr.shape[1]
                arr = np.vstack((arr, treatment))
    if cellType is not None:
        cols.append('cellType')
        arr = np.vstack((arr, [cellType] * arr.shape[1]))
    for k in adata.uns['rank_genes_groups'].keys():
        if k != 'params':
            arr = np.vstack((arr, adata.uns['rank_genes_groups'][k][adata.obs['batch'][-1]]))
            if k == 'names':
                cols.append('genes')
            else:
                cols.append(k)
        if counts:
            for t in treatment:
                cols.append(t)
                for gene in adata.uns['rank_genes_groups']['names'][adata.obs['batch'][-1]]:
                    meanCounts.append(np.mean(adata[(adata.obs['batch'] == t, gene)].X))

            else:
                arr = np.vstack((arr, meanCounts))

        if index:
            index = adata.uns['rank_genes_groups']['names'][adata.obs['batch'][-1]]
    return pd.DataFrame((arr.T), index=index, columns=cols)


def doDGE(adata, cellTypes):
    pass


def enrich(glist, species='Mouse', analyses=None, treatment='', title_suffix='', return_df=False):
    if analyses == None:
        analyses = [
         'KEGG_2019_Mouse', 'GO_Biological_Process_2021', 'WikiPathways_2019_Mouse']
    if treatment != '':
        treatment = treatment + '\n'
    dfs = []
    for analysis in analyses:
        enr_res = gseapy.enrichr(gene_list=glist, organism=species,
          gene_sets=analysis,
          cutoff=0.5)
        if len(enr_res.res2d[enr_res.res2d['Adjusted P-value'] <= 0.05]) > 0:
            print(len(enr_res.res2d[enr_res.res2d['Adjusted P-value'] <= 0.05]))
            gseapy.barplot((enr_res.res2d), title=(treatment + analysis + title_suffix))
        else:
            print('No Significant Stuffs')
        dfs.append(enr_res.results)
    if return_df:
        return dfs
