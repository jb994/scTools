
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


def dgeNorm(adata, useRaw=True, log=True, scale=False, reduceForScale=True, pca=False, n_comps=50):
    if useRaw:
        adata = adata.raw.to_adata()
    sc.pp.normalize_total(adata)
    if log:
        sc.pp.log1p(adata)
    if scale:
        if reduceForScale==True:
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            #sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=4000)
            adata=adata[:,adata.var.highly_variable].copy()
        elif type(reduceForScale)==int:
            #sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=reduceForScale)
            adata=adata[:,adata.var.highly_variable].copy()
        sc.pp.scale(adata)
    if pca:
        sc.tl.pca(adata, n_comps=n_comps)
    return adata


def tablize(adata, index=None, treatment=None, cellType=None, counts=False):
    arr = np.empty((0, len(adata.uns['rank_genes_groups']['names'][adata.obs['batch'][-1]])))
    cols = []
    if treatment is not None:
        arr, cols, treatment = formatTreatment(arr, cols, treatment)
    if cellType is not None:
        arr, cols = formatCellType(arr, cols, cellType)

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
        arr = np.vstack((arr, meanCounts))
    if index:
        index = adata.uns['rank_genes_groups']['names'][adata.obs['batch'][-1]]
    return pd.DataFrame((arr.T), index=index, columns=cols)


def removeSingles(adata, key='fineClusters'):
    ### Remove cells which are the only instance of a cluster in a certain treatment
    ### Do this because it raises an error for scanpy
    for batch in np.unique(adata.obs['batch']):
        for cellType in np.unique(adata.obs[key]):
            if adata[(adata.obs['batch']==batch) & (adata.obs[key]==cellType)].shape[0]==1:
                adata = adata[~((adata.obs[key]==cellType) & (adata.obs['batch']==batch))]
    return adata


def formatTreatment(arr, cols, treatment):
    if type(treatment) == list:
        cols.append('treatment')
        arr = np.vstack((arr, [f"{treatment[0]} x {treatment[1]}"] * arr.shape[1]))
        cols.append('baseTreatment')
        arr = np.vstack((arr, [treatment[0]] * arr.shape[1]))
        cols.append('contrastTreatment')
        arr = np.vstack((arr, [treatment[1]] * arr.shape[1]))
    elif type(treatment) == str:
        cols.append('treatment')
        treatment = [treatment] * arr.shape[1]
        arr = np.vstack((arr, treatment))
    return arr, cols, treatment

def formatCellType(arr, cols, cellType):
    cols.append('cellType')
    arr = np.vstack((arr, [cellType] * arr.shape[1]))
    return arr, cols

def checkClusterInAdatas(adataList, cellType, obsKey='fineClusters'):
    for adata in adataList:
        if adata[adata.obs[obsKey]==cellType].shape[0]<2:
            #print(adata[adata.obs[obsKey]==cellType].shape[0])
            return False
    return True
    #return adata[adata.obs[obsKey]==cellType].shape[0]>1

def doDGE(combinedData, treatments, obsKey='fineClusters'):
    ### Pass in large adata and treatments keys to access each batch
    ### Make sure to use 'dgeNorm' before calling this
    dfs=[]
    for i, treat1 in enumerate(treatments):
        print(treat1)
        for j, treat2 in enumerate(treatments):
            if j>i:
                Data1 = combinedData[combinedData.obs['batch']==treat1].copy()
                Data2 = combinedData[combinedData.obs['batch']==treat2].copy()
                cellTypes = sorted(( set(np.unique(Data1.obs[obsKey])) & set(np.unique(Data2.obs[obsKey])) ))
                adata = combinedData[(combinedData.obs['batch']==treat1) | (combinedData.obs['batch']==treat2)].copy()
                for celltype in cellTypes:
                    #if (Data1[Data1.obs[obsKey]==celltype].shape[0]>1) and (Data2[Data2.obs[obsKey]==celltype].shape[0]>1):
                    if checkClusterInAdatas([Data1,Data2], celltype, obsKey=obsKey):
                        dgeDF = getDEGsForCellType(adata, treat1, treat2, obsKey, celltype)
                        dfs.append(dgeDF)
    totalDF = pd.concat(dfs, ignore_index=True)
    return totalDF


def getDEGsForCellType(adata, treat1, treat2, obsKey, celltype):
    print(f"{treat1} {treat2} {celltype}")
    cellTypeAdata = adata[adata.obs[obsKey]==celltype].copy()
    sc.tl.rank_genes_groups(cellTypeAdata, 'batch', groups=[treat2], reference=treat1, use_raw=False, method='wilcoxon', max_iter=2000, pts=False)
    dgeDF = tablize(cellTypeAdata, index=True, treatment=[treat1,treat2], cellType=celltype, counts=False)
    return dgeDF


def logFoldByCell(adata, t1, t2, genes=None, foldDF=None):
    ### Get the logfoldchange of a cell(s) given a treatment
    ### foldDF should be 
    if genes is not None:
        adata = adata[:,genes]
        if foldDF:
            foldDF = foldDF[foldDF.genes]

    if foldDF is None:
        ### Calculate log fold changes on its own
        pass

    else:
        pass


###################################################################################################
###################################################################################################
###################################################################################################

def addReadCounts(dgeDF, combinedData, treatments, obsKey='fineClusters'):
    ### Adds the read counts for each treatment for each cellType into a dge matrix
    ### Also adds the number of cells for each cellType
    for i,treat in enumerate(treatments):
        print(f"\nOn treatment {treat}")
    
        adata = combinedData[combinedData.obs['batch']==treat]
        col=[]
        col2=[]
        
        print("Calculating EXP")
        adata.X=np.expm1(adata.X)
        
        print("Making meanMatrixDF")
        meanMatrixDF = makeAvgReadDF(dgeDF, adata, obsKey)
        
        print("Making cellMatrix")
        cellCountDF = getCellCounts(dgeDF, adata, obsKey)
        
        print("Appending reads to DF")
        dgeDF = appendReadsDF(dgeDF, meanMatrixDF, cellCountDF, treat)
        
        adata.X=np.log1p(adata.X)
        
    avgReads = [f'{t}_AvgReadsPerCell' for t in treatments]
    cells = [f'{t}_NumCells' for t in treatments]
    cols=dgeDF.columns[:9].tolist()+ avgReads + cells
    dgeDF = dgeDF[cols] 

    return dgeDF

def makeAvgReadDF(df, adata, obsKey='fineClusters'):
    arr=[]
    for cluster in np.unique(df['cellType']):
        #print(cluster)
        if cluster in np.unique(adata.obs[obsKey]):
            arr.append(np.asarray(adata[adata.obs[obsKey]==cluster].X.mean(axis=0))[0])
            #print(np.asarray(adata[adata.obs[obsKey]==cluster].X.mean(axis=0))[0])
        else:
            arr.append([np.nan]*df.shape[1])
            #print([np.nan]*df.shape[1])

    meanMatrixDF = pd.DataFrame(arr, index=np.unique(df['cellType']), columns=adata.var.index)
    return meanMatrixDF

def getCellCounts(df, adata, obsKey='fineClusters'):
    cellCounts=[]
    for cluster in np.unique(df['cellType']):
        cellCounts.append(adata[adata.obs[obsKey]==cluster].shape[0])
    cellCountDF = pd.DataFrame(cellCounts, index=np.unique(df['cellType']), columns=['cellCounts'])
    return cellCountDF

def appendReadsDF(df, meanMatrixDF, cellCountDF, treat):
    col1 = []
    col2 = []
    for index, row in df.iterrows():
        cluster = row['cellType']
        gene = row['genes']
        col1.append(meanMatrixDF.loc[cluster, gene])
        col2.append(cellCountDF.loc[cluster,'cellCounts'])
    df[f'{treat}_AvgReadsPerCell']=col1
    df[f'{treat}_NumCells']=col2
    return df
###################################################################################################
###################################################################################################
###################################################################################################
###Functions to calculate the number of Differentially Expressed Genes (DEGs) between treatments###
def filterDEGs(DGEDF, t1, t2, pvalThreshold=0.05, lfThreshold = 0.75):
    sigDF = DGEDF[DGEDF.pvals<=pvalThreshold]
    sigDF = sigDF[np.abs(sigDF.logfoldchanges)>=lfThreshold]
    sigDF = sigDF[sigDF['baseTreatment'] == t1]
    sigDF = sigDF[sigDF['contrastTreatment'] == t2]
    for i,t in enumerate([t1,t2]):
        avgMask=(sigDF[t+'_AvgReadsPerCell']>0)
        mask = (sigDF['baseTreatment'] != t)
        sigDF = sigDF[(mask)|(avgMask)]
        mask = (sigDF['contrastTreatment'] != t)
        sigDF = sigDF[(mask)|(avgMask)]
    return sigDF

def numDEGs(DGEDF, t1, t2, cellTypes=None, doFilter=True):
    if cellTypes is None:
        cellTypes = np.unique(DGEDF['cellType'])
    if doFilter:
        DGEDF = filterDEGs(DGEDF, t1, t2)
    DEGNumbers = []
    for cellType in cellTypes:
        numDEGs = DGEDF[DGEDF['cellType']==cellType].shape[0]
        DEGNumbers.append(numDEGs)
    return DEGNumbers



###################################################################################################
###################################################################################################
###################################################################################################
### Perform Pathway Enrichment Analysis on a given adata ###

def enrichAllCellTypes(adata, t1, t2, analysis, comparison, obsKey='fineClusters', savePics=False):
    df = []
    cellTypes = sorted(( set(np.unique(adata[adata.obs['batch']==t2].obs[obsKey])) & set(np.unique(adata[adata.obs['batch']==t1].obs[obsKey]) )))
    for cellType in cellTypes:
        print(cellType)
        tempAdata = adata[adata.obs[obsKey]==cellType]
        sc.tl.rank_genes_groups(tempAdata, 'batch', groups=[t2], reference=t1, use_raw=False, method='wilcoxon', max_iter=2000, pts=False, key_added = "wilcoxon")
        glist = sc.get.rank_genes_groups_df(tempAdata, group=t2, 
                                            key='wilcoxon', log2fc_min=0.25, 
                                            pval_cutoff=0.01)['names'].squeeze()
        if type(glist)!=str:
            glist = glist.str.strip().tolist()
        else:
              glist=[glist]

        if len(glist)>0:
            print(len(glist))
            result = enrich(glist, analyses=analysis, treatment=comparison, title_suffix='_'+cellType, return_df=True)[0]
            if savePics:
                plt.savefig(f"images/{comparison}_{cellType}_{analysis}.png", bbox_inches='tight')
            plt.show()
            result = result[result['Adjusted P-value']<0.1]
            result['Cell Type'] = [cellType]*result.shape[0]
            result = result[['Gene_set','Cell Type', 'Term','Overlap','Adjusted P-value','Genes']]

            df.append(result)
        else:
            print("No sig. genes")

    return df

def enrich(glist, species='Mouse', analyses=None, treatment='', title_suffix='', return_df=False):
    if analyses == None:
        analyses = ['KEGG_2019_Mouse', 'GO_Biological_Process_2021', 'WikiPathways_2019_Mouse']
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



