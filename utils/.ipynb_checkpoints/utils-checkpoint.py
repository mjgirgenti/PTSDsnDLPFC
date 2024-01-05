import pegasus as pg
import pegasusio
from pegasusio import UnimodalData, MultimodalData
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import subprocess
import random
import json
import csv


protein_coding = pd.read_csv('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/protein_coding.genes.with.chr.txt',sep='\t',header=None)[1].values

lake_genes = ['MRC1','TMEM119','CX3CR1','APBB1IP','CLDN5','FLT1','DUSP1','COBLL1','GLUL','SOX9','AQP4','GJA1','NDRG2','GFAP','ALDH1A1','ALDH1L1','VIM','SLC4A4','PDGFRA','PCDH15','OLIG2','OLIG1','PLP1','MAG','MOG','MOBP','MBP','SATB2','SLC17A7','GRM4','GAD2','GAD1','SLC32A1','SST','PVALB','STMN2','RBFOX3']

celltypes = ['EXC','INH','OLI','OPC','END','AST','MIC']

subtypes = ['CUX2','RORB','FEZF2','OPRK1','LAMP5','KCNG1','VIP','SST','PVALB','OLI','OPC','END','AST','MIC']

exc_subtypes = ['CUX2','RORB','FEZF2','OPRK1']

inh_subtypes = ['LAMP5','KCNG1','VIP','SST','PVALB']

celltype_colors = {
    'EXC':'#b22222',
    'INH':'#2E8B57',
    'OLI':'#5254a3',
    'OPC':'#aec7e8',
    'END':'#bc80bd',
    'AST':'#ffed6f',
    'MIC':'#7f7f7f'}

subtype_colors = {
    'CUX2': '#ad494a',
    'RORB': '#e7969c',
    'FEZF2': '#843c39',
    'OPRK1': '#d6616b',
    'LAMP5': '#637939',
    'KCNG1': '#a6d854',
    'VIP': '#cedb9c',
    'SST': '#8ca252',
    'PVALB': '#b5cf6b',
    'OLI': '#5254a3',
    'OPC': '#aec7e8',
    'END': '#bc80bd',
    'AST': '#ffed6f',
    'MIC': '#7f7f7f'}



# To convert h5ad to h5seurat for DEG analysis
def save_h5ad_and_metadata(data,data_folder,name):
    data.obs.to_csv(data_folder+f'{name}.csv')
    data.obs = data.obs[['Channel']]
    data.obsm['X_pca'] = data.obsm['X_pca_harmony']
    data.obsm = pegasusio.datadict.DataDict(dict((k,data.obsm[k]) for k in ('X_pca','X_umap')))
    data.varm = pegasusio.datadict.DataDict()
    data.uns = pegasusio.datadict.DataDict()
    pg.write_output(data, data_folder+f'{name}.h5ad')
    
def get_cluster(data,cluster):
    barcode_list = []
    for ii,name in enumerate(data.obs_names):
        if data.obs.anno[ii] in cluster:
            barcode_list.append(True)
        else:
            barcode_list.append(False)
    cluster = data[barcode_list].copy()
    cluster = MultimodalData(cluster)
    return cluster

def get_louvain(data,cluster):
    barcode_list = []
    for ii,name in enumerate(data.obs_names):
        if data.obs.louvain_labels[ii] in cluster:
            barcode_list.append(True)
        else:
            barcode_list.append(False)
    cluster = data[barcode_list].copy()
    cluster = MultimodalData(cluster)
    return cluster

def get_library_data(data,library_name):
    barcode_list = []
    for i in data.obs.Channel:
        if i == library_name:
            barcode_list.append(True)
        else:
            barcode_list.append(False)
    library_data = data[barcode_list].copy()
    library_data = MultimodalData(library_data)
    return library_data

def subcluster(data,n_clusters,markers,folder):
    leiden_label = data.obs.leiden_labels[0]
    name = data.obs.anno[0]
    pg.split_one_cluster(data,'leiden_labels',leiden_label,n_clusters,'split')
    pg.de_analysis(data, cluster='split', t=True)
    cell_type_dict = pg.infer_cell_types(data, markers=markers)
    cluster_names = pg.infer_cluster_names(cell_type_dict)
    pg.annotate(data, name='anno', based_on='split', anno_dict=cluster_names)
    pg.scatter(data, attrs='anno', basis='umap', dpi=300, panel_size=(3,3))
    plt.savefig(folder+'{}_subclusters.png'.format(name),bbox_inches='tight')

def get_remove_barcodes(cluster,subclusters_to_remove,cluster_str,folder):
    remove_list = []
    for ii,name in enumerate(cluster.obs_names):
        if cluster.obs.anno[ii] in subclusters_to_remove:
            remove_list.append(name)
    remove_series = pd.Series(remove_list)
    remove_series.to_csv(folder+'{}_remove.txt'.format(cluster_str),index=False)
    return remove_list

def subset_samples(data,n,n_iter):
    for i in range(n_iter):
        print(i)
        samples = list(data.obs.Channel.unique())
        random.shuffle(samples)
        samples = samples[:n]

        barcode_list = []
        for sample in data.obs.Channel:
            if sample in samples:
                barcode_list.append(True)
            else:
                barcode_list.append(False)
        data_subset = data[barcode_list].copy()
        data_subset = MultimodalData(data_subset)
        print(data_subset)
        print(list(data_subset.obs.Channel.unique()))
        #pg.write_output(data_subset,'/home/ah2428/girgenti_scratch60/scRNA/pegasus/data/subset{}.zarr'.format(i+1))
    
def efficient_subset(data,n):
    data = pg.read_input('/home/ah2428/scratch60/scRNA/pegasus/data/exc.zarr')

    lst = list(np.arange(0,data.X.shape[0]))
    random.shuffle(lst)
    idx = lst[:n]

    barcodes = [False]*data.X.shape[0]
    for i in idx:
        barcodes[i] = True
    data_subset = data[barcodes].copy()
    subset = MultimodalData(data_subset)
    return subset

def replace_subclass_name(df):
    df['subclass'][df['subclass'].isin(exc_subclasses)] = 'EXC'
    df['subclass'][df['subclass'].isin(inh_subclasses)] = 'INH'
    return df

def get_unique_genes_for_class(df):
    unique_genes_df = []
    for g in df.DEG.unique():
        df_new = df[df.DEG==g].reset_index(drop=True)
        if df_new.shape[0] >= 2:
            idx = df_new.fc.abs().idxmax()
            x = df_new.loc[idx]
            x = pd.DataFrame(x).T
        else:
            x = df_new
        unique_genes_df.append(x) 
    df_final = pd.concat(unique_genes_df)
    return df_final