import numpy as np
import pandas as pd
import pegasus as pg
import json
import os
import pegasusio
from pegasusio import UnimodalData, MultimodalData
import sys
from utils.utils import *

data = pg.read_input('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/RNA_FINAL.zarr')
print(data)

print('PTSD and CON subclass')
print(data.obs[data.obs.Condition.isin(['PTSD','CON'])]['subclass'].value_counts())

print('MDD and CON subclass')
print(data.obs[data.obs.Condition.isin(['MDD','CON'])]['subclass'].value_counts())

print('PTSD and CON class')
print(data.obs[data.obs.Condition.isin(['PTSD','CON'])]['class'].value_counts())

print('MDD and CON class')
print(data.obs[data.obs.Condition.isin(['MDD','CON'])]['class'].value_counts())

def save_h5ad_and_metadata(data,condition,data_folder,anno):
    data.obs.to_csv(f'{data_folder}{condition}_CON_{anno}_meta.txt', sep='\t')
    data.obs = data.obs[['Channel']]
    data.obsm['X_pca'] = data.obsm['X_pca_harmony']
    data.obsm = pegasusio.datadict.DataDict(dict((k,data.obsm[k]) for k in ('X_pca','X_umap')))
    data.varm = pegasusio.datadict.DataDict()
    data.uns = pegasusio.datadict.DataDict()
    pg.write_output(data, data_folder+f'{condition}_CON_{anno}.h5ad')

for celltype in subclasses_old+['EXC','INH']:
    if celltype not in ['EXC','INH']:
        for condition in ['PTSD','MDD']:
            subset = data[(data.obs['subclass']==celltype) & (data.obs['Condition'].isin([condition,'CON']))].copy()
            print(celltype,condition,subset.X.shape)
            save_h5ad_and_metadata(subset,condition,f'/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{condition}_vs_CON/',celltype)
    else:
        for condition in ['PTSD','MDD']:
            subset = data[(data.obs['class']==celltype) & (data.obs['Condition'].isin([condition,'CON']))].copy()
            print(celltype,condition,subset.X.shape)
            save_h5ad_and_metadata(subset,condition,f'/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{condition}_vs_CON/',celltype)
    
