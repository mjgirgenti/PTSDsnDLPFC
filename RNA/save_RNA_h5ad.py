import pegasus as pg
import pandas as pd
import numpy as np
import pegasusio
from pegasusio import UnimodalData, MultimodalData

data = pg.read_input('/home/ah2428/ShareZhangLab/PTSD/RNA/data/RNA_FINAL.zarr')
channels = data.obs.Channel.unique()
atac_genes = pd.read_csv('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/genes.csv',sep='\t',header=None)[0].values
data2 = data[:,atac_genes].copy()

data3 = data2[data2.obs.Channel.isin(channels[50:])].copy()

def save_h5ad_and_metadata(data):
    data.obs = data.obs[['Channel']]
    data.obsm['X_pca'] = data.obsm['X_pca_harmony']
    data.obsm = pegasusio.datadict.DataDict(dict((k,data.obsm[k]) for k in ('X_pca','X_umap')))
    data.varm = pegasusio.datadict.DataDict()
    data.uns = pegasusio.datadict.DataDict()
    data.obsp = pegasusio.datadict.DataDict()
    return data

data4 = save_h5ad_and_metadata(data3)
print(data4)
pg.write_output(data4,'/home/ah2428/ShareZhangLab/PTSD/RNA/data/RNA_FINAL2.h5ad')
