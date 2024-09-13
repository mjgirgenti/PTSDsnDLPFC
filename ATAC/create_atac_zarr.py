import numpy as np
import pandas as pd
import pegasus as pg
import json
import os
import pegasusio
from pegasusio import UnimodalData, MultimodalData
import matplotlib.pyplot as plt 
import scipy.sparse
from scripts.plot_utils import *
from scripts.utils import *
import matplotlib.patches as mpatches
import itertools
from matplotlib import cm
from mycolorpy import colorlist as mcp
import random
from scipy.spatial import ConvexHull
import scipy.io as sio
import muon as mu
from muon import MuData
import scanpy as sc, anndata as ad
import h5py

genes = pd.read_csv('/home/ah2428/girgenti_scratch60/atac/all_genes.txt',sep='\t')
genes = genes.x.values

atac_meta = pd.read_csv('/home/ah2428/project/PTSD_ATAC/atac_meta.csv',sep='\t')

atac_umap = pd.read_csv('/home/ah2428/project/PTSD_ATAC/umap.txt',sep='\t')

f1 = h5py.File('/home/ah2428/girgenti_scratch60/atac/proj1/mat1/assays.h5', 'r')
f2 = h5py.File('/home/ah2428/girgenti_scratch60/atac/proj2/mat2/assays.h5', 'r')

mat1 = f1['assay001'][:]
mat2 = f2['assay001'][:]

adata1 = sc.AnnData(mat1)
adata2 = sc.AnnData(mat2)

adata_concat = ad.concat([adata1,adata2])

adata_concat.obs_names = atac_meta.index.values
adata_concat.var_names = genes

atac = MultimodalData(adata_concat)

atac.var_names = genes
atac.obs_names = atac_meta.index.values

atac.obsm['X_umap'] = atac_umap.to_numpy()

atac.obs['anno'] = list(atac_meta.Int_Cluster.values)

pg.scatter(atac,attrs='anno')
plt.savefig('/home/ah2428/atac2.png')

pg.write_output(atac,'/home/ah2428/girgenti_scratch60/atac/atac2.zarr')