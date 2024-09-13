import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pickle as pkl
import scanpy as sc
from rnanorm import CPM
import conorm
import warnings
import scipy.stats as stats

data = sc.read_h5ad(f'/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/MDD_vs_CON/MDD_CON_EXC.h5ad')
#data = sc.read_h5ad(f'/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/PTSD_with_MDD_vs_CON/PTSD_CON_EXC.h5ad')

protein = pd.read_csv('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt',sep='\t',header=None)[1].values
protein_genes = list(set(data.var_names) & set(protein))
pd.DataFrame(protein_genes).to_csv('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/MDD_vs_CON/EXC_genes.csv',header=None,index=None)

nf = pd.read_csv('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/MDD_vs_CON/EXC_norm_factors.csv',sep='\t')
nf = nf[['norm.factors']]

mat = pd.DataFrame(data.raw.X.toarray())
mat.index = data.obs_names
mat.columns = data.var_names
mat = mat[protein_genes]

cpms = conorm.cpm(mat.T,norm_factors=nf)
print(cpms.shape)

cpms.to_csv('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/MDD_vs_CON/EXC_cpms.csv.gz',index=None,header=None,compression='gzip')
