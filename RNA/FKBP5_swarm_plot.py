import numpy as np
import pandas as pd
import pegasus as pg
import json
import os
import pegasusio
from pegasusio import UnimodalData, MultimodalData
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import seaborn as sns
import scanpy as sc


# RNA

data = pg.read_input('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/RNA_FINAL.zarr')
data.select_matrix('X')

gene = "FKBP5"
idx = list(data.var_names).index('FKBP5') 

expression_data = pd.DataFrame({
    "Sample": data.obs['Channel'],  # Sample info from metadata
    "Condition": data.obs['Condition'],
    "Celltype": data.obs['class'],
    "Expression": data[:, idx].X.toarray().flatten()  # Extract expression values
})

plt.rcParams['font.size']='24'
fig, ax = plt.subplots(1,8,figsize=(35,5),sharey=True)
ax = ax.flatten()
for i,celltype,celltype2 in zip(np.arange(8),['ALL','EXC','INH','OLI','OPC','END','AST','MIC'],['ALL','EXN','IN','OLG','OPC','END','AST','MG']):
    if celltype=='ALL':
        exp = expression_data.groupby(['Sample','Condition']).mean('Expression').reset_index()
    else:
        exp = expression_data[expression_data.Celltype==celltype]
        exp = exp.groupby(['Sample','Condition']).mean('Expression').reset_index()
    exp = exp.dropna()
    exp['Condition'] = exp['Condition'].cat.reorder_categories(['CON','MDD','PTSD'])
    sns.violinplot(data=exp, x="Condition", y="Expression", inner='points',palette=['tab:gray','tab:orange','tab:blue'],zorder=2, ax=ax[i])
    ax[i].set_title(f"{celltype2}")
    ax[i].set_xlabel(None)
    ax[i].set_ylabel(None)
    ax[i].grid(alpha=0.5)

ax[0].set_ylabel("FKBP5 Average\nGene Expression\n(RNA)", fontsize=24, labelpad=5)

plt.tight_layout()
plt.savefig('/home/ah2428/review2_figs/FKBP5_rna.pdf',bbox_inches='tight',dpi=300)
plt.show()


# Xenium

xenium = sc.read_h5ad('/home/ah2428/girgenti/ah2428/xenium/data_processed_resegment/snXenium.h5ad')

gene = "FKBP5"
idx = list(xenium.var_names).index('FKBP5') 

expression_data = pd.DataFrame({
    "Sample": xenium.obs['Channel'],  # Sample info from metadata
    "Condition": xenium.obs['Condition'],
    "Celltype": xenium.obs['celltype'],
    "Expression": xenium[:, idx].X.toarray().flatten()  # Extract expression values
})

plt.rcParams['font.size']='24'
fig, ax = plt.subplots(1,8,figsize=(35,5),sharey=True)
ax = ax.flatten()
for i,celltype,celltype2 in zip(np.arange(8),['ALL','EXC','INH','OLI','OPC','END','AST','MIC'],['ALL','EXN','IN','OLG','OPC','END','AST','MG']):
    if celltype=='ALL':
        exp = expression_data.groupby(['Sample','Condition']).mean('Expression').reset_index()
    else:
        exp = expression_data[expression_data.Celltype==celltype]
        exp = exp.groupby(['Sample','Condition']).mean('Expression').reset_index()
    exp = exp.dropna()
    exp['Condition'] = exp['Condition'].cat.reorder_categories(['CON','MDD','PTSD'])
    sns.violinplot(data=exp, x="Condition", y="Expression", inner='points',palette=['tab:gray','tab:orange','tab:blue'],zorder=2, ax=ax[i])
    ax[i].set_title(f"{celltype2}")
    ax[i].set_xlabel(None)
    ax[i].set_ylabel(None)
    ax[i].grid(alpha=0.5)

ax[0].set_ylabel("FKBP5 Average\nGene Expression\n(Xenium)", fontsize=24)

plt.tight_layout()
plt.savefig('/home/ah2428/review2_figs/FKBP5_xenium.pdf',bbox_inches='tight',dpi=300)
plt.show()