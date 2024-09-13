import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pegasus as pg
import squidpy as sq
import seaborn as sns
import sys
import os
import pandas as pd

brnum = sys.argv[1]

folder = f'/gpfs/gibbs/pi/girgenti/ah2428/xenium/single_sample_QC/'

if not os.path.exists(folder): 
    os.makedirs(folder)

meta = pd.read_csv('/gpfs/gibbs/pi/girgenti/ah2428/xenium/xenium_meta.csv',sep='\t')
condition = meta[meta.BrNum==int(brnum)].Condition.values[0]

data = pg.read_input(f'/gpfs/gibbs/pi/girgenti/ah2428/xenium/data/{brnum}/cell_feature_matrix_pg.h5')
data = data.to_anndata()
sc.pp.calculate_qc_metrics(data, percent_top=(10, 20, 50, 150), inplace=True)
sc.pp.filter_cells(data, min_counts=10)

fig, axs = plt.subplots(1, 2, figsize=(8, 4))
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    data.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    data.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)
fig.suptitle(f'BrNum {brnum} {condition}')
fig.tight_layout()
plt.savefig(folder+f'{brnum}_{condition}_qc.png',dpi=300)

data.layers["counts"] = data.X.copy()
sc.pp.normalize_total(data, inplace=True)
sc.pp.log1p(data)
sc.pp.pca(data,n_comps=50)
sc.pp.neighbors(data,n_neighbors=100)
sc.tl.umap(data,random_state=0)
sc.tl.leiden(data,resolution=1.0)
sc.pl.umap(data,color='leiden',legend_loc='on data',title=f'BrNum {brnum} {condition}')
plt.savefig(folder+f'{brnum}_{condition}_leiden.png',dpi=300)

sc.pl.dotplot(data,['ANO3','SLC17A7','GAD1','GAD2','PDGFRA','OLIG1','PLP1','MAG','AQP4','GJA1','FLT1','PECAM1','SLC11A1','PTPRC'],groupby='leiden',dendrogram=True)
plt.savefig(folder+f'{brnum}_{condition}_dotplot.png',dpi=300)

