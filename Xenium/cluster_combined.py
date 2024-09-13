import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pegasus as pg
import squidpy as sq
import seaborn as sns
import sys
import os

file = sys.argv[1]
savepath = sys.argv[2]
pca = sys.argv[3]
neigh = sys.argv[4]
res = sys.argv[5]


if not os.path.exists(savepath): 
    os.makedirs(savepath)

data = sc.read_h5ad(f'{file}')
print(data)

data.layers["counts"] = data.X.copy()
sc.pp.normalize_total(data, inplace=True)
sc.pp.log1p(data)
sc.pp.pca(data,n_comps=int(pca))
sc.external.pp.harmony_integrate(data, 'Channel')
data.obsm['X_pca'] = data.obsm['X_pca_harmony']
sc.pp.neighbors(data,n_neighbors=int(neigh))
sc.tl.umap(data,random_state=0)
sc.tl.leiden(data,resolution=float(res))
sc.pl.umap(data,color='leiden',legend_loc='on data')
plt.savefig(savepath+f'combined_{pca}_{neigh}_{res}_leiden.png',dpi=300)

#sc.pl.dotplot(data,['ANO3','SLC17A7','GAD1','GAD2','PDGFRA','OLIG1','PLP1','MAG','AQP4','GJA1','FLT1','PECAM1','SLC11A1','PTPRC'],groupby='leiden',dendrogram=True)
#plt.savefig(folder+f'combined2_{pca}_{neigh}_{res}_dotplot.png',dpi=300)

data.write(savepath+f'combined_{pca}_{neigh}_{res}_clustered.h5ad')
