import numpy as np
import pandas as pd
import pegasus as pg
import json
import os
import pegasusio
from pegasusio import UnimodalData, MultimodalData
import matplotlib.pyplot as plt 
import scipy.sparse
import itertools
import random
import scipy.io as sio
import scanpy as sc, anndata as ad
import h5py
import re
import subprocess
import multiprocessing as mp
import doubletdetection as dd
import scrublet as scr
import sys

sample = sys.argv[1]
file = f"/home/ah2428/girgenti_scratch60/scRNA/cellbender/{sample}/{sample}_cellbender_filtered.h5"
data = pg.read_input(file)

mito_df = pd.read_csv('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/ahyeon/PEC_scRNA_pipeline_code/mito_genes.csv')
mito_df = mito_df.loc[0:1135,'HumanGeneID':'Symbol']
mito_list = []
for i in range(mito_df.shape[0]):
    mito_list.append(mito_df.loc[i,'Symbol'])
    mito_list = set(mito_list)
    mito_list = list(mito_list)
    
for ii,gene in enumerate(data.var.index):
    if gene in mito_list:
        data.var = data.var.rename(index={gene:'Mito-{}'.format(gene)})
        
pg.qc_metrics(data,mito_prefix='Mito-',percent_mito=10,min_genes=200,min_umis=500)
pg.filter_data(data)
print('Mito filtered:', data.X.shape)

mito_genes = []
for i in data.var.index:
    if 'Mito' in i:
        mito_genes.append(i)
        
r = re.compile("^RP[SL]")
ribo_genes = list(filter(r.match, data.var.index))

for ii,gene in enumerate(data.var.index):
    if gene in ribo_genes:
        data.var = data.var.rename(index={gene:'Ribo-{}'.format(ii)})

ribo_genes = []
for i in data.var.index:
    if 'Ribo' in i:
        ribo_genes.append(i)
        
pg.qc_metrics(data,mito_prefix='Ribo-',percent_mito=2,min_genes=200,min_umis=500)
pg.filter_data(data)
print('Ribo filtered:', data.X.shape)

protein_coding = pd.read_csv('/home/ah2428/girgenti_project/pegasus/data/protein_coding.genes.with.chr.txt',sep='\t',header=None)
sex_genes = protein_coding[protein_coding[0]=='chrY'].iloc[:,1].values
remove_genes = list(itertools.chain(mito_genes, sex_genes, ribo_genes))

barcode_list = []
for i in data.var_names:
    if i in remove_genes:
        barcode_list.append(False)
    else:
        barcode_list.append(True)
data_final = data[:, barcode_list].copy()
data = MultimodalData(data_final)
print('Genes filtered:', data.X.shape)

pg.identify_robust_genes(data)
pg.log_norm(data)

barcode_list = []
for ii,gene in enumerate(data.var_names):
    if data.var.robust[ii]==True:
        barcode_list.append(True)
    else:
        barcode_list.append(False)
data = data[:, barcode_list].copy()
data = MultimodalData(data)
print('Robust genes:', data.X.shape)

print('Beginning scrublet')
counts_matrix = data.X
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
doublet = scrub.predicted_doublets_
data.obs["doublet"] = doublet
data_subset = data[data.obs["doublet"] == False,:].copy()
data = data_subset
data = MultimodalData(data)
print(data.X.shape)

print('Beginning DD')
clf = dd.BoostClassifier(n_iters=25, use_phenograph=False, standard_scaling=True)
doublets = clf.fit(data.X).predict(p_thresh=1e-16, voter_thresh=0.3)
doublet_score = clf.doublet_score()
data.obs["doublet"] = doublets
data.obs["doublet_score"] = doublet_score
data_subset = data[data.obs["doublet"] == 0,:].copy()
data = data_subset
data = MultimodalData(data)
print(data.X.shape)

pg.highly_variable_features(data, consider_batch=True)
pg.pca(data)
pca_key = pg.run_harmony(data)
pg.neighbors(data, rep=pca_key, full_speed=True)
pg.leiden(data, rep=pca_key)
pg.umap(data, rep=pca_key, full_speed=True)

pg.de_analysis(data, cluster='leiden_labels')
cell_type_dict = pg.infer_cell_types(data,markers='human_brain')
cluster_names = pg.infer_cluster_names(cell_type_dict)
pg.annotate(data, name='anno', based_on='leiden_labels', anno_dict=cluster_names)
pg.scatter(data, attrs='anno', basis='umap', legend_loc = 'on data', dpi=300)
plt.title(sample)
plt.savefig(f"/home/ah2428/girgenti_scratch60/scRNA/onesample_umaps/{sample}_umap.png",bbox_inches='tight')
