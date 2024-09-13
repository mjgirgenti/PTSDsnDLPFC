import numpy as np
import pandas as pd
import pegasus as pg
import json
import os
import pegasusio
from pegasusio import UnimodalData, MultimodalData
import matplotlib.pyplot as plt 
import scipy.sparse
import matplotlib.patches as mpatches
import itertools
from matplotlib import cm
import matplotlib
import random
from scipy.spatial import ConvexHull
import scipy.io as sio
import scanpy as sc, anndata as ad
import h5py
import re
from scripts.utils import *
from scipy.stats import wilcoxon, mannwhitneyu, ttest_ind
import doubletdetection as dd
import scrublet as scr

path = '/ysm-gpfs/pi/gerstein/jz435/Share-ZhangLab/scRNA/cellbender_merged/'
out_path = '/home/ah2428/girgenti_scratch60/scRNA/pegasus/data/ahyeon_mito_ribo_filtered/'
sample = sys.argv[1]
ribo = sys.argv[2]

# input cellbender file
file = path+'{}/{}_cellbender_filtered.h5'.format(sample,sample)
data = pg.read_input(file)
print(data.X.shape)

# filter cells greater than 10% mito
mito_df = pd.read_csv('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/ahyeon/PEC_scRNA_pipeline_code/mito_genes.csv')
mito_df = mito_df.loc[0:1135,'HumanGeneID':'Symbol']
mito_list = list(set(mito_df.Symbol))

for ii,gene in enumerate(data.var.index):
    if gene in mito_list:
        data.var = data.var.rename(index={gene:'Mito-{}'.format(gene)})
        
pg.qc_metrics(data,mito_prefix='Mito-',percent_mito=10,min_genes=200,min_umis=500)
pg.filter_data(data)
print(data.X.shape)

# filter cells greater than 2% ribo
r = re.compile("^RP[SL]")
ribo_genes = list(filter(r.match, data.var.index))
for ii,gene in enumerate(data.var.index):
    if gene in ribo_genes:
        data.var = data.var.rename(index={gene:'Ribo-{}'.format(ii)})
if ribo=='T':
    pg.qc_metrics(data,mito_prefix='Ribo-',percent_mito=2,min_genes=200,min_umis=500)
    pg.filter_data(data)
    print(data.X.shape)

# filter mito,ribo,sex genes
protein_coding = pd.read_csv('/home/ah2428/girgenti_project/pegasus/data/protein_coding.genes.with.chr.txt',sep='\t',header=None)
sex_genes = protein_coding[protein_coding[0]=='chrY'].iloc[:,1].values
mito_genes = [x for x in data.var.index if 'Mito' in x]
ribo_genes = [x for x in data.var.index if 'Ribo' in x]
remove_genes = list(itertools.chain(mito_genes, sex_genes, ribo_genes))

barcode_list = []
for i in data.var_names:
    if i in remove_genes:
        barcode_list.append(False)
    else:
        barcode_list.append(True)
data = data[:, barcode_list].copy()
data = MultimodalData(data)
print(data.X.shape)

# Doublet detection -- Scrublet
print("Beginning doublet detection - scrublet")
#data.select_matrix('raw.X')
counts_matrix = data.X
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
doublet = scrub.predicted_doublets_
if doublet is None:
    doublet = scrub.call_doublets(threshold=0.3)
data.obs["doublet"] = doublet
data_subset = data[data.obs["doublet"] == False,:].copy()
data = data_subset
data = MultimodalData(data)
print(data.X.shape)

# Doublet detection -- Doublet Detection
print("Beginning doublet detection - DD")
clf = dd.BoostClassifier(n_iters=25, clustering_algorithm="phenograph", standard_scaling=True)
#data.select_matrix('raw.X')
doublets = clf.fit(data.X).predict(p_thresh=1e-16, voter_thresh=0.3)
doublet_score = clf.doublet_score()
data.obs["doublet"] = doublets
data.obs["doublet_score"] = doublet_score
data_subset = data[data.obs["doublet"] == 0,:].copy()
data = data_subset
data = MultimodalData(data)
print(data.X.shape)

pg.write_output(data,out_path+'{}_filtered.zarr'.format(sample))