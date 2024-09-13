import pegasus as pg
import numpy as np
import pandas as pd
import os
import sys
import subprocess
import multiprocessing as mp
import doubletdetection as dd
import scrublet as scr
from pegasusio import UnimodalData, MultimodalData
import itertools

###Doublet detection -- Scrublet and Doublet Detection

path = '/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/26_samples_multiome/'

sample = sys.argv[1]

data = pg.read_input(path+sample+'/outs/filtered_feature_bc_matrix.h5')

mito_df = pd.read_csv('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/ah2428/PEC_scRNA_pipeline_code/mito_genes.csv')
mito_df = mito_df.loc[0:1135,'HumanGeneID':'Symbol']
mito_list = list(set(mito_df.Symbol))

for ii,gene in enumerate(data.var.index):
    if gene in mito_list:
        data.var = data.var.rename(index={gene:'Mito-{}'.format(gene)})

pg.qc_metrics(data,mito_prefix='Mito-',percent_mito=10,min_genes=200,min_umis=500)
pg.filter_data(data)
print(data.X.shape)

protein_coding = pd.read_csv('/home/ah2428/girgenti_project/ptsd/data/protein_coding.genes.with.chr.txt',sep='\t',header=None)
sex_genes = protein_coding[protein_coding[0]=='chrY'].iloc[:,1].values
mito_genes = [x for x in data.var.index if 'Mito' in x]

remove_genes = list(itertools.chain(mito_genes, sex_genes))

barcode_list = []
for i in data.var_names:
    if i in remove_genes:
        barcode_list.append(False)
    else:
        barcode_list.append(True)
data = data[:, barcode_list].copy()
data = MultimodalData(data)
print(data.X.shape)


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
clf = dd.BoostClassifier(n_iters=25, clustering_algorithm="phenograph", standard_scaling=True)
doublets = clf.fit(data.X).predict(p_thresh=1e-16, voter_thresh=0.3)
doublet_score = clf.doublet_score()
data.obs["doublet"] = doublets
data.obs["doublet_score"] = doublet_score
data_subset = data[data.obs["doublet"] == 0,:].copy()
data = data_subset
data = MultimodalData(data)
print(data.X.shape)

pd.DataFrame(data.obs_names).to_csv(f'/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/rna_doublets_filtered/{sample}_rna_filtered_barcodes.csv',index=None,header=None)
