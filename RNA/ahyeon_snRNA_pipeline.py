import pegasus as pg
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import subprocess
import time
import json
import csv
import multiprocessing as mp
from pegasusio import UnimodalData, MultimodalData

start_time = time.time()

input = sys.argv[1]
out_folder = sys.argv[2]
save_file = sys.argv[3]
res = float(sys.argv[4])

if not os.path.exists(out_folder):
    os.makedirs(out_folder)

print('----------LOAD INPUT FILE: DOUBLET REMOVED LIBRARIES AGGREGATED----------') 
data = pg.read_input(input)
print(data.X.shape)

print('----------IDENTIFY ROBUST GENES----------')
pg.identify_robust_genes(data)

barcode_list = []
for ii,gene in enumerate(data.var_names):
    if data.var.robust[ii]==True:
        barcode_list.append(True)
    else:
        barcode_list.append(False)
data = data[:, barcode_list].copy()
data = MultimodalData(data)
print(data.X.shape)

pg.log_norm(data)

print('----------BATCH CORRECTION & LEIDEN & UMAP----------')
pg.highly_variable_features(data, batch='Channel', n_top=2000)
pg.pca(data)
pca_key = pg.run_harmony(data, n_jobs=15)
pg.neighbors(data, rep=pca_key, n_jobs=15)
pg.leiden(data, resolution=res, rep=pca_key)
pg.umap(data, rep=pca_key, n_jobs=15)
print(data)

print('----------ANNOTATE CLUSTERS----------')
pg.de_analysis(data, cluster='leiden_labels')

hodge_markers = '/home/ah2428/girgenti_project/pegasus/scripts/markers/Hodge_three-layer_hierarchical_cortical_markers.json'
lake_markers = '/home/ah2428/girgenti_project/pegasus/scripts/markers/Lake_cortical_markers.json'
nenad_markers = '/home/ah2428/girgenti_project/pegasus/scripts/markers/Ma_Sestan_hierarchical_markers.json'

hodge_dict = pg.infer_cell_types(data,markers=hodge_markers,output_file=out_folder+'Hodge_infer_cell_types.txt')
cluster_names = pg.infer_cluster_names(hodge_dict)
pg.annotate(data, name='hodge_anno', based_on='leiden_labels', anno_dict=cluster_names)
pg.scatter(data, attrs='hodge_anno', basis='umap', legend_loc='on data', dpi=300)
plt.savefig(out_folder+'hodge_anno.png',bbox_inches='tight')

lake_dict = pg.infer_cell_types(data,markers=lake_markers,output_file=out_folder+'Lake_infer_cell_types.txt')
cluster_names = pg.infer_cluster_names(lake_dict)
pg.annotate(data, name='lake_anno', based_on='leiden_labels', anno_dict=cluster_names)
pg.scatter(data, attrs='lake_anno', basis='umap', legend_loc='on data', dpi=300)
plt.savefig(out_folder+'lake_anno.png',bbox_inches='tight')

nenad_dict = pg.infer_cell_types(data,markers=nenad_markers,output_file=out_folder+'Nenad_infer_cell_types.txt')
cluster_names = pg.infer_cluster_names(nenad_dict)
pg.annotate(data, name='nenad_anno', based_on='leiden_labels', anno_dict=cluster_names)
pg.scatter(data, attrs='nenad_anno', basis='umap', legend_loc='on data', dpi=300)
plt.savefig(out_folder+'nenad_anno.png',bbox_inches='tight')

print('----------SAVE OUTPUT FILE----------')
pg.write_output(data, save_file)

print('----------GENE EXPRESSION PLOTS----------')
pg.scatter(data, attrs=['lake_anno','n_genes'], wspace = 0.1, panel_size=(3,3),legend_loc=['on data','right margin'],basis='umap',dpi=300)
plt.savefig(out_folder+'ngenes.png',bbox_inches='tight')
pg.dotplot(data, genes=['SLC17A7','GAD1','FLT1','AQP4','MOBP','PDGFRA','APBB1IP'], groupby='hodge_anno')
plt.savefig(out_folder+'simple_dotplot.png', bbox_inches='tight')

hours = ((time.time()-start_time)/60)/60
print('Total time = {} hours'.format(hours))

