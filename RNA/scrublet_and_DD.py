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

###Doublet detection -- Scrublet and Doublet Detection

library_name = sys.argv[1]
print(library_name)

data = pg.read_input('/home/ah2428/girgenti_scratch60/scRNA/pegasus/data/111_samples/lognormed/{}_lognormed.zarr'.format(library_name))

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

pg.write_output(data,'/home/ah2428/girgenti_scratch60/scRNA/pegasus/data/111_samples/doublet_removed/{}_doublet_removed.zarr'.format(library_name))
