import pegasus as pg
import numpy as np

data = pg.read_input('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/RNA_FINAL.zarr')
data

mat = data.X.toarray()

idx = [list(data.obs_names).index(x) for x in data.obs[data.obs.fine_anno=='AST AQP4 SLC1A2'].index]

m = data.X[idx,:].toarray()

np.savez('/home/ah2428/AST_AQP4_SLC1A2.npz',m)

idx = [list(data.obs_names).index(x) for x in data.obs[data.obs.fine_anno=='AST GFAP FABP7'].index]
m = data.X[idx,:].toarray()

np.savez('/home/ah2428/AST_GFAP_FABP7.npz',m)
