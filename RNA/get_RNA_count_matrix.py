import pegasus as pg
import numpy as np

data = pg.read_input('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/RNA_FINAL.zarr')
data

data.select_matrix('raw.X')

mat = data.X.toarray()

np.savez('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/RNA_count_matrix.npz',mat=mat)

