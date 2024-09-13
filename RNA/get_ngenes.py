import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from utils.lists import *
from utils.plot_utils import *
import itertools
import warnings
warnings.filterwarnings("ignore")

data = pg.read_input('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/RNA_FINAL.zarr')
data

data.select_matrix('raw.X')

mat = data.X.toarray()

ngenes = mat.sum(axis=1)

data.obs['ngenes'] = ngenes
data.obs['log10ngenes'] = np.log10(ngenes)

pg.write_output(data,'/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/RNA_FINAL.zarr')

