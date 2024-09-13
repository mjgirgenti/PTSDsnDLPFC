import spatialdata as sd
import squidpy as sq
import spatialdata_io
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import spatialdata_plot
import matplotlib.colors as mcolors
import scanpy as sc

data = sc.read_h5ad('/gpfs/gibbs/pi/girgenti/ah2428/xenium/data_processed/xenium_final.h5ad')

sdata = spatialdata_io.xenium('/gpfs/gibbs/pi/girgenti/ah2428/xenium/data/5173/')
sc.pp.normalize_total(sdata.tables["table"])
sc.pp.log1p(sdata.tables["table"])
sc.pp.highly_variable_genes(sdata.tables["table"])
sdata.tables["table"].var.sort_values("means")
sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")

sq.gr.spatial_neighbors(data)
sq.gr.spatial_autocorr(data, mode="moran", genes=data.var_names)

for gene in data.uns['moranI'][:100].index:
    sdata.pl.render_shapes("cell_boundaries",color=gene).pl.show()
    plt.savefig(f'/home/ah2428/palmer_scratch/xenium_svg/{gene}_svg.png',dpi=300)


