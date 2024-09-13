import spatialdata as sd
import spatialdata_io
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import spatialdata_plot
import matplotlib.colors as mcolors
import scanpy as sc
import squidpy as sq
import sys
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--brnum', type=int)
parser.add_argument('--celltype', type=str)
parser.add_argument('--gene', type=str)
args = parser.parse_args()

brnum = args.brnum
celltype = args.celltype
gene = args.gene

folder = '/home/ah2428/palmer_scratch/xenium_slides/'

sdata = spatialdata_io.xenium(f'/gpfs/gibbs/pi/girgenti/ah2428/xenium/data/{brnum}/')
data = sc.read_h5ad('/gpfs/gibbs/pi/girgenti/ah2428/xenium/xenium_subclass.h5ad')
condition = data.obs[data.obs.Channel==int(brnum)].Condition.unique()[0]

sc.pp.normalize_total(sdata.tables["table"])
sc.pp.log1p(sdata.tables["table"])
sc.pp.highly_variable_genes(sdata.tables["table"])
sdata.tables["table"].var.sort_values("means")
sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")

sdata.tables["table"].obs['cellid'] = [x.split('-')[0] for x in sdata.tables["table"].obs['cell_id'].values]
df = data.obs[data.obs.cellid.isin(sdata.tables["table"].obs.cellid)][['cellid','subtype']]
subtypes = []
for cell in sdata.tables['table'].obs['cellid'].values:
    if cell in df.cellid.values:
        subtype = df[df.cellid==cell].subtype.values[0]
    else:
        subtype = 'NA'
    subtypes.append(subtype)
sdata.tables["table"].obs['subtype'] = subtypes
sdata.tables['table'].obs['subtype'] = sdata.tables['table'].obs['subtype'].astype('category')

# Plot 7 celltypes
#sdata.tables['table'].obs.Celltype = sdata.tables['table'].obs.Celltype.cat.reorder_categories(['EXC','INH','OLI','OPC','END','AST','MIC','NA'])
#sdata.pl.render_shapes("cell_boundaries",color='Celltype',groups=['AST', 'END', 'EXC', 'INH', 'MIC', 'NA', 'OLI', 'OPC'],palette=['#ffed6f','#bc80bd','#b22222','#2E8B57','#7f7f7f','white','#5254a3','#aec7e8']).pl.show()
#plt.title(f'{brnum} {condition} celltypes')
#plt.savefig(f'/home/ah2428/palmer_scratch/xenium_celltypes_fkbp5/{brnum}_{condition}_celltypes.png',dpi=300)


# FKBP5 in END
gene_idx = np.where(sdata.table.var['gene_ids'].index.values==gene)[0]
sdata.tables['table'].obs['gene_exp'] = sdata.tables['table'].X[:,gene_idx].toarray()
sdata.tables['table'].obs.loc[sdata.tables['table'].obs.subtype.isin([x for x in list(sdata.tables['table'].obs.subtype.unique()) if x != celltype]),'gene_exp'] = 0
sdata.pl.render_shapes('cell_boundaries',color='gene_exp',cmap='Reds').pl.show()
plt.title(f'{brnum} {condition} {gene} in {celltype}')
plt.savefig(f'{folder}/{gene}_{celltype}_{brnum}_{condition}.png',dpi=500,bbox_inches='tight')

