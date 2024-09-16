import pegasus as pg
import numpy as np
import scanpy as sc
import importlib
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

import sccoda
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
from sccoda.util import *


def run_scCODA(meta,channel_colname,celltype_colname):
    samples = meta[channel_colname].unique()
    celltypes = meta[celltype_colname].unique()
    counts = np.zeros((len(samples),len(celltypes)))
    for ii,s in enumerate(samples):
        for jj,c in enumerate(celltypes):
            count = meta[meta[channel_colname]==s][meta[celltype_colname]==c].shape[0]
            counts[ii,jj] = int(count)
    counts = pd.DataFrame(counts,index=samples,columns=celltypes)
    counts = counts.astype(int)
    cov = pd.DataFrame(meta[[channel_colname,'Condition']])
    cov = cov.drop_duplicates(channel_colname)
    cov = cov.set_index(channel_colname)
    counts['Condition'] = cov['Condition'].values
    counts.Condition = counts.Condition.astype('category')
    counts.Condition = counts.Condition.cat.reorder_categories(['PTSD','CON','MDD'])
    sccoda_data = sccoda.util.cell_composition_data.from_pandas(counts,covariate_columns=['Condition'])

    print('PTSD vs CON, MDD vs CON')
    model = mod.CompositionalAnalysis(sccoda_data, formula="C(Condition, Treatment('CON'))", reference_cell_type="automatic")
    sim_results = model.sample_hmc()
    sim_results.set_fdr(est_fdr=0.05)
    sim_results.summary()

    print('PTSD vs MDD, CON vs MDD')
    model = mod.CompositionalAnalysis(sccoda_data, formula="C(Condition, Treatment('MDD'))", reference_cell_type="automatic")
    sim_results = model.sample_hmc()
    sim_results.set_fdr(est_fdr=0.05)
    sim_results.summary()

# run scCODA

rna = pd.read_csv('/home/ah2428/ShareZhangLab/PTSD/RNA/data/RNA_FINAL_meta.csv',sep='\t')
run_scCODA(rna,'Channel','subclass')

atac = pd.read_csv('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL_meta.csv',sep='\t')
run_scCODA(atac,'Sample','Int_Cluster')

multiome = pd.read_csv('/home/ah2428/ShareZhangLab/PTSD/Multiome/data/MULTIOME_FINAL_meta.csv',sep='\t')
run_scCODA(multiome,'Sample','subtype')

xenium = pd.read_csv('/home/ah2428/ShareZhangLab/PTSD/xenium/data/XENIUM_FINAL_meta.csv',sep='\t')
run_scCODA(xenium,'Sample','celltype')


def get_cell_proportions(meta,channel_colname,celltype_colname,celltypes):
    props = []
    for condition in ['CON','MDD','PTSD']:
        data_condition = meta[meta.Condition==condition]
        n_channels = len(data_condition[channel_colname].unique())
        prop = np.zeros((n_channels,len(celltypes)))
        for i,channel in enumerate(data_condition[channel_colname].unique()):
            data_channel = data_condition[data_condition[channel_colname]==channel]
            for j,celltype in enumerate(celltypes):
                data_celltype = data_channel[data_channel[celltype_colname]==celltype]
                prop[i,j] = data_celltype.shape[0]/data_channel.shape[0]
        props.append(prop)
    return props

rna_props = get_cell_proportions(rna,'Channel','subclass',['CUX2','RORB','FEZF2','OPRK1','LAMP5','KCNG1','VIP','SST','PVALB','OLI','OPC','END','AST','MIC'])
atac_props = get_cell_proportions(atac,'Sample','Int_Cluster',['Ex','In','Oli','OPC','End','Ast','Mic'])
multiome_props = get_cell_proportions(multiome,'Sample','subtype',['CUX2','RORB','FEZF2','OPRK1','LAMP5','KCNG1','VIP','SST','PVALB','OLI','OPC','END','AST','MIC'])
xenium_props = get_cell_proportions(xenium,'Sample','celltype',['EXN','IN','OLG','OPC','END','AST','MG'])

# scCODA plot

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size']='16'

fig, ax = plt.subplots(2,2,figsize=(16,10),sharey=True)
ax = ax.flatten()

bp1 = ax[0].boxplot(rna_props[0],positions=np.arange(1,len(subtypes)*2,2),vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:gray'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp2 = ax[0].boxplot(rna_props[1],positions=np.arange(1,len(subtypes)*2,2)+0.5,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:orange'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp3 = ax[0].boxplot(rna_props[2],positions=np.arange(1,len(subtypes)*2,2)+1,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:blue'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
ax[0].set_xticks(np.arange(1.5,len(subtypes)*2+1.5,2),subtypes,rotation=45)
ax[0].legend([bp1["boxes"][0], bp2["boxes"][0], bp3["boxes"][0]],['CON','MDD','PTSD'],frameon=False)
ax[0].set_title('RNA (38 CON, 32 MDD, 35 PTSD)')
ax[0].text(7,0.53,'** FDR < 0.05',color='blue')
ax[0].get_xticklabels()[4].set_weight("bold")
ax[0].axhline(0.55,19.1/28,19.7/28,color='blue')
ax[0].text(19.55,0.55,'**',color='blue',fontsize='small')
ax[0].axhline(0.6,18.3/28,19.7/28,color='blue')
ax[0].text(19.1,0.6,'**',color='blue',fontsize='small')
ax[0].axhline(0.55,24.4/28,25/28,color='blue')

bp1 = ax[1].boxplot(atac_props[0],positions=np.arange(1,len(celltypes)*2,2),vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:gray'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp2 = ax[1].boxplot(atac_props[1],positions=np.arange(1,len(celltypes)*2,2)+0.5,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:orange'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp3 = ax[1].boxplot(atac_props[2],positions=np.arange(1,len(celltypes)*2,2)+1,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:blue'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
ax[1].set_xticks(np.arange(1.5,len(celltypes)*2+1.5,2),celltypes,rotation=45)
ax[1].set_title('ATAC (34 CON, 31 MDD, 29 PTSD)')
ax[1].get_xticklabels()[3].set_weight("bold")

bp1 = ax[2].boxplot(multiome_props[0],positions=np.arange(1,len(subtypes)*2,2),vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:gray'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp2 = ax[2].boxplot(multiome_props[1],positions=np.arange(1,len(subtypes)*2,2)+0.5,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:orange'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp3 = ax[2].boxplot(multiome_props[2],positions=np.arange(1,len(subtypes)*2,2)+1,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:blue'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
ax[2].set_xticks(np.arange(1.5,len(subtypes)*2+1.5,2),subtypes,rotation=45)
ax[2].set_title('Multiome (9 CON, 4 MDD, 12 PTSD)')
ax[2].get_xticklabels()[4].set_weight("bold")

bp1 = ax[3].boxplot(xenium_props[0],positions=np.arange(1,len(celltypes)*2,2),vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:gray'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp2 = ax[3].boxplot(xenium_props[1],positions=np.arange(1,len(celltypes)*2,2)+0.5,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:orange'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
bp3 = ax[3].boxplot(xenium_props[2],positions=np.arange(1,len(celltypes)*2,2)+1,vert=True,patch_artist=True,showfliers=True,boxprops=dict(facecolor='tab:blue'),medianprops=dict(color='k'),flierprops={'marker':'.','markersize':1,'markerfacecolor':'k'})
ax[3].set_xticks(np.arange(1.5,len(celltypes)*2+1.5,2),celltypes,rotation=45)
ax[3].set_title('Xenium (10 CON, 4, MDD, 4 PTSD)')
ax[3].get_xticklabels()[5].set_weight("bold")

fig.supylabel('Proportion')
fig.supxlabel('Cell type')
fig.tight_layout()

for i, label in enumerate(('a', 'b', 'c', 'd')):
    ax[i].text(0.02, 1.15, label, transform=ax[i].transAxes, fontweight='bold', va='top', ha='right', fontsize=20)

plt.savefig('/home/ah2428/palmer_scratch/figs/r2r/scCODA_plot.png',bbox_inches='tight',dpi=300)


# scCODA table

data = ['RNA','RNA']
fdr = ['<0.05','<0.05']
condition = ['PTSD vs CON','PTSD vs MDD']
celltype = ['OLG','OLG']
log2FC = [-0.273407,-0.282037]

df = pd.DataFrame({'Data':data,'FDR':fdr,'Condition':condition,'Celltype':celltype,'log2FC':log2FC})
dfi.export(df.style.hide(axis='index'), '/home/ah2428/palmer_scratch/figs/r2r/scCODA_table.png',table_conversion = 'matplotlib',dpi=300)

