import numpy as np
import pandas as pd
import json
import os
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns

from scipy.sparse import issparse
from pandas.api.types import is_numeric_dtype, is_categorical_dtype, is_list_like
from scipy.stats import zscore
from sklearn.metrics import adjusted_mutual_info_score
from natsort import natsorted

import anndata
import scanpy as sc
import pegasus as pg
from pegasusio import UnimodalData, MultimodalData
from typing import List, Tuple, Union, Optional, Callable
import logging
logger = logging.getLogger(__name__)
from pegasus.tools import X_from_rep, slicing, largest_variance_from_random_matrix

import spatialdata as sd
import spatialdata_plot
import spatialdata_io
from spatialdata import bounding_box_query



def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def deg_counts_plot(deg_ups,deg_downs,celltypes):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    
    minColor = .45
    maxColor = 0
    blue_t = truncate_colormap(plt.get_cmap("coolwarm"), minColor, maxColor)
    red_t = truncate_colormap(plt.get_cmap("coolwarm_r"), minColor, maxColor)

    degs_class_ptsd_count = pd.DataFrame({'Down':deg_downs,'Up':deg_ups})

    degs_class_ptsd = np.log(degs_class_ptsd_count)
    
    f = plt.figure(figsize=(5,5))
    plt.rcParams['font.size']='12'
    gs1 = GridSpec(3, 3)
    gs1.update(left=0, right=0.3, wspace=0)
    ax1 = plt.subplot(gs1[:, 1])
    ax2 = plt.subplot(gs1[:, 2])

    ax1.imshow(pd.DataFrame(degs_class_ptsd['Down']), cmap=blue_t, aspect='auto')
    ax2.imshow(pd.DataFrame(degs_class_ptsd['Up']), cmap=red_t, aspect='auto')

    ax1.set_yticks(np.arange(len(celltypes)),minor=False)
    ax1.set_yticklabels(celltypes)

    x_labels = ['Down','Up']
    ax1.set(xticks=np.arange(1),xticklabels=['Down'])
    ax2.set(xticks=np.arange(1),xticklabels=['Up'])
    ax2.set_yticks([])
    ax1.set_aspect(aspect=0.6)
    ax2.set_aspect(aspect=0.6)

    for i in range(len(celltypes)):
        color='k'
        ax1.text(ax1.get_xticks(), i, np.array(degs_class_ptsd_count['Down'])[i], ha="center", va="center", color=color)
        ax2.text(ax2.get_xticks(), i, np.array(degs_class_ptsd_count['Up'])[i], ha="center", va="center", color=color)

    ax1.tick_params(axis=u'both', which=u'both',length=0)
    ax2.tick_params(axis=u'both', which=u'both',length=0)

    minor_ticks = np.arange(0, len(celltypes)-1, 0.5)
    ax1.set_yticks(minor_ticks, minor=True)
    ax1.yaxis.grid(which='minor',color='k',linewidth=1)
    minor_ticks = np.arange(0.5, len(celltypes)-1, 1)
    ax2.set_yticks(minor_ticks, minor=True)
    ax2.yaxis.grid(which='minor',color='k',linewidth=1)

    ax1.axvline(x=0.5,color='k',linewidth=1)
    ax2.axvline(x=-0.5,color='k',linewidth=1)
    ax1.spines['right'].set_edgecolor('k')
    ax2.spines['left'].set_edgecolor('k')
    
    f.tight_layout()

def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def cluster_ngenes_ncells_across_samples_plot(data,cluster):
    df = data.obs[data.obs.anno==cluster].sort_values(['condition_merge','Channel'])
    counts = df[['Channel','condition_merge']].drop_duplicates('Channel').reset_index(drop=True)
    counts['ncells'] = df.groupby('Channel').size().reindex(df.Channel.unique()).values #reindex is correct
    counts['ngenes'] = df.groupby('Channel').n_genes.sum().reindex(df.Channel.unique()).values
    total = counts.groupby('condition_merge').agg(sum)
    total.columns = ['ncells_total','ngenes_total']
    counts = pd.merge(counts,total,on='condition_merge')
    counts['ngenes_percent'] = round((counts.ngenes/counts.ngenes_total)*100,2)
    counts = counts.set_index('Channel')
    if len(counts)<92:
        missing = list(set(data.obs.Channel.unique()) - set(df.Channel.unique()))
        for x in missing:
            counts.loc[x] = 0
    fig,ax = plt.subplots(2,1,figsize=(30,10),gridspec_kw={'hspace': 0.7})
    ax[0]= df.boxplot(column='n_genes',by=['condition_merge','Channel'],rot=90,fontsize=14,grid=False,showfliers=False,patch_artist=True,ax=ax[0])
    ax[0].axhline(df.n_genes.mean(), c='r')
    ax[0].set_ylabel('n_genes',fontsize=14)
    ax[0].set_xticks(np.arange(1,93))
    xlabels = [counts.index[i]+' '+str(counts.ngenes_percent[i])+'%' for i in range(len(counts))]
    ax[0].set_xticklabels(xlabels)
    ax[0].set_title('')
    ax[1] = counts.ncells.plot.bar(ax=ax[1],fontsize=14)
    ax[1].set_ylabel('n_cells',fontsize=14)
    ax[1].set_xticklabels(counts.index)
    for j in range(2):
        ax[j].set_xlabel('')
        for i in range(92):
            if counts.condition_merge[i]=='CON':
                ax[j].findobj(matplotlib.patches.Patch)[i].set_facecolor("red")
            elif counts.condition_merge[i]=='MDD':
                ax[j].findobj(matplotlib.patches.Patch)[i].set_facecolor("green")
            elif counts.condition_merge[i]=='PTSD':
                ax[j].findobj(matplotlib.patches.Patch)[i].set_facecolor("blue")
            else:
                pass
    con = matplotlib.patches.Patch(color='r', label='CON')
    mdd = matplotlib.patches.Patch(color='g', label='MDD')
    ptsd = matplotlib.patches.Patch(color='b', label='PTSD')
    plt.legend(handles=[con,mdd,ptsd])
    plt.suptitle('{} n_genes and n_cells across samples'.format(cluster),fontsize=30)
    
def plot_TF_footprint(motif,motif_name,celltype,data):
    con = data[data[data.columns[-1]]==f'{celltype}_CON']
    ptsd = data[data[data.columns[-1]]==f'{celltype}_PTSD']

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.size']= '20'
    plt.rc('axes', unicode_minus=False)

    fig,ax = plt.subplots(1,1,figsize=(7,5))
    ax.plot(con[motif_name+'.x'],con[motif_name+'.mean'],color='gray',label=f'CON',linewidth=2)
    ax.plot(ptsd[motif_name+'.x'],ptsd[motif_name+'.mean'],color='red',label=f'PTSD',linewidth=2)
    ax.set_xlim([-250,250])
    ax.set_xticks([-250,0,250])
    #ax.set_ylim([-0.1,data[f'{motif_name}.mean'].max().max()+0.05])
    ax.set_ylim([-0.5,1.3])
    ax.set_yticks([0,0.5,1])
    ax.set_ylabel('Tn5 Bias Subtracted\n Normalized Insertions')
    ax.set_xlabel('Distance to motif center (bp)')
    ax.set_title(motif,x=0.2,y=0.85)
    plt.legend(frameon=False,loc='upper right',fontsize=16)

def plot_hiC_track(chrom,res,a,b,tad1,tad2):
    matrix_object = hic_pos.getMatrixZoomData(chrom, chrom, "observed", "KR", "BP", res)
    numpy_matrix = matrix_object.getRecordsAsMatrix(a,b,a,b)
    tad = matrix_object.getRecordsAsMatrix(tad1,tad2,tad1,tad2)

    fig, ax = plt.subplots()
    ax.imshow(numpy_matrix,cmap='Greys',vmin=0,vmax=5)
    
    length = int((b-a)/res)
    tri = np.triu(numpy_matrix)
    for i in range(length):
        tri[i,int(length * (2/5))+i:] = 0
    plt.matshow(tri, cmap='Greys', vmin=0, vmax=10)
    plt.axis('off')

    x = (b-a)/res
    x1 = int(length * ((tad1-a)/(b-a)))
    x2 = int(length * ((tad2-a)/(b-a)))

    plt.hlines(y=x1, xmin=x1, xmax=x2, color='r', alpha=0.5)
    plt.vlines(x=x2, ymin=x1, ymax=x2, color='r', alpha=0.5)

def snXenium_violin(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    attrs: Union[str, List[str]],
    groupby: str,
    hue: Optional[str] = None,
    matkey: Optional[str] = None,
    stripplot: Optional[bool] = False,
    stripsize: int = 1,
    inner: Optional[str] = None,
    scale: Optional[str] = 'width',
    panel_size: Optional[Tuple[float, float]] = (8, 0.5),
    palette: Optional[List[str]] = None,
    left: Optional[float] = 0.15,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.1,
    ylabel: Optional[str] = None,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:

    if not is_list_like(attrs):
        attrs = [attrs]

    if not isinstance(data, anndata.AnnData):
        cur_matkey = data.current_matrix()
    if matkey is not None:
        assert not isinstance(data, anndata.AnnData)
        data.select_matrix(matkey)

    # Filter out attributes not existing in the data
    #attrs = _get_valid_attrs(data, attrs)

    nrows = len(attrs)
    #fig, axes = _get_subplot_layouts(nrows=nrows, ncols=1, panel_size=panel_size, dpi=dpi, left=left, bottom=bottom, wspace=wspace, hspace=0, squeeze=False, sharey=False)

    obs_keys = []
    genes = []

    for key in attrs:
        if key in data.obs:
            assert is_numeric_dtype(data.obs[key])
            obs_keys.append(key)
        else:
            genes.append(key)

    df_list = [pd.DataFrame({"label": data.obs[groupby].values})]

    if hue is not None:
        df_list.append(pd.DataFrame({hue: data.obs[hue].values}))
        stripplot = False
        kwargs['hue'] = hue
        kwargs['split'] = True
    else:
        kwargs['hue'] = "label"
        kwargs['legend'] = False
        kwargs['split'] = False

    if len(obs_keys) > 0:
        df_list.append(data.obs[obs_keys].reset_index(drop=True))
    if len(genes) > 0:
        expr_mat = slicing(data[:, genes].X)
        df_list.append(pd.DataFrame(data=expr_mat, columns=genes))
    df = pd.concat(df_list, axis = 1)
    return df

def plot_xenium_slide(brnum,gene,color,savepath):
    sdata = spatialdata_io.xenium(f'/gpfs/gibbs/pi/girgenti/ah2428/xenium/data/{brnum}_resegment/outs/')
    data = sc.read_h5ad('/gpfs/gibbs/pi/girgenti/ah2428/xenium/data_processed_resegment/snXenium.h5ad')
    condition = data.obs[data.obs.Channel==int(brnum)].Condition.unique()[0]
    
    sc.pp.normalize_total(sdata.tables["table"])
    sc.pp.log1p(sdata.tables["table"])
    sc.pp.highly_variable_genes(sdata.tables["table"])
    sdata.tables["table"].var.sort_values("means")
    sdata.tables["table"].obs["region"] = "cell_boundaries"
    sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")
    
    sdata.tables["table"].obs['cellid'] = [x.split('-')[0] for x in sdata.tables["table"].obs['cell_id'].values]
    df = data.obs[data.obs.cellid.isin(sdata.tables["table"].obs.cellid)][['cellid','celltype']]
    subtypes = []
    for cell in sdata.tables['table'].obs['cellid'].values:
        if cell in df.cellid.values:
            subtype = df[df.cellid==cell].celltype.values[0]
        else:
            subtype = 'NA'
        subtypes.append(subtype)
    sdata.tables["table"].obs['celltype'] = subtypes
    sdata.tables['table'].obs['celltype'] = sdata.tables['table'].obs['celltype'].astype('category')
    
    cells = sdata.tables['table'].obs[sdata.tables['table'].obs.celltype.isin(['EXC','INH','OLI','OPC','END','AST','MIC'])].cell_labels.values
    cells = [int(c)-1 for c in cells]
    sdata_new = sdata.tables['table'][cells,:]
    sdata.tables['table'] = sdata_new
    cellids = sdata.tables['table'].obs.cell_id.values
    cell_boundaries = sdata.shapes['cell_boundaries'] 
    cell_circles = sdata.shapes['cell_circles']
    nucleus_boundaries = sdata.shapes['nucleus_boundaries']
    sdata.shapes['cell_boundaries'] = cell_boundaries[cell_boundaries.index.isin(cellids)]
    sdata.shapes['cell_circles'] = cell_circles[cell_circles.index.isin(cellids)]
    sdata.shapes['nucleus_boundaries'] = nucleus_boundaries[nucleus_boundaries.index.isin(cellids)]
    
    sc.pp.normalize_total(sdata.tables["table"])
    sc.pp.log1p(sdata.tables["table"])
    sc.pp.highly_variable_genes(sdata.tables["table"])
    sdata.tables["table"].var.sort_values("means")
    sdata.tables["table"].obs["region"] = "cell_boundaries"
    sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")
    
    sdata.tables['table'].obs.celltype = sdata.tables['table'].obs.celltype.cat.reorder_categories(['EXC','INH','OLI','OPC','END','AST','MIC'])
    
    sdata.pl.render_shapes("cell_boundaries",
                           color='celltype',
                           groups=['EXC','INH','OLI','OPC','END','AST','MIC'],
                           palette=['#b22222', '#2E8B57', '#5254a3', '#aec7e8', '#bc80bd', '#ffed6f', '#7f7f7f'],
                           scale=5, 
                           outline=False, 
                           fill_alpha=0.1).pl.render_points('transcripts',color='feature_name',groups=gene,palette=color,size=1).pl.show(figsize=(10,20))
    
    plt.legend('',frameon=False)
    plt.title('')
    plt.axis('off')
    plt.savefig(f'{savepath}/{brnum}_{gene}_xenium_slide.pdf',bbox_inches='tight',dpi=1000)
    plt.show()
    return sdata

def plot_xenium_gene_exp_barplot(gene):
    data = sc.read_h5ad('/gpfs/gibbs/pi/girgenti/ah2428/xenium/data_processed_resegment/snXenium.h5ad')
    data_pg = MultimodalData(data)
    pg.log_norm(data_pg,base_matrix='raw')
    
    con_cells = data.obs[data.obs.Condition=='CON'].index.values
    ptsd_cells = data.obs[data.obs.Condition=='PTSD'].index.values
    data_ptsd_con = data[np.append(con_cells,ptsd_cells),:]
    
    df = snXenium_violin(data_ptsd_con,attrs=gene,groupby='celltype',hue='Condition')
    df['label'] = df['label'].astype('category')
    df['label'] = df['label'].cat.reorder_categories(['EXC','INH','OLI','OPC','END','AST','MIC'])
    x = pd.DataFrame(df.groupby(['label','Condition'])[gene].mean())
    x = x.loc[['EXC','INH','OLI','OPC','END','AST','MIC']]
    x = x.reset_index()
    
    plt.rcParams['font.size']='16'
    fig, ax = plt.subplots(figsize=(8,2))
    ax.bar(np.arange(7),x[x.Condition=='CON'][gene].values,width=0.4,label='CON',color='lightgray')
    ax.bar(np.arange(7)+0.4,x[x.Condition=='PTSD'][gene].values,width=0.4,label='PTSD',color='tab:blue')
    ax.set_xticks(np.arange(7)+0.2,['EXN','IN','OLG','OPC','END','AST','MG'])
    plt.legend(frameon=False)
    plt.ylabel(f'{gene} exp')
    plt.show()

def plot_colortable(colors, *, ncols=4, sort_colors=True):
    if sort_colors is True:
        names = sorted(
            colors, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c))))
    else:
        names = list(colors)
    return names

def create_metadata_df(meta,colname,name):
    arr = []
    for s in meta[colname].cat.categories:
        vals = list(meta[meta[colname]==s][name].value_counts()[sorted(meta[name].unique())]/meta[meta[colname]==s][name].value_counts()[sorted(meta[name].unique())].sum())
        arr.append(vals)
    df = pd.DataFrame(arr)
    df.insert(0, colname, meta[colname].cat.categories, True)
    df.columns = [colname]+sorted(meta[name].unique())
    return df 

def create_metadata_df_bins(meta,colname,name):
    arr = []
    for s in meta[colname].cat.categories:
        temp = pd.cut(meta[meta[colname]==s][name],bins=10)
        vals = list(temp.value_counts()[sorted(temp.cat.categories)]/meta[meta[colname]==s][name].value_counts().sum())
        arr.append(vals)
    df = pd.DataFrame(arr)
    df.insert(0, colname, meta[colname].cat.categories, True)
    df.columns = [colname]+sorted(temp.cat.categories)
    return df 

def metadata_barplot(df,colname,cmap,savefile):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    if cmap is not None:
        df.plot(colname,kind='bar',stacked=True,rot=45,fontsize=16,figsize=(7,5),legend=False,cmap=cmap)
    else:
        df.plot(colname,kind='bar',stacked=True,rot=45,fontsize=16,figsize=(7,5),legend=False)
    plt.xticks(np.arange(len(df[colname].values)), df[colname].values, ha='right',rotation_mode='anchor')
    plt.ylabel('Proportion',fontsize=16)
    plt.xlabel(' ')
    plt.savefig(savefile,bbox_inches='tight',dpi=300)



