import numpy as np
import pandas as pd
import pegasus as pg
import json
import os
from pegasusio import UnimodalData, MultimodalData
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import hicstraw

def plot_colortable(colors, *, ncols=4, sort_colors=True):
    # Sort colors by hue, saturation, value and name.
    if sort_colors is True:
        names = sorted(
            colors, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c))))
    else:
        names = list(colors)
    return names

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def DEG_counts_plot(deg_ups,deg_downs,celltypes,title):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

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

def plot_hiC_track(hic_data,chrom,res,a,b,tad1,tad2):
    matrix_object = hic_data.getMatrixZoomData(chrom, chrom, "observed", "KR", "BP", res)
    numpy_matrix = matrix_object.getRecordsAsMatrix(a,b,a,b)
    tad = matrix_object.getRecordsAsMatrix(tad1,tad2,tad1,tad2)

    #fig, ax = plt.subplots()
    #ax.imshow(numpy_matrix,cmap='Greys',vmin=0,vmax=5)
    
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

def create_metadata_df(meta,name,colname):
    arr = []
    for s in meta[colname].cat.categories:
        print(s)
        vals = list(meta[meta[colname]==s][name].value_counts()[sorted(meta[name].unique())]/meta[meta[colname]==s][name].value_counts()[sorted(meta[name].unique())].sum())
        arr.append(vals)
    df = pd.DataFrame(arr)
    df.insert(0, "subclass", meta[colname].cat.categories, True)
    df.columns = ['subclass']+sorted(meta[name].unique())
    return df 

def create_metadata_df_bins(meta,name,colname):
    arr = []
    for s in meta[colname].cat.categories:
        print(s)
        temp = pd.cut(meta[meta[colname]==s][name],bins=10)
        vals = list(temp.value_counts()[sorted(temp.cat.categories)]/meta[meta[colname]==s][name].value_counts().sum())
        arr.append(vals)
    df = pd.DataFrame(arr)
    df.insert(0, "subclass", meta[colname].cat.categories, True)
    df.columns = ['subclass']+sorted(temp.cat.categories)
    return df 

def metadata_barplot(meta,df,savefile,cmap,colname):
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    #df.plot('subclass',kind='bar',stacked=True,rot=45,fontsize=16,figsize=(7,5),legend=False)
    if cmap is not None:
        df.plot('subclass',kind='bar',stacked=True,rot=45,fontsize=16,figsize=(7,5),legend=False,cmap=cmap)
    else:
        df.plot('subclass',kind='bar',stacked=True,rot=45,fontsize=16,figsize=(7,5),legend=False)
    plt.xticks(np.arange(len(meta[colname].cat.categories)), meta[colname].cat.categories, ha='right',rotation_mode='anchor')
    plt.ylabel('Proportion',fontsize=16)
    plt.xlabel(' ')
    plt.savefig(savefile,bbox_inches='tight',dpi=300)
    
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