import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.stats import pearsonr
from matplotlib.colors import LinearSegmentedColormap


# ED Fig. 10b

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n)))
    return new_cmap


minColor = .45
maxColor = 0
blue_t = truncate_colormap(plt.get_cmap("coolwarm"), minColor, maxColor)
red_t = truncate_colormap(plt.get_cmap("coolwarm_r"), minColor, maxColor)

data = pd.read_csv('pip_vs_pval_data.csv')
data = data.sort_values(['Genename','GWAS','Celltype'])
df = data.drop_duplicates(['Genename','GWAS','Celltype'])

genes = df.Genename.values
gwass = df.GWAS.values
celltypes = df.Celltype.values


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

plt.rcParams['font.size']='18'
fig, ax = plt.subplots(2,5,figsize=(18,8))
ax = ax.flatten()
for i in np.arange(len(genes)):
    gene = genes[i]
    gwas = gwass[i]
    celltype = celltypes[i]
    current_data = data[data.Genename==gene][data.GWAS==gwas][data.Celltype==celltype]
    r, _ = pearsonr(current_data.PIP, current_data.log_p)
    ax[i].scatter(current_data.PIP, current_data.log_p, c=current_data.R2, cmap=blue_t,label=None)
    sns.regplot(x=current_data.PIP, y=current_data.log_p, scatter=False, color='red', line_kws={'lw': 1}, ci=None,ax=ax[i])
    ax[i].text(min(current_data.PIP), max(current_data.log_p), f'r = {r:.3f}', fontsize=18, color='red', verticalalignment='top')
    if gene=='ELFN1':
        gene = 'MAD1L1'
    else:
        pass
    ax[i].set_title(f'{gene}\n {gwas} {celltype}',fontsize=18)
    ax[i].set_ylabel('-log(pval)')

plt.xlabel("PIP")

plt.tight_layout()
plt.savefig('/home/ah2428/review2_figs/pip_vs_pval.pdf',bbox_inches='tight',dpi=300)
plt.show()