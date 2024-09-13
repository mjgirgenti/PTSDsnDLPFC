import pandas as pd
import numpy as np

p2g = pd.read_csv('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/p2glinks_maxDist_2M.csv',sep='\t')

enh = p2g[p2g.Correlation > 0.4][p2g.FDR < 0.05]
num_links = []
for ii,g in enumerate(enh.geneName.unique()):
    df = enh[enh.geneName==g]
    num_link = len(df.peakName.unique())
    num_links.append(num_link)
    if ii % 1000 == 0:
        print(ii)

pd.DataFrame(num_links).to_csv('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/nCRE_for_hist.csv',header=None,index=None)
    
