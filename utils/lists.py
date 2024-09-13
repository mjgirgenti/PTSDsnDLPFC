import pandas as pd

lake_genes = ['MRC1','TMEM119','CX3CR1','APBB1IP','CLDN5','FLT1','DUSP1','COBLL1','GLUL','SOX9','AQP4','GJA1','NDRG2','GFAP','ALDH1A1','ALDH1L1','VIM','SLC4A4','PDGFRA','PCDH15','OLIG2','OLIG1','PLP1','MAG','MOG','MOBP','MBP','SATB2','SLC17A7','GRM4','GAD2','GAD1','SLC32A1','SST','PVALB','STMN2','RBFOX3']

celltypes = ['EXN','IN','OLG','OPC','END','AST','MG']

subtypes = ['CUX2','RORB','FEZF2','OPRK1','LAMP5','KCNG1','VIP','SST','PVALB','OLG','OPC','END','AST','MG']

exn_subtypes = ['CUX2','RORB','FEZF2','OPRK1']

in_subtypes = ['LAMP5','KCNG1','VIP','SST','PVALB']

protein_coding = pd.read_csv('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/data/protein_coding.genes.with.chr.txt',sep='\t',header=None)[1].values

celltype_colors = {
    'CUX2': '#ad494a',
    'RORB': '#e7969c',
    'FEZF2': '#843c39',
    'OPRK1': '#d6616b',
    'LAMP5': '#637939',
    'KCNG1': '#a6d854',
    'VIP': '#cedb9c',
    'SST': '#8ca252',
    'PVALB': '#b5cf6b',
    'OLI': '#5254a3',
    'OPC': '#aec7e8',
    'END': '#bc80bd',
    'AST': '#ffed6f',
    'MIC': '#7f7f7f',
    'EXC':'#b22222',
    'INH':'#2E8B57'}

celltype_colors2 = {'EXN':'#b22222',
                'IN':'#2E8B57',
                'OLG':'#5254a3',
                'OPC':'#aec7e8',
                'END':'#bc80bd',
                'AST':'#ffed6f',
                'MG':'#7f7f7f'}



