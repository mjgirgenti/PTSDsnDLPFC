library(MAST)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(EnhancedVolcano)
library(glue)
library(data.table)
library(stringr)
library(edgeR)
library(dplyr)
library(DESeq2)
library(R.utils)
library(nebula)
#library(argparse)

args = commandArgs(trailingOnly=TRUE)
#parser <- ArgumentParser()
#parser$add_argument('--celltype', type='character', help='celltype')
#parser$add_argument('--method', type='character')
#parser$add_argument('--model', type='character')
#parser$add_argument('--kappa', type='character')
#parser$add_arguemnt('--cell', type='character')
#args <- parser$parse_args()
#print('here')
celltype <- args[1]
method <- args[2]
model <- args[3]
kappa <- as.integer(args[4])
cell <- as.integer(args[5])

cond <- 'PTSD'
data <- LoadH5Seurat(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{cond}_vs_CON/{cond}_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
meta <- read.delim2(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{cond}_vs_CON/{cond}_CON_{celltype}_meta.txt'))

protein <- read.table('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt')$V2
protein_genes <- intersect(rownames(data),protein)
data <- data[protein_genes,]
print(dim(data))

counts <- data$RNA@counts
bin_counts <- as.matrix((counts > 0) + 0)
num_con <- floor(dim(meta[meta$Condition=='CON',])[1]*0.05)
num_cond <- floor(dim(meta[meta$Condition==cond,])[1]*0.05)
con_sums <- rowSums(bin_counts[, which(meta$Condition=='CON')])
cond_sums <- rowSums(bin_counts[, which(meta$Condition==cond)])
df <- data.frame(con=con_sums,cond=cond_sums)
filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
data2 <- data[filtered_genes,]

meta2 = meta[meta$barcodekey %in% colnames(data2),]
meta2 <- meta2[order(meta$Channel),]
data2@meta.data$Sample <- meta2$Channel
data2 <- data2[,meta2$barcodekey]
meta2 <- meta2[,c('AgeDeath','Sex','Race','PMI','RIN','Condition')]
cdr <- colSums(data2$RNA@counts>0)

meta2$cdr <- as.numeric(cdr)
meta2$AgeDeath <- as.numeric(meta2$AgeDeath)
meta2$RIN <- as.numeric(meta2$RIN)
meta2$PMI <- as.numeric(meta2$PMI)
meta2$Race <- as.factor(meta2$Race)
meta2$Sex <- as.factor(meta2$Sex)
meta2$Condition <- as.factor(meta2$Condition)

df = model.matrix(~cdr+AgeDeath+Sex+Race+PMI+RIN+Condition, data=meta2)
re = nebula(data2$RNA@counts,data2@meta.data$Sample,pred=df,method=method,model=model,kappa=kappa,cutoff_cell=cell,ncore=5)
deg <- re$summary
rownames(deg) <- rownames(data2$RNA@counts)

deg$padj <- p.adjust(deg$p_ConditionPTSD,method='BH')
write.table(deg,glue('/home/ah2428/palmer_scratch/nebula/{celltype}_{method}_{model}_{kappa}_{cell}.csv'),sep='\t')

sig_deg = deg[ which(deg$padj < 0.2 & abs(deg$logFC_ConditionPTSD) > log2(1.2)), ]

num_up <- nrow(sig_deg[sig_deg$logFC_ConditionPTSD > 0,])
num_down <- nrow(sig_deg[sig_deg$logFC_ConditionPTSD < 0,])

EnhancedVolcano(deg,subtitle=NULL,
                lab=rownames(deg),x='logFC_ConditionPTSD',y='padj',pCutoff=0.2,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", nrow(deg), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
ggsave(glue('/home/ah2428/palmer_scratch/nebula/{celltype}_{method}_{model}_{kappa}_{cell}.png'))

