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

args = commandArgs(trailingOnly=TRUE)

cluster <- args[1] # celltype
cond <- args[2] # condition (PTSD,MDD)

folder <- glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/RNA/results/062522_DEG/subclass/wilcox/{cond}_vs_CON')
data <- LoadH5Seurat(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/RNA/data/DEG_data_subclass_062522/{cond}_vs_CON/{cond}_CON_{cluster}.h5seurat'),assays='RNA')
meta <- read.delim2(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/RNA/data/DEG_data_subclass_062522/{cond}_vs_CON/{cond}_CON_{cluster}_meta.txt'))

# if exceed R memory 2^31
if (dim(data)[2] > 76745) {
    print('True')
    mat <- data$RNA@counts
    tmp <- matrix(data = 0, nrow=nrow(mat), ncol=ncol(mat))
    row_pos <- mat@i+1
    col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
    val <- mat@x
    for (i in seq_along(val)){
        tmp[row_pos[i],col_pos[i]] <- val[i]
    }
    rownames(tmp) <- rownames(mat)
    colnames(tmp) <- colnames(mat)
    dge <- DGEList(counts = tmp)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- cpm(dge)
} else {
    dge <- DGEList(counts = data$RNA@counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- cpm(dge)
}


num_con <- floor(dim(meta[meta$Condition=='CON',])[1]*0.05)
num_cond <- floor(dim(meta[meta$Condition==cond,])[1]*0.05)

counts <- data$RNA@counts
if (dim(data)[2] > 76745) {
    bin_counts <- as.matrix((tmp > 0) + 0)
} else {
    bin_counts <- as.matrix((counts > 0) + 0)
}
con_sums <- rowSums(bin_counts[, which(meta$Condition=='CON')])
cond_sums <- rowSums(bin_counts[, which(meta$Condition==cond)])
df <- data.frame(con=con_sums,cond=cond_sums)
filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
data <- data[filtered_genes,]

data$condition <- meta$Condition
Idents(data) <- 'condition'

markers <- FindMarkers(data, ident.1=glue('{cond}'), ident.2='CON', slot='data', test.use='wilcox', logfc.threshold=0, min.pct=0)
sig_deg = markers[ which(markers$p_val_adj < 0.01 & abs(markers$avg_log2FC) > log2(1.2)), ]
write.table(markers,glue('{folder}/{cluster}_DEG.csv'),sep='\t')
write.table(sig_deg,glue('{folder}/{cluster}_SIG_DEG.csv'),sep='\t')

mat <- cpms

data@meta.data$condition <- meta$Condition
df <- data@meta.data['condition']
rownames(df) <- 1:nrow(df)
name2 = 'avg_log2FC'

con_ind <- which(df$condition=='CON')
cond_ind <- which(df$condition==cond)

# check consistency of direction 
print('check consistency')

sig_down <- sig_deg[sig_deg[[name2]] < 0,]
sig_up <- sig_deg[sig_deg[[name2]] > 0,]

count <- 0
for (i in 1:dim(sig_up)[1]) {
    gene <- rownames(sig_up)[i]
    df[[i]] <- mat[gene,]
    con_mean <- mean(df[[i]][con_ind])
    cond_mean <- mean(df[[i]][cond_ind])
    if ((sig_up[rownames(sig_up)==gene,][[name2]] > 0) & (cond_mean - con_mean > 0)) {
        count <- count + 1
    }
}
consistent_up <- round(count/(dim(sig_up)[1]),3)
print(glue('Consistent up: ', {consistent_up}))

count <- 0
for (i in 1:dim(sig_down)[1]) {
    gene <- rownames(sig_down)[i]
    df[[i]] <- mat[gene,]
    con_mean <- mean(df[[i]][con_ind])
    cond_mean <- mean(df[[i]][cond_ind])
    if ((sig_down[rownames(sig_down)==gene,][[name2]] < 0) & (cond_mean - con_mean < 0)) {
        count <- count + 1
    }
}
consistent_down <- round(count/(dim(sig_down)[1]),3)
print(glue('Consistent down: ', {consistent_down}))

EnhancedVolcano(markers,lab=rownames(markers),x='avg_log2FC',y='p_val_adj',pCutoff=0.01,title=glue("{cluster} {cond} vs CON Wilcox"),labSize=4,FCcutoff=log2(1.2),caption=glue("genes in at least 5% of cells = ", nrow(markers), " sig degs = ", nrow(sig_deg)),legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))

ggsave(glue('{folder}/{cluster}_volcano.png'))
