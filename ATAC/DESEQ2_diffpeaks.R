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
library(ArchR)
library(ArchRtoSignac)
library(EnsDb.Hsapiens.v86)

addArchRThreads(threads = 10)

args = commandArgs(trailingOnly=TRUE)
celltype <- args[1]

data <- LoadH5Seurat(glue('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/celltype_condition/{celltype}_signac.h5seurat'))
DefaultAssay(data) <- 'peaks'
Idents(data) <- data@meta.data$Cluster_Cond

peaks <- read.table('/home/ah2428/palmer_scratch/celltype_condition_union_peaks.csv',sep='\t')
celltype_peaks <- peaks[peaks$celltype %in% c(glue('{celltype}_PTSD'),glue('{celltype}_CON')),]$peak
data2 <- data[celltype_peaks,]
counts <- data2$peaks@counts
bin_counts <- as.matrix((counts > 0) + 0)

meta <- data2@meta.data
#cond <- 'PTSD'
#num_con <- floor(dim(meta[meta$Condition=='CON',])[1]*0.05)
#num_cond <- floor(dim(meta[meta$Condition==cond,])[1]*0.05)
#con_sums <- rowSums(bin_counts[, which(meta$Condition=='CON')])
#cond_sums <- rowSums(bin_counts[, which(meta$Condition==cond)])
#df <- data.frame(con=con_sums,cond=cond_sums)
#filtered_peaks <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
#data2 <- data2[filtered_peaks,]

counts <- data.frame(data2$peaks@counts)
strings <- colnames(counts)
new_colnames <- sapply(strsplit(strings,'[_]'),`[`,1)
colnames(counts) <- new_colnames
counts_df <- sapply(split.default(counts, colnames(counts)), rowSums)
counts_df <- counts_df+1

meta_new <- meta[,c('Sample','Condition')]
meta_new <- meta_new[!duplicated(meta_new[ , c('Sample')]),]
rownames(meta_new) <- meta_new$Sample
meta_new <- meta_new[order(meta_new$Sample),]

dds <- DESeqDataSetFromMatrix(countData=counts_df, 
                              colData=meta_new, 
                              design=~Condition)
dds <- DESeq(dds)
res <- results(dds)
write.table(res,glue('/home/ah2428/palmer_scratch/diffpeaks/{celltype}_diffpeaks.csv'),sep='\t')
sig_deg = res[ which(res$pvalue < 0.05 & abs(res$log2FoldChange) > log2(1.2)), ]
num_up <- nrow(sig_deg[sig_deg$log2FoldChange > 0,])
num_down <- nrow(sig_deg[sig_deg$log2FoldChange < 0,])

EnhancedVolcano(res,subtitle=NULL,title=glue("{celltype} DESEQ2 diffpeaks"),
                lab=rownames(res),x='log2FoldChange',y='pvalue',pCutoff=0.05,labSize=4,FCcutoff=log2(1.2),
                caption=glue("TOTAL = ", nrow(res), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
ggsave(glue('/home/ah2428/palmer_scratch/diffpeaks/{celltype}_volcano.png'))

