library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(ggpubr)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(HDF5Array)
library(SummarizedExperiment)
library(Signac)

addArchRThreads(threads = 20) 

args = commandArgs(trailingOnly=TRUE)
cluster <- args[1]

enrichMotifs <- loadHDF5SummarizedExperiment('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL/EnrichMotifs')
dev <- loadHDF5SummarizedExperiment('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL/MotifMatrix')

motifs <- data.frame(mlog10Padj=assays(enrichMotifs[,glue({cluster})])$mlog10Padj[,1])
motifs$TF <- rownames(motifs)
motifs <- motifs[order(-motifs$mlog10Padj),]
num_motifs <- dim(motifs)[1]
motifs$rank <- c(1:num_motifs)

cluster_dev <- dev[, dev$Cluster_Cond %in% c(glue('{cluster}_PTSD'),glue('{cluster}_CON'))]

pvals <- c()
diffs <- c()
ranks <- c()
enrichments <- c()

for (i in rownames(motifs[1:num_motifs,])) {
    rank <- motifs[i,]$rank
    enrichment <- motifs[i,]$mlog10Padj

    z_score <- assays(cluster_dev)$z[i,]
    df <- data.frame(z_score)
    df$Cluster_Cond <- colData(cluster_dev)$Cluster_Cond
    df['Cluster_Cond'][df['Cluster_Cond'] == glue('{cluster}_PTSD')] <- 'PTSD'
    df['Cluster_Cond'][df['Cluster_Cond'] == glue('{cluster}_CON')] <- 'CON'
    res <- compare_means(z_score ~ Cluster_Cond, data=df, method='wilcox.test', ref.group='PTSD')
    pval <- res$p
    pvals <- append(pvals, pval)
    diff <- mean(df[df$Cluster_Cond=='PTSD',]$z_score) - mean(df[df$Cluster_Cond=='CON',]$z_score)
    diffs <- append(diffs,diff)
    ranks <- append(ranks,rank)
    enrichments <- append(enrichments,enrichment)
}
pvals_adj <- p.adjust(pvals)

out <- data.frame(motif=rownames(motifs[1:num_motifs,]),pval=pvals,padj=pvals_adj,rank=ranks,enrichment=enrichments,meandiff=diffs)
out <- out[order(out$padj),]

write.table(out,glue('/home/ah2428/tf_ranks/{cluster}.csv'),sep='\t',quote=F,row.names=F)
