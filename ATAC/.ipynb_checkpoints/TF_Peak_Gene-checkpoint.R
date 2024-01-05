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

proj <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL')

p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF[c('chr', 'start', 'end')] <- str_split_fixed(p2geneDF$peakName, '_', 3)
p2geneDF <- p2geneDF[,c('idxATAC','peakName','geneName','Correlation','FDR')]

matches <- readRDS('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC_copy/data/ATAC_FINAL/Annotations/Motif-Matches-In-Peaks.rds')

df <- as.data.frame(as.matrix(assay(matches)))
df$idxATAC <- c(1:dim(df)[1])

# Top 20 TFs to create TF-Peak-Gene links
for (celltype in c('EXC','INH','OLI','OPC','AST','END','MIC')) {
    print(celltype)
    celltype_tfs <- read.table(glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/TF_ranks/{celltype}.csv'),sep='\t',header=1)[1:20,]$motif
    celltype_degs_up <- read.table(glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/mast_wilcox_intersected_DEG/{celltype}_PTSD_UP.csv'))$V1
    celltype_degs_down <- read.table(glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/RNA/mast_wilcox_intersected_DEG/{celltype}_PTSD_DOWN.csv'))$V1
    
    tf_peak_gene_dfs <- c()
    for (tf in celltype_tfs) {
        #print(tf)
        tf_peaks <- unique(df[,c(tf,'idxATAC')][df[tf]==TRUE,]$idxATAC)
        tf_name <- str_split(tf,'_')[[1]][1]
        tf_peak_gene_df <- p2geneDF[p2geneDF$idxATAC %in% tf_peaks,]
        tf_peak_gene_df$TF <- tf_name
        tf_peak_gene_df <- tf_peak_gene_df[,c('TF','idxATAC','peakName','geneName','Correlation','FDR')]

        tf_peak_gene_dfs <- rbind(tf_peak_gene_dfs,tf_peak_gene_df)
    }
    print(dim(tf_peak_gene_dfs))
    write.table(tf_peak_gene_dfs,glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/TF_Peak_Gene_2M/{celltype}_TF_Peak_Gene.csv'),quote=F,sep='\t',row.names=F)
}

