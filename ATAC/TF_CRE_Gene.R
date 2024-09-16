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

proj <- loadArchRProject('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL')

matches <- readRDS('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL/Annotations/Motif-Matches-In-Peaks.rds')
df <- as.data.frame(as.matrix(assay(matches)))
df$idxATAC <- c(1:dim(matches)[1])

p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF[c('chr', 'start', 'end')] <- str_split_fixed(p2geneDF$peakName, '_', 3)
p2geneDF <- p2geneDF[,c('idxATAC','peakName','geneName','Correlation','FDR')]
p2geneDF <- p2geneDF[p2geneDF$Correlation > 0.4 & p2geneDF$FDR < 0.05,]

deg <- read.table('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_results_df/snRNA/PTSD_vs_CON/MAST_WILCOX_intersect.csv',sep='\t',header=1)


# Top 20 TF-Peak-Gene
list1 <- c('EXC','INH','OLI','OPC','AST','END','MIC')
list2 <- c('Ex','In','Oli','OPC','End','Ast','Mic')
for (i in seq_along(list1)) {
    celltype <- list1[i]
    celltype2 <- list2[i]
    print(celltype)
    celltype_tfs <- read.table(glue('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/TF_ranks/{celltype}.csv'),sep='\t',header=1)[1:20,]$motif
    celltype_deg <- deg[deg$Celltype==celltype,]
    celltype_degs_up <- celltype_deg[celltype_deg$MAST.log2FC > 0,]$Genename
    celltype_degs_down <- celltype_deg[celltype_deg$MAST.log2FC < 0,]$Genename
    
    tf_peak_gene_dfs <- c()
    for (tf in celltype_tfs) {
        df2 <- df[,c(tf,'idxATAC')][df[tf]==TRUE,]
        df2$celltype <- rownames(df2)
        celltype_peaks <- grepl(celltype2, df2$celltype)
        matching_strings <- df2$celltype[celltype_peaks]

        tf_peaks <- unique(df2[df2$celltype %in% matching_strings,]$idxATAC)
        tf_name <- str_split(tf,'_')[[1]][1]
        tf_peak_gene_df <- p2geneDF[p2geneDF$idxATAC %in% tf_peaks,]
        tf_peak_gene_df$TF <- tf_name
        tf_peak_gene_df <- tf_peak_gene_df[,c('TF','idxATAC','peakName','geneName','Correlation','FDR')]

        tf_peak_gene_dfs <- rbind(tf_peak_gene_dfs,tf_peak_gene_df)
    }
    print(dim(tf_peak_gene_dfs))
    write.table(tf_peak_gene_dfs,glue('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/TF_CRE_GENE/{celltype}_TF_CRE_GENE.csv'),quote=F,sep='\t',row.names=F)
}

