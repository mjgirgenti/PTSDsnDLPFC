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

proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL')

markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Subclass_Int_Cluster",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

genes <- c('SLC17A7','SV2B','SATB2','GAD1','GAD2','SLC32A1','FLT1','EBF1','CLDN5','AQP4','GJA1','GFAP','MOBP','MOG','MAG','PDGFRA','OLIG1','OLIG2','APBB1IP','CX3CR1','MRC1')


heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.2", 
  labelMarkers = genes,
  transpose = TRUE,
  scaleRows=TRUE,
  returnMatrix=TRUE
)

write.csv(heatmapGS,'/home/ah2428/girgenti_project/ptsd/data/atac_marker_genes_heatmap_scaled_subclass.csv')


