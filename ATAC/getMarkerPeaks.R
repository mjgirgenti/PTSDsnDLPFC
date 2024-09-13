library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(HDF5Array)
library(SummarizedExperiment)

addArchRThreads(threads = 10)

proj <- loadArchRProject('/home/ah2428/palmer_scratch/ATAC_FINAL3_condition')

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Cluster_Cond",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveHDF5SummarizedExperiment(markersPeaks,'/home/ah2428/palmer_scratch/MarkersPeaks_Cluster_Cond')

saveArchRProject(proj)
