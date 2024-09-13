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

addArchRThreads(threads = 10)

proj <- loadArchRProject('/home/ah2428/project/PTSD_ATAC/final_proj_copy2')

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Int_Cluster",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveHDF5SummarizedExperiment(markersPeaks,'/home/ah2428/project/PTSD_ATAC/final_proj_copy2/MarkersPeaks')

saveArchRProject(proj)
