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

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(proj,load=F)
