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

proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/final_proj_coaccess_all_cells')

meta <- getCellColData(proj)

#ptsd_cells <- rownames(meta[meta$Condition=='PTSD',])

proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)

saveArchRProject(proj,load=F)
