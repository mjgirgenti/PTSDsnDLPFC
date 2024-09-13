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

addArchRThreads(threads = 20) 

#args = commandArgs(trailingOnly=TRUE)
#cluster <- args[1]

#proj <- loadArchRProject(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/celltype/{cluster}'))
proj <- loadArchRProject('/home/ah2428/palmer_scratch/celltype_condition/EXC_old')

new_peakset <- proj@peakSet[names(proj@peakSet) %in% c('EXC_PTSD','EXC_CON'),]

proj <- addPeakSet(
  ArchRProj = proj,
  peakSet = new_peakset,
  genomeAnnotation = getGenomeAnnotation(proj),
  force = TRUE
)

proj <- addPeakMatrix(proj)


saveArchRProject(proj,load=F)



