library(Seurat)
library(SeuratDisk)
library(ArchR)
library(glue)
library(parallel)

addArchRThreads(threads = 20)

proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL')

proj<- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)

saveArchRProject(proj,load=FALSE)
