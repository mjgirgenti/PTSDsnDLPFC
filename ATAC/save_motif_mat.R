library(Seurat)
library(SeuratDisk)
library(ArchR)
library(glue)
library(parallel)
library(HDF5Array)
library(SummarizedExperiment)

addArchRThreads(threads = 20) 

proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL')

motifmatrix <- getMatrixFromProject(proj,useMatrix='MotifMatrix')

saveHDF5SummarizedExperiment(motifmatrix,'/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL/MotifMatrix')
