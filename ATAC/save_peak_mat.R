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
library(SingleCellExperiment)
library(zellkonverter)
library(anndata)

addArchRThreads(threads = 20) 

proj1 <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/proj1')

mat <- getMatrixFromProject(proj1,useMatrix='PeakMatrix')
print('got cell x peak matrix')
print(mat)

sce <- SingleCellExperiment(assays = list(counts = mat))
print(sce)
adata <- SCE2AnnData(sce)

writeH5AD(adata,'/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/proj1_peak_matrix.h5ad')
#atac_ad <- Convert(from = sce, to = "anndata", filename = "/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/proj1_peak_matrix.h5ad")
