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
library(rliger)

proj1 <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/proj1') 
proj2 <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/proj2')

peakmat1 <- getMatrixFromProject(proj1,useMatrix='PeakMatrix')
peakmat2 <- getMatrixFromProject(proj2,useMatrix='PeakMatrix')

mat_all <- cbind(peakmat1,peakmat2)


rownames(mat_all) <- rowData(mat_all)$name
atac_peaks <- rowData(mat_all)
mat <- as.matrix(assays(mat_all)$PeakMatrix)
atac <- CreateSeuratObject(counts=mat,assay='ATAC')

SaveH5Seurat(atac, '/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/peakmat.h5seurat')
