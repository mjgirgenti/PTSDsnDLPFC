library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(ggpubr)
library(Seurat)
library(SeuratDisk)
library(HDF5Array)
library(SummarizedExperiment)
library(Signac)

proj = loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL')

cells <- getCellNames(proj)

subsetArchRProject(proj,cells[1:200000],outputDirectory='/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/proj1')

subsetArchRProject(proj,cells[200001:473033],outputDirectory='/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/proj2')


