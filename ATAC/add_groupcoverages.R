library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(ggpubr)
library(Seurat)
library(SeuratDisk)

addArchRThreads(threads = 20) 

proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/celltype_recall_peaks/EXC')


proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Cluster_Cond", minCells=10000, maxCells=20000, force=TRUE)

saveArchRProject(proj,load=F)


