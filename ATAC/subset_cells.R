library(ArchR)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)

addArchRThreads(threads = 10) 

args = commandArgs(trailingOnly=TRUE)
cluster <- args[1]
 
proj <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL2')

subset <- subsetArchRProject(proj,cells=getCellNames(proj)[proj$Int_Cluster==cluster],outputDirectory=glue('/home/ah2428/palmer_scratch/celltype/{cluster}'),force=TRUE)
subset <- subsetCells(proj,cellNames=getCellNames(proj)[proj$Int_Cluster==cluster])
print(subset)
#subset <- addGroupCoverages(ArchRProj = subset, groupBy = "Clusters_Cond", force=TRUE)
#subset <- addMotifAnnotations(ArchRProj = subset, motifSet = "cisbp", name = "Motif", force=TRUE)

#saveArchRProject(subset,load=FALSE)




