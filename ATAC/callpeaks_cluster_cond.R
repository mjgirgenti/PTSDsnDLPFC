library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)


addArchRThreads(threads = 20)

args = commandArgs(trailingOnly=TRUE)
 
proj <- loadArchRProject('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL2')

print(unique(proj$Cluster_Cond))

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Cluster_Cond", force=TRUE)

pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Cluster_Cond", 
    pathToMacs2 = pathToMacs2,
    force=T,
    threads=1
)

saveArchRProject(proj)

proj <- addPeakMatrix(proj)

saveArchRProject(proj)


