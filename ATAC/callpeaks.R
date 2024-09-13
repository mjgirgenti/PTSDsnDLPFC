library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)

addArchRThreads(threads = 20)

 
proj <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL')

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Int_Cluster", force=TRUE)

pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Int_Cluster", 
    pathToMacs2 = pathToMacs2
)

saveArchRProject(proj)

proj <- addPeakMatrix(proj)
