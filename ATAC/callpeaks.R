library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)

addArchRThreads(threads = 20)
 
proj <- loadArchRProject('/home/ah2428/girgenti_scratch60/atac/final_proj_copy2')

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Int_Cluster", force=TRUE)

pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Int_Cluster", 
    pathToMacs2 = pathToMacs2,
    cutOff = 0.00001,
    maxPeaks = 10000
)

saveArchRProject(proj)

proj <- addPeakMatrix(proj)
