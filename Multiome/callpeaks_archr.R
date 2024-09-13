library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)

addArchRThreads(threads = 20)

args = commandArgs(trailingOnly=TRUE)

group <- args[1]
 
proj <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples_archr')

proj <- addGroupCoverages(ArchRProj = proj, groupBy = group, force=TRUE)

pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = group, 
    pathToMacs2 = pathToMacs2
)

saveArchRProject(proj)

proj <- addPeakMatrix(proj)


