library(ArchR)
library(glue)
library(parallel)

addArchRThreads(threads = 20) 

args = commandArgs(trailingOnly=TRUE)
#cluster <- args[1]
#numcells <- as.integer(args[2])

#proj <- loadArchRProject(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/celltype_recall_peaks/{cluster}'))

proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL')

getGroupBW(
  ArchRProj = proj,
  groupBy = "Int_Cluster",
  normMethod = "ReadsInTSS",
  tileSize = 1000,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
