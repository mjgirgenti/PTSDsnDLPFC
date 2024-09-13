library(ArchR)
library(glue)
library(parallel)

addArchRThreads(threads = 20) 


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
