library(ArchR)
library(stringr)
library(parallel)


set.seed(1)

addArchRThreads(threads = 20)
addArchRGenome("hg38")

path <- '/ysm-gpfs/pi/gerstein/jz435/Share-ZhangLab/scATAC/frag'
inputFiles <- Sys.glob(paste0(path,'/*.fragments.tsv.gz'))
remove <- str_detect(inputFiles,"MS0031LL|MS0071AA|MS0073CC|MS0077GG|MS0079II|MS0083MM|MS0085OO|MS0151EE")
inputFiles <- inputFiles[!remove]
sampleNames <- sapply(strsplit(inputFiles,'/'),`[`,9)
sampleNames <- sapply(strsplit(sampleNames,'[.]'),`[`,1)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 2, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force=T
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/home/ah2428/girgenti_scratch60/atac/96_samples",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
proj <- filterDoublets(ArchRProj = proj)

proj <- saveArchRProject(ArchRProj = proj)

