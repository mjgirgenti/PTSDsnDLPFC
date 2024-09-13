library(ArchR)
library(stringr)
library(parallel)

set.seed(1)

addArchRThreads(threads = 20) 
addArchRGenome("hg38")

path <- '/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples_frag'
inputFiles <- Sys.glob(paste0(path,'/*atac_fragments.tsv.gz'))

sampleNames <- sapply(strsplit(inputFiles,'/'),`[`,length(strsplit(path,'/')[[1]])+1)
sampleNames <- sapply(strsplit(sampleNames,'[.]'),`[`,1)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 2, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force=T,
  threads=1
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples_tss2",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
#proj <- filterDoublets(ArchRProj = proj)

proj <- saveArchRProject(ArchRProj = proj)
