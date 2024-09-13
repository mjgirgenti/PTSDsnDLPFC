library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(HDF5Array)
library(SummarizedExperiment)
library(Signac)
library(tibble)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 20)

proj <- loadArchRProject('/home/ah2428/palmer_scratch/ATAC_FINAL3_condition')

print(proj)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=T)

proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
ggsave('/home/ah2428/palmer_scratch/motif_variability.png')

saveArchRProject(proj,load=F)

