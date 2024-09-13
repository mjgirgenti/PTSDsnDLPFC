library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(ggpubr)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(HDF5Array)
library(SummarizedExperiment)

addArchRThreads(threads = 20) 

args = commandArgs(trailingOnly=TRUE)
cluster <- args[1]

proj <- loadArchRProject(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/celltype/{cluster}'))


if (cluster=='EXC') {peak_cluster <- 'Ex'}
if (cluster=='INH') {peak_cluster <- 'In'}
if (cluster=='MIC') {peak_cluster <- 'Mic'}
if (cluster=='END') {peak_cluster <- 'End'}
if (cluster=='OLI') {peak_cluster <- 'Oli'}
if (cluster=='AST') {peak_cluster <- 'Ast'}
if (cluster=='OPC') {peak_cluster <- 'OPC'}

new_peakset <- proj@peakSet[names(proj@peakSet)==peak_cluster,]

proj <- addPeakSet(
  ArchRProj = proj,
  peakSet = new_peakset,
  genomeAnnotation = getGenomeAnnotation(proj),
  force = TRUE
)

proj <- addPeakMatrix(proj)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=T)

saveArchRProject(proj,load=F)

mat <- getMatrixFromProject(proj,useMatrix='PeakMatrix')

print(mat)

markerTest <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_Cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = glue("{cluster}_PTSD"),
  bgdGroups = glue("{cluster}_CON"),
  maxCells=15000
)

pv <- plotMarkers(seMarker = markerTest, name = glue("{cluster}_PTSD"), cutOff = "FDR < 0.05 & abs(Log2FC) >= log2(1.2)", plotAs = "Volcano")
pv
ggsave(glue('/home/ah2428/diffpeaks2/{cluster}_volcano.png'))


