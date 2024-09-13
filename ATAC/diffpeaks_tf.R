library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(ggpubr)
library(Seurat)
library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
cluster <- args[1]

addArchRThreads(threads = 20)

proj <- loadArchRProject(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/celltype/{cluster}'))

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

#df <- data.frame(peak = rownames(markerTest), Log2FC = assays(markerTest)$Log2FC, FDR = assays(markerTest)$FDR)
#colnames(df) <- c('peak','Log2FC','FDR')
#df <- df[order(df$FDR),]

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR < 0.05 & Log2FC > log2(1.2)"
  )

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

write.table(df,glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL/DiffPeaks_TF/{cluster}_Diffpeaks_TF.csv'),sep='\t')

