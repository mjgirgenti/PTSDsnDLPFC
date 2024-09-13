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

ast <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/AST')

ast <- addMotifAnnotations(ArchRProj = ast, motifSet = "cisbp", name = "Motif", force=T)

markerTest <- getMarkerFeatures(
  ArchRProj = ast, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters_Cond",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AST_PTSD",
  bgdGroups = "AST_CON",
  maxCells=10000
)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = ast,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Ast-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = ast, addDOC = FALSE)

saveArchRProject(ast,load=F)


