library(ArchR)
library(stringr)
library(ggplot2)
library(parallel)

addArchRThreads(threads = 20) 

proj <- loadArchRProject('/home/ah2428/girgenti_scratch60/atac/94_samples_tss4_no_harmony_final_clusters3')

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters3",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

proj <- addArchRAnnotations(ArchRProj = proj, collection = "EncodeTFBS")

enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

proj <- addArchRAnnotations(ArchRProj = proj, collection = "ATAC")
enrichATAC <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "ATAC",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

proj <- addArchRAnnotations(ArchRProj = proj, collection = "Codex")
enrichCodex <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Codex",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapCodex, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapCodex, name = "Codex-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

