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
library(rtracklayer)

addArchRThreads(threads = 20) 

args = commandArgs(trailingOnly=TRUE)
celltype <- args[1]

proj <- loadArchRProject(glue('/home/ah2428/girgenti/ah2428/PTSD/data/celltype_condition/{celltype}'))

peaks <- read.table('/home/ah2428/palmer_scratch/celltype_condition_union_peaks.csv',sep='\t')
diffpeaks <- read.table(glue('/home/ah2428/palmer_scratch/diffpeaks/{celltype}_diffpeaks.csv'),sep='\t')
diffpeaks <- diffpeaks[ which(diffpeaks$pvalue < 0.1 & abs(diffpeaks$log2FoldChange) > log2(1.2)), ]
up_peaks <- rownames(diffpeaks[diffpeaks$log2FoldChange > 0,])
down_peaks <- rownames(diffpeaks[diffpeaks$log2FoldChange < 0,])
celltype_ptsd_peaks <- peaks[peaks$celltype==glue('{celltype}_PTSD'),]
celltype_con_peaks <- peaks[peaks$celltype==glue('{celltype}_CON'),]
celltype_ptsd_peaks <- paste0(celltype_ptsd_peaks$seqnames,'-',celltype_ptsd_peaks$start,'-',celltype_ptsd_peaks$end)
celltype_con_peaks <- paste0(celltype_con_peaks$seqnames,'-',celltype_con_peaks$start,'-',celltype_con_peaks$end)
up_diffpeaks <- intersect(up_peaks,celltype_ptsd_peaks)
down_diffpeaks <- intersect(down_peaks,celltype_con_peaks)

up <- diffpeaks[rownames(diffpeaks) %in% up_diffpeaks,]
up[c('seqnames','start','end')] <- str_split_fixed(rownames(up), '-', 3)
up <- up[,c('seqnames','start','end')]
up_gr <- GRanges(up[,c('seqnames','start','end')])
up_gr <- setNames(up_gr, rep(glue('{celltype}_PTSD'),dim(up)[1]))

down <- diffpeaks[rownames(diffpeaks) %in% down_diffpeaks,]
down[c('seqnames','start','end')] <- str_split_fixed(rownames(down), '-', 3)
down <- down[,c('seqnames','start','end')]
down_gr <- GRanges(down[,c('seqnames','start','end')])
down_gr <- setNames(down_gr, rep(glue('{celltype}_CON'),dim(down)[1]))

gr <- GRanges(c(up_gr,down_gr))
names(gr) <- gr$name

proj <- proj[proj$Cluster_Cond %in% c(glue('{celltype}_PTSD'),glue('{celltype}_CON')),]

proj <- addPeakSet(
  ArchRProj = proj,
  peakSet = gr,
  genomeAnnotation = getGenomeAnnotation(proj),
  force = TRUE
)

proj <- addPeakMatrix(proj)
print('here')
saveArchRProject(proj,load=F)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Cluster_Cond",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells=10000
)
pv <- plotMarkers(seMarker = markersPeaks, name = glue("{celltype}_PTSD"), cutOff = "FDR <= 1 & abs(Log2FC) > log2(1.2)", plotAs = "Volcano")
ggsave(glue('/home/ah2428/palmer_scratch/celltype_condition/{celltype}_archr_volcano.png'))

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=T)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 1 & Log2FC > log2(1.2)"
  )

saveHDF5SummarizedExperiment(enrichMotifs,glue('/home/ah2428/palmer_scratch/celltype_condition/{celltype}_enrichMotifs'))

df <- data.frame(TF = rownames(enrichMotifs), mlog10Padj = assay(enrichMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
write.table(df,glue('/home/ah2428/palmer_scratch/celltype_condition/{celltype}_TF_rank.csv',sep='\t'))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(glue('/home/ah2428/palmer_scratch/celltype_condition/{celltype}_TF_rank_plot.png'))

saveArchRProject(proj,load=F)






