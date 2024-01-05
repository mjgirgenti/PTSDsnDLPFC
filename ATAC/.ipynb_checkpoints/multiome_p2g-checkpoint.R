library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)

# Extended Data Figure 7a: Multiome P2G Heatmap

proj <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/MULTIOME_P2G')

genes <- getGenes(proj)$symbol

rna <- LoadH5Seurat('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/data/multiome_RNA.h5seurat')
rna_mat <- rna[['RNA']]@data
inter_genes <- intersect(rownames(rna_mat),genes)
rna_mat_subset <- rna_mat[inter_genes,]
rna_exp <- SummarizedExperiment(assays = rna_mat_subset, colData=meta)

meta <- getCellColData(proj)
rna_exp@metadata$subclass = meta$subclass

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    seRNA = rna_exp,
    groupATAC = 'subclass',
    groupRNA = 'subclass',
    force = TRUE
)

proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    maxDist = 2000000,
    reducedDims = "IterativeLSI")

p <- plotPeak2GeneHeatmap(
    ArchRProj = proj, 
    groupBy = "subclass", k=20, 
    palGroup=ArchRPalettes$subtype, 
    limitsATAC = c(-1.5, 1.5), 
    limitsRNA = c(-1.5, 1.5))

plotPDF(p, 
    name = "multiome_p2g_heatmap_14.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 7, height = 10)