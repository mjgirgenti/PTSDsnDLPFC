library(ArchR)
library(stringr)
library(parallel)

addArchRThreads(threads = 10)
addArchRGenome("hg38")

proj <- loadArchRProject('/home/ah2428/girgenti_scratch60/atac/94_samples_tss4_100K_3_removed_3_5_7_13_16_17_18')
print(proj)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",force=TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",sampleCells=100000,force=TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force=TRUE)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj <- addImputeWeights(proj)

proj <- saveArchRProject(proj)

Micro <- c("MRC1","TMEM119","CX3CR1","APBB1IP","TREM1")
Endo <- c("CLDN5","FLT1","DUSP1","COBLL1")
Astros <- c("GLUL","SOX9","AQP4","GJA1","NDRG2","GFAP","ALDH1A1","ALDH1L1","VIM","SLC4A4")
OPC <- c("PDGFRA", "PCDH15", "OLIG2","OLIG1")
Oligos <- c("PLP1","MAG","MOG","MOBP","MBP")
Exc <- c("SATB2","SLC17A7","GRM4")
Inh <- c("GAD2","GAD1",'SLC32A1',"SST","PVALB","VIP")
Neu <- c("STMN2","RBFOX3")
markerGenes <- rev(c(Micro,Endo,Astros,OPC,Oligos,Exc,Inh,Neu))
p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p,
    name = "genes.pdf",
    ArchRProj = proj,
    addDOC = FALSE, width = 5, height = 5)

