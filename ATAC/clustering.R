library(ArchR)
library(stringr)
library(parallel)

addArchRThreads(threads = 30)
addArchRGenome("hg38")

proj <- loadArchRProject('/home/ah2428/girgenti_scratch60/atac/94_samples')

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
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

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 5),p2))
ggsave('/home/ah2428/girgenti_scratch60/atac/94_samples/Plots/genes.pdf')


