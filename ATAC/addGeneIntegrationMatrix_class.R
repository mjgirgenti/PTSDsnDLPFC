library(Seurat)
library(SeuratDisk)
library(ArchR)
library(glue)
library(parallel)

addArchRThreads(threads = 10)

atac_proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL')
#atac_proj <- addImputeWeights(atac_proj)
#saveArchRProject(atac_proj)

rna_proj <- LoadH5Seurat('/home/ah2428/girgenti_project/pegasus/data/RNA_6_samples_for_archr_integration.h5seurat')
rna_meta <- read.csv('/home/ah2428/girgenti_project/pegasus/data/RNA_6_samples_for_archr_integration_meta.csv',sep=',')
rna_proj@meta.data['celltype'] = rna_meta$class

proj <- addGeneIntegrationMatrix(
    ArchRProj = atac_proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = rna_proj,
    addToArrow = FALSE,
    groupRNA = "celltype",
    nameCell = "Int_Cell",
    nameGroup = "Int_Cluster",
    nameScore = "Int_Score",
    force=TRUE
)
getAvailableMatrices(proj)

saveArchRProject(proj)

pal <- paletteDiscrete(values = unique(rna_proj@meta.data$celltype))

p1 <- plotEmbedding(
    proj,
    colorBy = "cellColData",
    name = "Int_Cluster",
    pal = pal
)

plotPDF(p1, name = "Plot-UMAP-RNA-Integration-FINAL.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)





