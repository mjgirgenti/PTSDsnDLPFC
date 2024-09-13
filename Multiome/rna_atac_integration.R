library(Seurat)
library(SeuratDisk)
library(ArchR)
library(glue)
library(parallel)

addArchRThreads(threads = 1)

atac_proj <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL')

rna_proj <- LoadH5Seurat('/home/ah2428/girgenti_project/pegasus/data/RNA_6_samples_for_archr_integration.h5seurat')
rna_meta <- read.csv('/home/ah2428/girgenti_project/pegasus/data/RNA_6_samples_for_archr_integration_meta.csv',sep=',')
rna_proj@meta.data['celltype'] = rna_meta$class

cM <- as.matrix(confusionMatrix(atac_proj$Clusters, atac_proj$Clusters3))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) 

clustAst <- rownames(cM)[grep('Ast', preClust)]
clustEnd <- rownames(cM)[grep('End', preClust)]
clustEx <- rownames(cM)[grep('Ex', preClust)]
clustIn <- rownames(cM)[grep('In', preClust)]
clustMic <- rownames(cM)[grep('Mic', preClust)]
clustOli <- rownames(cM)[grep('Oli', preClust)]
clustOPC <- rownames(cM)[grep('OPC', preClust)]

rnaAst <- colnames(rna_proj)[rna_proj@meta.data$celltype %in% 'AST']
rnaEnd <- colnames(rna_proj)[rna_proj@meta.data$celltype %in% 'END']
rnaEx <- colnames(rna_proj)[rna_proj@meta.data$celltype %in% 'EXC']
rnaIn <- colnames(rna_proj)[rna_proj@meta.data$celltype %in% 'INH']
rnaMic <- colnames(rna_proj)[rna_proj@meta.data$celltype %in% 'MIC']
rnaOli <- colnames(rna_proj)[rna_proj@meta.data$celltype %in% 'OLI']
rnaOPC <- colnames(rna_proj)[rna_proj@meta.data$celltype %in% 'OPC']

groupList <- SimpleList(
    Ast = SimpleList(ATAC = atac_proj$cellNames[atac_proj$Clusters %in% clustAst],RNA = rnaAst),
    End = SimpleList(ATAC = atac_proj$cellNames[atac_proj$Clusters %in% clustEnd],RNA = rnaEnd),
    Ex = SimpleList(ATAC = atac_proj$cellNames[atac_proj$Clusters %in% clustEx],RNA = rnaEx),
    In = SimpleList(ATAC = atac_proj$cellNames[atac_proj$Clusters %in% clustIn],RNA = rnaIn),
    Mic = SimpleList(ATAC = atac_proj$cellNames[atac_proj$Clusters %in% clustMic],RNA = rnaMic),
    Oli = SimpleList(ATAC = atac_proj$cellNames[atac_proj$Clusters %in% clustOli],RNA = rnaOli),
    OPC = SimpleList(ATAC = atac_proj$cellNames[atac_proj$Clusters %in% clustOPC],RNA = rnaOPC)
)   

proj <- addGeneIntegrationMatrix(
    ArchRProj = atac_proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = rna_proj,
    addToArrow = TRUE, 
    groupList = groupList,
    groupRNA = "celltype",
    nameCell = "Int_Cell",
    nameGroup = "Int_Cluster",
    nameScore = "Int_Score",
    threads=1,
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

plotPDF(p1, name = "Plot-UMAP-RNA-Constrained-Integration.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)





