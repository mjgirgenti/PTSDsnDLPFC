library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 20)

proj <- loadArchRProject('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/data/ATAC_FINAL')

args = commandArgs(trailingOnly=TRUE)
motif <- args[1]

motif <- paste0(motif,'_')
motif <- getFeatures(proj, select = paste(motif, collapse="|"), useMatrix = "MotifMatrix")
motif <- grep("z:", motif, value = TRUE)
motif
motif_name <- str_split(motif,':')[[1]][2]

# plot umap

p <- plotEmbedding(
        ArchRProj = proj, 
        colorBy = "MotifMatrix", 
        name = motif, 
        embedding = "UMAP",
        imputeWeights = getImputeWeights(proj),
        pal = ArchRPalettes$whiteBlue
    )

plotPDF(p, name = glue('/TF_umap/{motif_name}_TF_umap'), width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# footprint data

celltype <- args[2]

cluster_proj <- loadArchRProject(glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/celltype/{celltype}'))

motifPositions <- getPositions(cluster_proj,name='Motif',annoName=motif_name)
seFoot <- getFootprints(ArchRProj = proj, positions = motifPositions[motif_name], groupBy = "Cluster_Cond", useGroups=c(glue('{celltype}_PTSD'),glue('{celltype}_CON')))

source("/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/misc/girgenti_project/ArchR/R/plot_footprints.R")
fp_data = plotFootprintsManual(seFoot = seFoot,ArchRProj = cluster_proj,pal=c('#272E6A','#D51F26'),normMethod = "Subtract",addDOC = FALSE,smoothWindow = 5,width=10,plot=F)

write.csv(fp_data,glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/TF_data/{motif_name}_{celltype}_footprint_data.txt'),quote=F)
