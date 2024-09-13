library(ArchR)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)

addArchRThreads(threads = 10) 

proj <- loadArchRProject('/home/ah2428/girgenti_scratch60/atac/94_samples_tss4_no_harmony_final_clusters3')

#proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
motifPositions <- getPositions(proj)

#motifs <- c("JUNB", "NHLH2", "MSC", "NFIC", "KLF5", "SPIB", "MECP2", "SOX9" )
motifs <- 'NFIC'
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
print(markerMotifs)

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters_Cond")

saveArchRProject(proj,load=FALSE)

seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters_Cond",
  useGroups = c('Ast_PTSD','Ast_CON')
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

saveArchRProject(proj)



