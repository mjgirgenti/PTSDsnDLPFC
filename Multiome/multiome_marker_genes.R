library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(harmony)
library(future)

set.seed(1234)
save_dir <- '/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples_multiome_clustered/'
data <- readRDS('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples_multiome_clustered/multiome_annotated_subset.rds')

data <- RunUMAP(data, reduction = "harmony", dims = 1:50, assay = 'SCT')

DimPlot(object = data, label = TRUE)
ggsave(paste0(save_dir,'multiome_umap_clusters.pdf'))

saveRDS(object=data, file=paste0(save_dir,'multiome_clustering2.rds'))

genes <- c('SATB2','CUX2','RORB','FEZF2','OPRK1','GAD1','SST','VIP','PVALB','LAMP5','KCNG1','PDGFRA','PCDH15','PLP1','MAG','MRC1','APBB1IP','CLDN5','FLT1','AQP4','GJA1')
DotPlot(object = data, features = genes) + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(save_dir,'dotplot.pdf'))

FeaturePlot(data, genes, reduction = 'umap')
ggsave(paste0(save_dir,'genes.pdf'))

