library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(harmony)
library(future)

set.seed(1234)

root_dir <- '/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/multiome_analysis_outputs/'
save_dir <- '/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples_multiome_clustered/'
ids <- read.table('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples.txt')$V1

seurat_objs <- list()
for (x in ids) {
    tmp <- readRDS(paste0(root_dir, x, '/seurat_obj_filtered.rds'))
    seurat_objs[[x]] <- tmp
}

# merge samples
# note: it will take very long for this step to finish

merged <- seurat_objs[[ids[1]]]
i <- 1
for (x in ids[2:length(ids)]) {
    print(paste0('Reading ', x))
    merged <- merge(merged, seurat_objs[[x]])
}

saveRDS(object = merged, file = paste0(save_dir,'multiome_merged.rds'))

# run PCA for merged

options(future.globals.maxSize = 8000 * 1024^2)
options(warn=-1)
DefaultAssay(merged) <- "RNA"
merged <- SCTransform(merged)
merged <- RunPCA(merged)

# run LSI for merged

DefaultAssay(merged) <- "ATAC"

merged <- RunTFIDF(merged)
merged <- FindTopFeatures(merged, min.cutoff = 'q0')
merged <- RunSVD(merged)

# run joint embedding for merged

# build a joint neighbor graph using both assays
merged <- FindMultiModalNeighbors(
  object = merged,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 1:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
merged <- RunUMAP(
  object = merged,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# plot merged colored by samples

options(repr.plot.height = 8, repr.plot.width = 12)
DimPlot(object = merged, label = TRUE, group.by = 'dataset')
ggsave(paste0(save_dir,'samples.pdf'))

# run harmony

DefaultAssay(merged) <- "SCT"
options(repr.plot.height = 2.5, repr.plot.width = 6)
merged <- merged %>% 
    RunHarmony("dataset", plot_convergence = TRUE, assay.use = 'SCT')

# plot embedding after harmony

merged <- RunUMAP(merged, reduction = "harmony", dims = 1:50, assay = 'SCT')
options(repr.plot.height = 8, repr.plot.width = 12)
DimPlot(object = merged, label = TRUE, group.by = 'dataset')
ggsave(paste0(save_dir,'samples_harmony.pdf'))

saveRDS(object = merged, file = paste0(save_dir, 'multiome_merged_harmony.rds'))


