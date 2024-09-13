library(stringr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(glue)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

name <- args[1]
data <- LoadH5Seurat(glue('/home/ah2428/scratch60/wgcna_subclass_data/{name}.h5seurat'),assays='RNA')
meta <- read.delim2(glue('/home/ah2428/scratch60/wgcna_subclass_data/{name}.csv'))

source("create_metacells.R")

seurat <- FindNeighbors(data, dims=1:20, k.param=50, reduction='pca', return.neighbor=TRUE)
items <- construct_metacell_new(seurat)
barcodes <- colnames(seurat)
write.table(barcodes,glue('/home/ah2428/scratch60/wgcna_subclass_data/metacells/{name}_barcodes.txt'),sep='\t',quote=FALSE)
write.table(items[[2]],glue('/home/ah2428/scratch60/wgcna_subclass_data/metacells/{name}_metacells.txt'),sep='\t',quote=FALSE)