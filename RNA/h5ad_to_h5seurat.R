library(Seurat)
library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)

Convert(args[1], dest = args[2], overwrite = TRUE)


