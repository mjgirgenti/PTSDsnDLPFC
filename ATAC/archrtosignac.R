library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(HDF5Array)
library(SummarizedExperiment)
library(Signac)
library(ArchRtoSignac)
library(EnsDb.Hsapiens.v86)

addArchRThreads(threads = 20) 

args = commandArgs(trailingOnly=TRUE)
celltype <- args[1]

#proj <- loadArchRProject(glue('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/celltype_condition/{celltype}'))

proj <- loadArchRProject('/home/ah2428/palmer_scratch/ATAC_FINAL3_condition')
proj <- proj[proj$Cluster_Cond %in% c(glue('{celltype}_PTSD'),glue('{celltype}_CON')),]
mat <- getPeakMatrix(proj)

annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38")

fragments_dir <- "/gpfs/slayman/pi/gerstein/jz435/Share-ZhangLab/scATAC/frag/" # please see the fragments format provided by ArchR examples

# Conversion function
seurat_atac <- ArchR2Signac(
  ArchRProject = proj,
  samples = unique(proj$Sample), # Provide a list of unique sample
  fragments_dir = fragments_dir, # the folder that contains all fragments samples in '.fragments.tsv.gz' or '.tsv.gz'
  pm = mat, # geting peak martix
  fragments_fromcellranger = "NO",
  fragments_file_extension = '.fragments.tsv.gz',
  refversion = 'hg38', # write the EnsDb version
  annotation = annotations
)

print(seurat_atac)
SaveH5Seurat(seurat_atac, glue('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/celltype_condition/{celltype}_signac.h5seurat'), overwrite = TRUE)



