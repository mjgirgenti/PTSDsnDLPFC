library(MAST)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(glue)
library(data.table)
library(stringr)
library(edgeR)
library(dplyr)
library(DESeq2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--celltype', type='character', help='celltype')
parser$add_argument('--condition', type='character', help='condition')
args <- parser$parse_args()

celltype <- args$celltype
cond <- args$condition

data <- LoadH5Seurat(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{cond}_vs_CON/{cond}_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
meta <- read.delim2(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{cond}_vs_CON/{cond}_CON_{celltype}_meta.txt'))

protein <- read.table('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt')$V2
protein_genes <- intersect(rownames(data),protein)
data <- data[protein_genes,]

counts <- data$RNA@counts
bin_counts <- as.matrix((counts > 0) + 0)

num_con <- floor(dim(meta[meta$Condition=='CON',])[1]*0.05)
num_cond <- floor(dim(meta[meta$Condition==cond,])[1]*0.05)

con_sums <- rowSums(bin_counts[, which(meta$Condition=='CON')])
cond_sums <- rowSums(bin_counts[, which(meta$Condition==cond)])
df <- data.frame(con=con_sums,cond=cond_sums)
filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
data2 <- data[filtered_genes,]

bin_counts2 <- as.matrix((data2$RNA@counts > 0) + 0)
binary_cell_counts <- colSums(bin_counts2)
write.table(binary_cell_counts,glue('/home/ah2428/ShareZhangLab/PTSD/RNA/ngenes/{celltype}_{cond}_binary_cell_counts_robust.txt'),quote=F)


