library(edgeR)
library(MAST)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(EnhancedVolcano)
library(glue)
library(data.table)
library(stringr)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--celltype', type='character', help='celltype')
args <- parser$parse_args()

celltype <- args$celltype
cond <- 'PTSD'

data <- LoadH5Seurat(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/PTSD_vs_CON/PTSD_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
meta <- read.delim2(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/PTSD_vs_CON/PTSD_CON_{celltype}_meta.txt'))
meta2 <- meta[,c('Channel','Condition')]
meta2 <- meta2[!duplicated(meta2$Channel), ]

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

counts <- data.frame(data2$RNA@counts)

strings <- colnames(counts)
new_colnames <- sapply(strsplit(strings,'[.]'),`[`,1)

colnames(counts) <- new_colnames
samples <- unique(colnames(counts))

counts_df <- sapply(split.default(counts, colnames(counts)), rowSums)

dge <- DGEList(counts = counts_df, group = factor(meta2$Condition))

dge <- calcNormFactors(object = dge)

dge <- estimateDisp(y = dge)

et <- exactTest(object = dge)

top_degs = topTags(object = et, n = "Inf")

deg <- top_degs$table
write.table(deg,glue('/home/ah2428/palmer_scratch/edger/{celltype}_edger_DEG.csv'),sep='\t')

sig_deg = deg[ which(deg$FDR < 0.1 & abs(deg$logFC) > log2(1.2)), ]

num_up <- nrow(sig_deg[sig_deg$logFC > 0,])
num_down <- nrow(sig_deg[sig_deg$logFC < 0,])

EnhancedVolcano(deg,subtitle=NULL,
                lab=rownames(deg),x='logFC',y='FDR',pCutoff=0.1,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", nrow(deg), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))

ggsave(glue('/home/ah2428/palmer_scratch/edger/{celltype}_edger_volcano.png'))

