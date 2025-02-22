library(MAST)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(EnhancedVolcano)
library(glue)
library(data.table)
library(stringr)
library(edgeR)
library(dplyr)
library(DESeq2)
library(R.utils)
library(argparse)
library(glmGamPoi)
library(compositions)

parser <- ArgumentParser()
parser$add_argument('--condition', type='character')
parser$add_argument('--celltype', type='character')
args <- parser$parse_args()

cond <- args$condition
celltype <- args$celltype

data <- LoadH5Seurat(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{cond}_vs_CON/{cond}_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
meta <- fread(glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{cond}_vs_CON/{cond}_CON_{celltype}_meta.txt'))

protein <- fread('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt',header=FALSE)$V2
protein_genes <- intersect(rownames(data),protein)
data <- data[protein_genes,]
print(dim(data))

counts <- data$RNA@counts
counts <- as(counts, "dgCMatrix")
bin_counts <- (counts > 0) 
num_con <- floor(sum(meta$Condition == 'CON') * 0.2)
num_cond <- floor(sum(meta$Condition == cond) * 0.2)
con_sums <- rowSums(bin_counts[, which(meta$Condition=='CON')])
cond_sums <- rowSums(bin_counts[, which(meta$Condition==cond)])
df <- data.frame(con=con_sums,cond=cond_sums)
filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
data2 <- data[filtered_genes,]

data2@meta.data$Sample <- sapply(strsplit(rownames(data2@meta.data), "-"), head, 1)
samples <- unique(data2@meta.data$Sample)
remove <- samples[table(data2@meta.data$Sample) < 30]
Idents(data2) <- data2@meta.data$Sample
keep = setdiff(samples,remove)
data2 <- subset(data2,idents = keep)
print(dim(data2))

counts <- data.frame(data2$RNA@counts)
strings <- colnames(counts)
new_colnames <- sapply(strsplit(strings,'[.]'),`[`,1)
colnames(counts) <- new_colnames
samples <- unique(colnames(counts))
counts_df <- sapply(split.default(counts, colnames(counts)), rowSums)

meta2 <- meta[meta$Channel %in% samples,]
meta2 <- meta2[,c('Channel','Condition','Sex','PMI','RIN','Race','AgeDeath')]
meta2 <- meta2[!duplicated(meta2$Channel), ]
meta2 <- meta2[,c('Condition','Sex','PMI','RIN','Race','AgeDeath')]

meta2$Condition <- relevel(as.factor(meta2$Condition),'CON')
meta2$Sex<- relevel(as.factor(meta2$Sex),'M')
meta2$Race <- as.factor(meta2$Race)
meta2$AgeDeath <- scale(as.numeric(meta2$AgeDeath))
meta2$RIN <- scale(as.numeric(meta2$RIN))
meta2$PMI <- scale(as.numeric(meta2$PMI))

cdr <- log2(colSums(counts_df))
meta2$cdr <- scale(as.numeric(cdr))

dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = meta2, design = ~ cdr + AgeDeath + PMI + RIN + Sex + Race + Condition)
dds <- DESeq(dds, test = "LRT", useT = TRUE, reduced = ~ cdr + AgeDeath + PMI + RIN + Sex + Race, 
             minmu = 1e-6, minReplicatesForReplace = Inf, fitType = 'glmGamPoi')

res <- results(dds)
write.table(res,glue('/home/ah2428/ShareZhangLab/ah2428/PTSD/deseq2_final/{celltype}_{cond}_DEG.csv'),sep='\t')

deg <- res
sig_deg = deg[ which(deg$padj < 0.1 & abs(deg$log2FoldChange) > log2(1.2)), ]

num_up <- nrow(sig_deg[sig_deg$log2FoldChange > 0,])
num_down <- nrow(sig_deg[sig_deg$log2FoldChange < 0,])

EnhancedVolcano(deg,subtitle=NULL,
                lab=rownames(deg),x='log2FoldChange',y='padj',pCutoff=0.1,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", nrow(deg), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
ggsave(glue('/home/ah2428/ShareZhangLab/ah2428/PTSD/deseq2_final/{celltype}_{cond}_volcano.png'))

