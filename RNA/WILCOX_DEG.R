library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(EnhancedVolcano)
library(glue)
library(data.table)
library(stringr)
library(edgeR)
library(dplyr)
library(argparse)
library(presto)

parser <- ArgumentParser()
parser$add_argument('--celltype', type='character', help='celltype')
parser$add_argument('--condition', type='character', help='condition')
parser$add_argument('--sex', type='character', help='sex', default='F')
parser$add_argument('--datapath', type='character', help='datapath')
parser$add_argument('--savepath', type='character', help='savepath')
parser$add_argument('--cdr', type='character', help='cdr', default='F')
args <- parser$parse_args()

celltype <- args$celltype
cond <- args$condition
sex <- args$sex
datapath <- args$datapath
savepath <- args$savepath
cdr <- args$cdr
print(celltype)
print(cond)

if (sex=='T') {
    data <- LoadH5Seurat(glue('{datapath}/CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
    meta <- read.delim2(glue('{datapath}/CON_{celltype}_meta.txt'))
    protein <- read.table('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt')$V2
    protein_genes <- intersect(rownames(data),protein)
    data <- data[protein_genes,]
    print(dim(data))
    data <- NormalizeData(data, normalization.method = "RC", scale.factor = 1e6)
    num_con <- floor(dim(meta[meta$Sex=='F',])[1]*0.05)
    num_cond <- floor(dim(meta[meta$Sex==cond,])[1]*0.05)
    counts <- data$RNA@counts
    bin_counts <- as.matrix((counts > 0) + 0)
    con_sums <- rowSums(bin_counts[, which(meta$Sex=='F')])
    cond_sums <- rowSums(bin_counts[, which(meta$Sex==cond)])
    df <- data.frame(con=con_sums,cond=cond_sums)
    filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
    data <- data[filtered_genes,]
    print(dim(data))
    data$sex <- meta$Sex
    Idents(data) <- 'sex'
    deg <- FindMarkers(data, ident.1=glue('{cond}'), ident.2='F', slot='data', test.use='wilcox', logfc.threshold=0, min.pct=0)
    title <- glue('{celltype} CON M vs CON F WILCOX')
} else {
    data <- LoadH5Seurat(glue('{datapath}/{cond}_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
    meta <- read.delim2(glue('{datapath}/{cond}_CON_{celltype}_meta.txt'))
    protein <- read.table('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt')$V2
    protein_genes <- intersect(rownames(data),protein)
    data <- data[protein_genes,]
    print(dim(data))
    data <- NormalizeData(data, normalization.method = "RC", scale.factor = 1e6)
    num_con <- floor(dim(meta[meta$Condition=='CON',])[1]*0.05)
    num_cond <- floor(dim(meta[meta$Condition==cond,])[1]*0.05)
    counts <- data$RNA@counts
    bin_counts <- as.matrix((counts > 0) + 0)
    con_sums <- rowSums(bin_counts[, which(meta$Condition=='CON')])
    cond_sums <- rowSums(bin_counts[, which(meta$Condition==cond)])
    df <- data.frame(con=con_sums,cond=cond_sums)
    filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
    data <- data[filtered_genes,]
    print(dim(data))
    data$condition <- meta$Condition
    Idents(data) <- 'condition'
    data$cdr <- scale(data$nCount_RNA)
    if (cdr=='T') {
        data$RNA@data <- data$RNA@data/data$cdr
        deg <- FindMarkers(data, ident.1=glue('{cond}'), ident.2='CON', slot='data', test.use='wilcox', logfc.threshold=0, min.pct=0)
    } else {
        deg <- FindMarkers(data, ident.1=glue('{cond}'), ident.2='CON', slot='data', test.use='wilcox', logfc.threshold=0, min.pct=0)
    }
    title <- glue('{celltype} {cond} vs CON WILCOX')
}

sig_deg = deg[ which(deg$p_val_adj < 0.01 & abs(deg$avg_log2FC) > log2(1.2)), ]
write.table(deg,glue('{savepath}/{celltype}_WILCOX_DEG.csv'),sep='\t')

num_up <- nrow(sig_deg[sig_deg$avg_log2FC > 0,])
num_down <- nrow(sig_deg[sig_deg$avg_log2FC < 0,])

EnhancedVolcano(deg,subtitle=NULL,title=title,
                lab=rownames(deg),x='avg_log2FC',y='p_val_adj',pCutoff=0.01,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", length(filtered_genes), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))

ggsave(glue('{savepath}/{celltype}_WILCOX_volcano.png'))
