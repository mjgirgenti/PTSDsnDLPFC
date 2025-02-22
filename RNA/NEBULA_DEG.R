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
library(R.utils)
library(nebula)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--cond', type='character')
parser$add_argument('--celltype', type='character')
parser$add_argument('--ratio', type='character')
parser$add_argument('--model', type='character')
args <- parser$parse_args()

cond <- args$cond
celltype <- args$celltype
ratio <- args$ratio
model <- 'NBGMM'
kappa <- 800
cell <- 20

savepath <- glue('/home/ah2428/ShareZhangLab/ah2428/PTSD/nebula_final/')
if (!dir.exists(savepath)) {
  dir.create(savepath, recursive=TRUE)
  cat("Directory created:", savepath, "\n")
} else {
  cat("Directory already exists:", savepath, "\n")
}

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

data2@meta.data$Sample <- meta$Channel
sample_counts <- table(data2$Sample)
valid_samples <- names(sample_counts[sample_counts >= 30])
filtered_data <- subset(data2, subset = Sample %in% valid_samples)
data2 <- filtered_data

meta2 = meta[meta$barcodekey %in% colnames(data2),]
meta2 <- meta2[order(data2$Sample),]
data2@meta.data$Sample <- meta2$Channel
data2 <- data2[,meta2$barcodekey]
meta2 <- meta2[,c('AgeDeath','Sex','Race','PMI','RIN','Condition')]

meta2$AgeDeath <- as.numeric(meta2$AgeDeath)
meta2$RIN <- as.numeric(meta2$RIN)
meta2$PMI <- as.numeric(meta2$PMI)
meta2$Race <- as.factor(meta2$Race)
meta2$Sex <- as.factor(meta2$Sex)
meta2$Condition <- as.factor(meta2$Condition)

data2@meta.data$offset <- ratio * data2@meta.data$nCount_RNA + (1-ratio) * data2@meta.data$nFeature_RNA

df = model.matrix(~AgeDeath+Sex+Race+PMI+RIN+Condition, data=meta2)

re = nebula(data2$RNA@counts,data2@meta.data$Sample,pred=df,offset=data2@meta.data$offset,model=model,kappa=kappa,cutoff_cell=cell,ncore=5)
deg <- re$summary
rownames(deg) <- rownames(data2$RNA@counts)

if (cond=='PTSD') {
    deg$padj <- p.adjust(deg$p_ConditionPTSD,method='BH')
    write.table(deg,glue('{savepath}/{celltype}_{cond}_{model}_{ratio}.csv'),sep='\t')
    sig_deg = deg[ which(deg$padj < 0.1 & abs(deg$logFC_ConditionPTSD) > log2(1.2)), ]
    num_up <- nrow(sig_deg[sig_deg$logFC_ConditionPTSD > 0,])
    num_down <- nrow(sig_deg[sig_deg$logFC_ConditionPTSD < 0,])
    EnhancedVolcano(deg,subtitle=NULL,
                lab=rownames(deg),x='logFC_ConditionPTSD',y='padj',pCutoff=0.1,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", nrow(deg), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
    ggsave(glue('{savepath}/{celltype}_{cond}_{model}_{ratio}.png'))
    
} else {
    deg$padj <- p.adjust(deg$p_ConditionMDD,method='BH')
    write.table(deg,glue('{savepath}/{celltype}_{cond}_{model}_{ratio}.csv'),sep='\t')
    sig_deg = deg[ which(deg$padj < 0.1 & abs(deg$logFC_ConditionMDD) > log2(1.2)), ]
    num_up <- nrow(sig_deg[sig_deg$logFC_ConditionMDD > 0,])
    num_down <- nrow(sig_deg[sig_deg$logFC_ConditionMDD < 0,])
    EnhancedVolcano(deg,subtitle=NULL,
                lab=rownames(deg),x='logFC_ConditionMDD',y='padj',pCutoff=0.1,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", nrow(deg), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
    ggsave(glue('{savepath}/{celltype}_{cond}_{model}_{ratio}.png'))
}



