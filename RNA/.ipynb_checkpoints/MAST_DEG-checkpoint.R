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

args = commandArgs(trailingOnly=TRUE)

cluster <- args[1] # celltype
cond <- args[2] # condition (PTSD,MDD)

folder <- glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/RNA/results/062522_DEG/subclass/mast/{cond}_vs_CON')
data <- LoadH5Seurat(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/RNA/data/DEG_data_subclass_062522/{cond}_vs_CON/{cond}_CON_{cluster}.h5seurat'),assays='RNA')
meta <- read.delim2(glue('/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/RNA/data/DEG_data_subclass_062522/{cond}_vs_CON/{cond}_CON_{cluster}_meta.txt'))

print(cluster)
print(cond)

# if exceed R memory 2^31
if (dim(data)[2] > 76745) {
    print('True')
    mat <- data$RNA@counts
    tmp <- matrix(data = 0, nrow=nrow(mat), ncol=ncol(mat))
    row_pos <- mat@i+1
    col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
    val <- mat@x
    for (i in seq_along(val)){
        tmp[row_pos[i],col_pos[i]] <- val[i]
    }
    rownames(tmp) <- rownames(mat)
    colnames(tmp) <- colnames(mat)
    dge <- DGEList(counts = tmp)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- cpm(dge)
    bin_counts <- as.matrix((tmp > 0) + 0)
} else {
    counts <- data$RNA@counts
    dge <- DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- cpm(dge)
    bin_counts <- as.matrix((counts > 0) + 0)
}

logcpms <- log2(cpms+1)

data@meta.data$wellKey <- rownames(data@meta.data)
genes <- data.frame(rownames(data$RNA@counts))
colnames(genes) <- c('primerid')
sca <- FromMatrix(logcpms,data@meta.data,genes)

num_con <- floor(dim(meta[meta$Condition=='CON',])[1]*0.05)
num_cond <- floor(dim(meta[meta$Condition==cond,])[1]*0.05)

con_sums <- rowSums(bin_counts[, which(meta$Condition=='CON')])
cond_sums <- rowSums(bin_counts[, which(meta$Condition==cond)])
df <- data.frame(con=con_sums,cond=cond_sums)
filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
sca <- sca[filtered_genes,]

meta <- transform(meta, AgeDeath = as.numeric(AgeDeath))
meta <- transform(meta, PMI = as.numeric(PMI))
meta <- transform(meta, RIN = as.numeric(RIN))

colData(sca)$cdr <- scale(colSums(assay(sca)>0))
colData(sca)$condition <- as.factor(meta$Condition)
colData(sca)$age <- as.numeric(meta$AgeDeath)
colData(sca)$pmi <- as.numeric(meta$PMI)
colData(sca)$rin <- as.numeric(meta$RIN)
colData(sca)$sex <- as.factor(meta$Sex)
colData(sca)$race <- as.factor(meta$Race)

print('starting MAST')
zlm_output <- zlm(~ condition + cdr + age + pmi + rin + sex + race, sca)
print('finished MAST')

summaryCond <- summary(zlm_output, doLRT=glue('condition{cond}')) 
summaryDt <- summaryCond$datatable
deg <- merge(summaryDt[contrast==glue('condition{cond}') & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==glue('condition{cond}') & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

deg[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

deg <- merge(deg[fdr<1 & abs(coef)>0], as.data.table(mcols(sca)), by='primerid')
setorder(deg, fdr)
write.table(deg,glue('{folder}/{cluster}_DEG.csv'),sep='\t')

sig_deg = deg[ which(deg$fdr < 0.01 & abs(deg$coef) > log2(1.2)), ]
write.table(sig_deg,glue('{folder}/{cluster}_SIG_DEG.csv'),sep='\t')


# check consistency of direction 
print('check consistency')

mat <- cpms
data@meta.data$condition <- meta$Condition
df <- data@meta.data['condition']
rownames(df) <- 1:nrow(df)
con_ind <- which(df$condition=='CON')
cond_ind <- which(df$condition==cond)
sig_down <- sig_deg[sig_deg$coef < 0,]
sig_up <- sig_deg[sig_deg$coef > 0,]

count <- 0
for (i in 1:dim(sig_up)[1]) {
    gene <- sig_up$primerid[i]
    df[[i]] <- mat[gene,]
    con_mean <- mean(df[[i]][con_ind])
    cond_mean <- mean(df[[i]][cond_ind])
    if ((sig_up[sig_up$primerid==gene,]$coef > 0) & (cond_mean - con_mean > 0)) {
        count <- count + 1
    }
}
consistent_up <- round(count/(dim(sig_up)[1]),3)
print(glue('Consistent up: ', {consistent_up}))

count <- 0
for (i in 1:dim(sig_down)[1]) {
    gene <- sig_down$primerid[i]
    df[[i]] <- mat[gene,]
    con_mean <- mean(df[[i]][con_ind])
    cond_mean <- mean(df[[i]][cond_ind])
    if ((sig_down[sig_down$primerid==gene,]$coef < 0) & (cond_mean - con_mean < 0)) {
        count <- count + 1
    }
}
consistent_down <- round(count/(dim(sig_down)[1]),3)
print(glue('Consistent down: ', {consistent_down}))


EnhancedVolcano(deg,lab=deg$primerid,x='coef',y='fdr',pCutoff=0.01,title=glue("{cluster} {cond} vs CON MAST"),labSize=4,FCcutoff=log2(1.2),caption=glue("genes in at least 5% of cells = ", nrow(deg), " sig degs = ", nrow(sig_deg), "\nconsistent down = ",consistent_down, " consistent up = ", consistent_up),legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))

ggsave(glue('{folder}/{cluster}_volcano.png'))

