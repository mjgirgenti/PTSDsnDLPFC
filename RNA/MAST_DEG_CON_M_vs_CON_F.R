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
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--celltype', type='character', help='celltype')
parser$add_argument('--condition', type='character', help='condition', default='M')
parser$add_argument('--datapath', type='character', help='datapath')
parser$add_argument('--savepath', type='character', help='savepath')
args <- parser$parse_args()

celltype <- args$celltype
cond <- args$condition
datapath <- args$datapath
savepath <- args$savepath
print(celltype)
print(cond)

data <- LoadH5Seurat(glue('{datapath}/{cond}_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
meta <- fread(glue('{datapath}/{cond}_CON_{celltype}_meta.txt'))

protein <- fread('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt',header=FALSE)$V2
protein_genes <- intersect(rownames(data),protein)
data <- data[protein_genes,]
print(dim(data))

counts <- data$RNA@counts
counts <- as(counts, "dgCMatrix")
dge <- DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
cpms <- cpm(dge)
bin_counts <- (counts > 0) 
logcpms <- log2(cpms+1)

data@meta.data$wellKey <- rownames(data@meta.data)
genes <- data.frame(rownames(data$RNA@counts))
colnames(genes) <- c('primerid')
sca <- FromMatrix(logcpms,data@meta.data,genes)

num_con <- floor(dim(meta[meta$Sex=='F',])[1]*0.05)
num_cond <- floor(dim(meta[meta$Sex==cond,])[1]*0.05)

con_sums <- rowSums(bin_counts[, which(meta$Sex=='F')])
cond_sums <- rowSums(bin_counts[, which(meta$Sex==cond)])
df <- data.frame(con=con_sums,cond=cond_sums)
filtered_genes <- rownames(df[df$con >= num_con | df$cond >= num_cond,])
sca <- sca[filtered_genes,]
print(dim(sca))

colData(sca)$cdr <- scale(colSums(assay(sca)>0))
colData(sca)$condition <- as.factor(meta$Condition)
colData(sca)$age <- as.numeric(meta$AgeDeath)
colData(sca)$pmi <- as.numeric(meta$PMI)
colData(sca)$rin <- as.numeric(meta$RIN)
colData(sca)$sex <- as.factor(meta$Sex)
colData(sca)$sex <-relevel(colData(sca)$sex,'F')
colData(sca)$race <- as.factor(meta$Race)

print('starting MAST')
zlm_output <- zlm(~ sex + cdr + age + pmi + rin + race, sca)
print('finished MAST')

summaryCond <- summary(zlm_output, doLRT=glue('sex{cond}')) 
summaryDt <- summaryCond$datatable
deg <- merge(summaryDt[contrast==glue('sex{cond}') & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==glue('sex{cond}') & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

deg[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

deg <- merge(deg[fdr<1 & abs(coef)>0], as.data.table(mcols(sca)), by='primerid')
setorder(deg, fdr)
write.table(deg,glue('{savepath}/{celltype}_MAST_DEG.csv'),sep='\t')

sig_deg = deg[ which(deg$fdr < 0.01 & abs(deg$coef) > log2(1.2)), ]

num_up <- nrow(sig_deg[sig_deg$coef > 0,])
num_down <- nrow(sig_deg[sig_deg$coef < 0,])

EnhancedVolcano(deg,subtitle=NULL,title=glue("{celltype} CON M vs CON F MAST"),
                lab=deg$primerid,x='coef',y='fdr',pCutoff=0.01,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", length(filtered_genes), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))

ggsave(glue('{savepath}/{celltype}_MAST_volcano.png'))