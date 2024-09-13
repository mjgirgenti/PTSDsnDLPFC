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
parser$add_argument('--condition', type='character', help='condition')
parser$add_argument('--datapath', type='character', help='datapath')
parser$add_argument('--savepath', type='character', help='savepath')
parser$add_argument('--cpm', type='character', help='cpm', default='F')
parser$add_argument('--covs', type='character', help='covs', default='all')
args <- parser$parse_args()

celltype <- args$celltype
cond <- args$condition
datapath <- args$datapath
savepath <- args$savepath
covs <- args$covs
print(celltype)
print(cond)

data <- LoadH5Seurat(glue('{datapath}/{cond}_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
meta <- read.delim2(glue('{datapath}/{cond}_CON_{celltype}_meta.txt'))

protein <- read.table('/home/ah2428/ShareZhangLab/CommonData/protein_coding.genes.with.chr.txt')$V2
protein_genes <- intersect(rownames(data),protein)
data <- data[protein_genes,]
print(dim(data))

data_size <- as.double(dim(data)[1]) * as.double(dim(data)[2])
lim <- 2^31
print(data_size < lim)

# if exceed R memory 2^31
if (args$cpm == 'T') {
    genes = read.csv(glue('{datapath}/EXC_genes.csv'),header=0)$V1
    data <- data[genes,]
    cpms = fread(glue('{datapath}/EXC_cpms.csv.gz'))
    cpms <- as.matrix(cpms)
} else {
    counts <- data$RNA@counts
    dge <- DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- cpm(dge)
    bin_counts <- as.matrix((counts > 0) + 0)
}

counts <- data$RNA@counts
bin_counts <- as.matrix((counts > 0) + 0)
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
colData(sca)$condition <- relevel(colData(sca)$condition,'CON')
colData(sca)$age <- as.numeric(meta$AgeDeath)
colData(sca)$pmi <- as.numeric(meta$PMI)
colData(sca)$rin <- as.numeric(meta$RIN)
colData(sca)$sex <- as.factor(meta$Sex)
colData(sca)$race <- as.factor(meta$Race)

print('starting MAST')
if (covs=='all') {
    print('Running with all covariates')
    zlm_output <- zlm(~ condition + cdr + age + pmi + rin + sex + race, sca)
} else if (covs=='cdr') {
    print('Running with cdr covariate only')
    zlm_output <- zlm(~ condition + cdr, sca)
} else if (covs=='age') {
    print('Running with age covariate only')
    zlm_output <- zlm(~ condition + age, sca)
} else if (covs=='pmi') {
    print('Running with pmi covariate only')
    zlm_output <- zlm(~ condition + pmi, sca)
} else if (covs=='rin') {
    print('Running with rin covariate only')
    zlm_output <- zlm(~ condition + rin, sca)
} else if (covs=='sex') {
    print('Running with sex covariate only')
    zlm_output <- zlm(~ condition + sex, sca)
} else if (covs=='race') {
    print('Running with race covariate only')
    zlm_output <- zlm(~ condition + race, sca)
} else {
    print('Running without covariates')
    zlm_output <- zlm(~ condition, sca)
}
print('finished MAST')

summaryCond <- summary(zlm_output, doLRT=glue('condition{cond}')) 
summaryDt <- summaryCond$datatable
deg <- merge(summaryDt[contrast==glue('condition{cond}') & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==glue('condition{cond}') & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

deg[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

deg <- merge(deg[fdr<1 & abs(coef)>0], as.data.table(mcols(sca)), by='primerid')
setorder(deg, fdr)
write.table(deg,glue('{savepath}/{celltype}_MAST_DEG.csv'),sep='\t')

sig_deg = deg[ which(deg$fdr < 0.01 & abs(deg$coef) > log2(1.2)), ]
write.table(sig_deg,glue('{savepath}/{celltype}_MAST_SIG_DEG.csv'),sep='\t')

num_up <- nrow(sig_deg[sig_deg$coef > 0,])
num_down <- nrow(sig_deg[sig_deg$coef < 0,])

EnhancedVolcano(deg,subtitle=NULL,title=glue("{celltype} {cond} vs CON MAST"),
                lab=deg$primerid,x='coef',y='fdr',pCutoff=0.01,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", length(filtered_genes), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))

ggsave(glue('{savepath}/{celltype}_MAST_volcano.png'))

