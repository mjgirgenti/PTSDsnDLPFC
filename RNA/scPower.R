library(scPower)
library(ggplot2)
library(reshape2)
library(reticulate)
library(Seurat)
library(SeuratDisk)
library(glue)
library(data.table)
library(stringr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--celltype', type='character', help='celltype')
parser$add_argument('--condition', type='character', help='condition')
args <- parser$parse_args()

celltype <- args$celltype
cond <- args$condition
print(celltype)
print(cond)

datapath <- glue('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_data/{cond}_vs_CON')
data <- LoadH5Seurat(glue('{datapath}/{cond}_CON_{celltype}.h5seurat'),assays='RNA',meta.data=FALSE)
meta <- read.delim2(glue('{datapath}/{cond}_CON_{celltype}_meta.txt'))
allmeta <- read.csv('/home/ah2428/ShareZhangLab/PTSD/RNA/data/RNA_FINAL_meta.csv',sep='\t',header=1)

nsamples <- length(unique(meta$Channel))
countMatrix <- data$RNA@counts
data@meta.data$channel <- meta$Channel
data@meta.data$celltype <- celltype
data@meta.data$cell <- rownames(data@meta.data)
colnames(data@meta.data) <- c('individual','cell_type','cell')

expressed.genes.df <- NULL

count.matrix <- countMatrix
meta <- data@meta.data

condition <- unique(meta$condition)
print(all(colnames(count.matrix) == rownames(meta)))
#Reformat count matrix into pseudobulk matrix
pseudo.bulk <- create.pseudobulk(count.matrix, meta, colName='cell_type')
#Calculate expressed genes in the pseudobulk matrix
expressed.genes <- calculate.gene.counts(pseudo.bulk, min.counts=3, perc.indiv = 0.5)
#Get the number of expressed genes
num.expressed.genes <- nrow(expressed.genes)

#Save expressed genes
expressed.genes.df <- rbind(expressed.genes.df, data.frame(
            num.cells=ncol(count.matrix), expressed.genes=num.expressed.genes))
print(expressed.genes.df)
print(dim(pseudo.bulk))

norm.mean.values<-NULL
disp.param<-NULL

temp <- nbinom.estimation(as.matrix(countMatrix), sizeFactorMethod = "poscounts")

#Save the normalized mean values
norm.mean.values.temp <- temp[[1]]
norm.mean.values <- rbind(norm.mean.values,norm.mean.values.temp)

#Save the parameter of the mean-dispersion function
disp.param.temp <- temp[[3]]
disp.param <- rbind(disp.param,disp.param.temp)

#Number of cells per cell type as censoring point
censoredPoint <- 1 / ncol(countMatrix)
gamma.fits <- mixed.gamma.estimation(norm.mean.values$mean, num.genes.kept = 28000, censoredPoint = censoredPoint)
print(gamma.fits)
if (is.na(gamma.fits$p1)) {
    gamma.fits$p1 <- 1e-10
}
#Estimate the mean umi values per cell for each matrix
umi.values <- meanUMI.calculation(as.matrix(countMatrix))

gamma.fits <- merge(gamma.fits,umi.values)
gamma.fits <- convert.gamma.parameters(gamma.fits)

#Fit relationship between gamma parameters and UMI values
#gamma.linear.fit.new <- umi.gamma.relation(gamma.fits)
disp.fun.general.new <- dispersion.function.estimation(disp.param)

gamma.fits$ct <- celltype
disp.fun.general.new$ct <- celltype

freq <- dim(meta)[1]/dim(allmeta)[1]
powers <- c()
for (n in c(100,200,500,1000,2000,3000,4000,5000)) {
    print(n)
    power<-power.sameReadDepth.restrictedDoublets(nSamples=nsamples,nCells=n,
    ct.freq=freq,type="de",
    ref.study=scPower::de.ref.study,
    ref.study.name="Blueprint (CLL) iCLL-mCLL",
    cellsPerLane=20000,
    gamma.parameters = gamma.fits,
    ct=celltype,
    disp.fun.param=disp.fun.general.new,
    mappingEfficiency = 0.8,
    min.UMI.counts = 3,
    perc.indiv.expr = 0.5,
    sign.threshold = 0.01,
    MTmethod="FDR")
    powers <- append(powers,power$power)
}
write.table(powers,glue('/home/ah2428/ShareZhangLab/PTSD/RNA/scPower/{celltype}_{cond}_power.txt'))


