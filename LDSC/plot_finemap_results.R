library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(Seurat)
library(GenomicRanges)
library(HDF5Array)
library(SummarizedExperiment)
library("tidyr")


# Figure 6c-f and Extended Data Figure 8

args = commandArgs(trailingOnly=TRUE)

addArchRThreads(threads = 20) 

proj <- loadArchRProject('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL')

union_peaks <- read.table('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/peaks/union_peaks_meta.bed',sep='\t',header=1)

gene <- args[1]
gwas <- args[2]
celltype <- args[3]
chromstart <- as.integer(args[4])
chromend <- as.integer(args[5])
tad1 <- as.integer(args[6])
tad2 <- as.integer(args[7])
corr <- as.numeric(args[8])


# Xiangyu's results
dataset <- read.table(glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/finemap_results_0.1/{gene}.{gwas}.{celltype}.0.1.result.txt'), header = T)

gene_df = data.frame(getGenes(proj))
gene_df = na.omit(gene_df[gene_df$symbol==gene,])

chrom <- as.character(gene_df$seqnames)
#chromstart <- min(dataset$start)
#chromend <- max(dataset$end)

# get enhancer loops anchored to TSS of gene
loops <- getPeak2GeneLinks(proj,resolution = NULL, corCutOff=corr, FDRCutOff=0.05, returnLoops = T)
loops <- data.frame(loops$Peak2GeneLinks)
loops <- loops[loops$seqnames==chrom & loops$start > chromstart & loops$end < chromend,]

loops <- loops[,c('seqnames','start','end','value')]
dir <- gene_df$strand
if (dir=='-') {
    loops <- rbind(loops[loops$start==gene_df$end,],loops[loops$end==gene_df$end,])
} else {
    loops <- rbind(loops[loops$start==gene_df$start,],loops[loops$end==gene_df$start,])
}

# get peaks linked to gene
p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF[c('chr', 'start', 'end')] <- str_split_fixed(p2geneDF$peakName, '_', 3)
peaks <- p2geneDF[p2geneDF$geneName==gene & p2geneDF$Correlation > corr & p2geneDF$FDR < 0.05 & p2geneDF$start > chromstart & p2geneDF$end < chromend & p2geneDF$VarQATAC > 0.25 & p2geneDF$VarQRNA > 0.25,]

# get celltype specific union peaks
df <- union_peaks[,c('peak','celltype')]
colnames(df)[1] <- 'peakName'
peaks <- merge(peaks,df,on='peakName')
if (celltype=='EXC') {
    peaks <- peaks[peaks$celltype=='Ex',]
} else if (celltype=='INH') {
    peaks <- peaks[peaks$celltype=='In',]
} else {
    invisible()
}

# Xiangyu's results
r1 <- dataset$SuSiE_PIP_O # CHANGE HERE
r2 <- dataset$SuSiE_PIP_A # CHANGE HERE
startposition <- dataset$start # CHANGE HERE
snp <- dataset$SNP_ID
whichid <- which(dataset$SNP_ID==snp)

datapip <- data.frame(pip = c(r1,r2), snp = snp, position = rep(startposition,2),
                      method = c(rep("original",length(r1)),rep("atacprior",length(r1))))
datapoint <- datapip[c(whichid, whichid+length(r1)),]

datapip_atacprior = datapip[datapip$method == "atacprior",]
datapoint_atacprior = datapoint[datapoint$method == "atacprior",]

# get peaks with variants, get variant positions
pip_peaks <- c()
pip_pos <- c()
for (pos in datapip_atacprior[datapip_atacprior$pip > 0.1,]$position) {
    p <- peaks[peaks$start < pos & peaks$end > pos,]$peakName
    if (length(p) >= 1) {
        print(datapip_atacprior[datapip_atacprior$position==pos,])
        pip_pos <- append(pip_pos,pos)
        pip_peaks <- append(pip_peaks,p)
    }
}
pip_peaks <- data.frame(pip_peaks)
pip_peaks[c('chr','start', 'end')] <- str_split_fixed(pip_peaks$pip_peaks, '_', 3)

# get loops that connect TSS of gene to peaks with variants
ls <- c()
for (e in loops$end) {
    if (dim(pip_peaks[pip_peaks$start < e & pip_peaks$end > e,])[1] >= 1) {
        print(pip_peaks[pip_peaks$start < e & pip_peaks$end > e,])
        ls <- rbind(ls,loops[loops$end==e,])
    }
}
ls2 <- c()
for (e in loops$start) {
    if (dim(pip_peaks[pip_peaks$start < e & pip_peaks$end > e,])[1] >= 1) {
        print(pip_peaks[pip_peaks$start < e & pip_peaks$end > e,])
        ls2 <- rbind(ls2,loops[loops$start==e,])
    }
}
ls_both <- rbind(ls,ls2)

sapply(list.files('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/misc/girgenti_project/ArchR/R',full.name=T),source)
source('../utils/ArchRBrowser-Variant.R')

if (celltype=='EXC') {
    coverage_color <- '#b22222'
} else if (celltype=='INH') {
    coverage_color <- '#2E8B57'
} else {
    coverage_color <- 'NONE'
}

# plot track 
p <- plotBrowserTrackVariant(
    ArchRProj = proj, 
    groupBy = "Int_Cluster", 
    useGroups = c(glue('{celltype}')),
    region = GRanges(seqnames = chrom, ranges = IRanges(start=chromstart,end=chromend)),
    variants1 = TRUE,
    pip_pos = pip_pos,
    gwas = gwas,
    gene = gene,
    celltype = celltype,
    chromstart = chromstart,
    chromend = chromend,
    tad1 = tad1,
    tad2 = tad2,
    loops = GRanges(ls_both),
    ylim2 = c(0, 1),
    plotSummary = c("bulkTrack","loopTrack","variantTrack1","geneTrack"),
    sizes = c(1, 1, 1, 1),
    tileSize = 3000,
    verbose = FALSE,
    pal = c(coverage_color,'gray'),
    baseSize = 14,
    tickWidth = 1,
    borderWidth = 1
)

plotPDF(p, 
    name = glue("{gene}_with_tads"), 
    ArchRProj = proj, 
    addDOC = FALSE, width = 8, height = 3)
