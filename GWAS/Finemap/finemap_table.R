library(ArchR)
library(stringr)
library(ggplot2)
library(dplyr)
library(glue)
library(parallel)
library(GenomicRanges)
library("tidyr")

args = commandArgs(trailingOnly=TRUE)
gene <- args[1]
celltype <- args[2]
gwas <- args[3]

namelist <- read.table("/home/ah2428/palmer_scratch/PTSD/FINEMAP/namelist_new2.csv")

chromstart <- namelist[namelist$genes==gene,]$starts
chromend <- namelist[namelist$genes==gene,]$ends

addArchRThreads(threads = 20) 

proj <- loadArchRProject('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/ATAC_FINAL')


r2_data = read.table(glue('/home/ah2428/palmer_scratch/PTSD/FINEMAP/finemap_final/{gene}.{gwas}.{celltype}.0.1.r2.txt'))

union_peaks <- read.table('/home/ah2428/ShareZhangLab/PTSD/ATAC/data/peaks/union_peaks_meta.bed',sep='\t',header=1)

gene_df = data.frame(getGenes(proj))
gene_df = na.omit(gene_df[gene_df$symbol==gene,])
chr <- as.character(gene_df$seqnames)

loops <- getPeak2GeneLinks(proj,resolution = NULL, corCutOff=0.4, FDRCutOff=0.05, returnLoops = T)
loops <- data.frame(loops$Peak2GeneLinks)
loops <- loops[loops$seqnames==chr & loops$start > chromstart & loops$end < chromend,]

loops <- loops[,c('seqnames','start','end','value')]
dir <- gene_df$strand
if (dir=='-') {
    loops <- rbind(loops[loops$start==gene_df$end,],loops[loops$end==gene_df$end,])
} else {
    loops <- rbind(loops[loops$start==gene_df$start,],loops[loops$end==gene_df$start,])
}
loops

p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF[c('chr', 'start', 'end')] <- str_split_fixed(p2geneDF$peakName, '_', 3)
peaks <- p2geneDF[p2geneDF$geneName==gene & p2geneDF$Correlation > 0.4 & p2geneDF$FDR < 0.05 & p2geneDF$chr==chr & p2geneDF$start > chromstart & p2geneDF$end < chromend & p2geneDF$VarQATAC > 0.25 & p2geneDF$VarQRNA > 0.25,]

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

dataset <- read.table(glue('/home/ah2428/palmer_scratch/PTSD/FINEMAP/finemap_final/{gene}.{gwas}.{celltype}.0.1.result.txt'),header=T)

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

pip_peaks <- c()
pip_pos <- c()
for (pos in datapip_atacprior[datapip_atacprior$pip > 0.05,]$position) {
    p <- peaks[peaks$start < pos & peaks$end > pos,]$peakName
    if (length(p) >= 1) {
        #print(datapip_atacprior[datapip_atacprior$position==pos,])
        pip_pos <- append(pip_pos,pos)
        pip_peaks <- append(pip_peaks,p)
    }
}
pip_peaks <- data.frame(pip_peaks)
pip_peaks[c('chr','start', 'end')] <- str_split_fixed(pip_peaks$pip_peaks, '_', 3)

df <- datapip_atacprior[datapip_atacprior$pos %in% pip_pos,]
df <- df[,c('snp','pip','position')]
df$peakName <- pip_peaks$pip_peaks

rescaled <- read.table('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/peaks_for_xiangyu_0.4/UNION.peaks.score.rescaled.corr0.4.bed',sep='\t')

df2 <- c()
df3 <- c()
scores <- c()
cres <- df$peakName
for (cre in cres) {
    x <- union_peaks[union_peaks$peak==cre,]
    x <- x[,c('peak','celltype','peakType','score')]
    df2 <- rbind(df2,x)
    
    y <- peaks[peaks$peakName==cre,]
    y <- y[,c('peakName','Correlation','FDR')]
    df3 <- rbind(df3,y)
    
    score <- rescaled[rescaled$V4==cre,]$V6
    scores <- append(scores,score)
}

colnames(df2)[colnames(df2)=='peak']='peakName'
df3 <- data.frame(df3)


df2 <- cbind(df,df2)
df3[order(match(df3$peakName,df2$peakName)),]

snps_meta <- cbind(df2,df3[,c('Correlation','FDR')])
snps_meta['rescaled'] = scores
snps_meta['gene'] = gene
snps_meta['gwas'] = gwas
snps_meta <- snps_meta[,c('gene','celltype','gwas','snp','pip','position','peakName','peakType','score','rescaled','Correlation','FDR')]

lead_snp <- snps_meta[which.max(snps_meta$pip), ]$snp

snp_r2s <- c()
for (snp in snps_meta$snp) {
    snp_r2 <- r2_data[lead_snp,snp]
    print(snp_r2)
    snp_r2s <- append(snp_r2s,snp_r2)
}

snps_meta$r2 <- snp_r2s
snps_meta <- snps_meta[order(-snps_meta$pip),]

write.table(snps_meta,glue('/home/ah2428/palmer_scratch/PTSD/FINEMAP/finemap_tables/{gene}_{celltype}_{gwas}_SNPs_meta.csv'),row.names=F,quote=F,sep='\t')


