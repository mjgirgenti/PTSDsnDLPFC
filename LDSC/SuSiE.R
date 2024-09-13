library(genio)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(susieR)
library(ggplot2)
library(glue)

### Credit: Xiangyu Zhang ###

args = commandArgs(trailingOnly=TRUE)
gene <- args[1]
gwas <- args[2]
celltype <- args[3]
weight <- as.numeric(args[4]) 


namelist <- read.table("/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/genedata/namelist.csv")

bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm=T)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


totalPCL <- read.table(glue("/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/GWAS/GWAS_PTSD/dbGAP_{gwas}_eur_N"), header=T)
# read in the summary statistics for totalPCL


if (celltype=='EXC') {
    cellt <- 'Ex'
} else if (celltype=='INH') {
    cellt <- 'In'
} else {
    invisible()
}

i <- as.integer(rownames(namelist[namelist$genes==gene,]))
#i <- 14 # ELFN1 for example

chr <- namelist[i,2] # chromosome for ELFN1
gene <- namelist[i,1]
begin <- namelist[i,3]  # Starting position for ELFN1
end <- namelist[i,4]


path <- paste('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/peaks_for_xiangyu_0.4/UNION_intersect_',
              celltype,'.peaks.score.rescaled.corr0.4.bed',sep='')
ctpeaks <- read.table(path, sep="\t")
# read in the ATAC peak information for certain cell type



path <- paste('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/genedata/',gene,'pos.txt',sep='')
snplist <- read.table(path)
## read in the name and position of SNPs in ELFN1:


colnames(snplist) <- c("start","end","SNP_ID")
snplist1 <- merge(snplist,totalPCL,by="SNP_ID")
## get the common SNPs of summary statistics and risk gene


ct <- ctpeaks[(ctpeaks$V1==chr)&(ctpeaks$V5==cellt),]
ctgene <- ct[ct$V2>(begin-1)&ct$V3<(end+1),]

index <- numeric(0)
for (i in 1:nrow(snplist1)) {
  pos <- snplist1$start[i]
  s <- sum((pos>(ctgene$V2-1))&(pos<(ctgene$V3)))
  index <- append(index,s)
}
peaksnp <- snplist1$SNP_ID[index==1]
# select all the SNPs in ELFN1 that also lay in ATAC peaks, with SNP names saved in peaksnp



path <- paste('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/genedata/',gene,'1000g',sep='')
# the path for reference panel
dat <- read_plink(path)
fam <- dat$fam
bim <- dat$bim
Xmatrix <- bedNA(t(dat$X))
Xmatrix <- scale(Xmatrix)

# readin and process the genotype data of reference panel


commonlist <- intersect(snplist$SNP_ID, totalPCL$SNP_ID)
commonlist1 <- intersect(bim$id, commonlist) 
# find the common SNPs of reference panel, summary statistics and risk gene, saved in the commonlist1.


#whichid = which(commonlist1==snp)
#length(commonlist1)
# get the index for the leading GWAS signal


snporder <- match(commonlist1, snplist$SNP_ID)
snplistnew <- snplist[snporder,]
startposition <- snplist$start[snporder]
endposition <- snplist$end[snporder]

indexall <- match(commonlist1, totalPCL$SNP_ID)
summaryfinal <- totalPCL[indexall,]
# make sure the snp names, positions, and summary statistics has the same order with commonlist1



peakcommon <- intersect(peaksnp,commonlist1)
peakindex <- match(peakcommon, commonlist1)
print(length(peakindex))



peakstart <- snplistnew[peakindex,]
peakstart <- peakstart$start
prior <- numeric(length(commonlist1)) + weight
# where you select the prior weights for SNPs outside ATAC peaks


if(length(peakindex)!=0){
  for (j in 1:length(peakindex)) {
    iii <- which((ct$V2<(peakstart[j]+1))&(ct$V3>peakstart[j]))
    prior[peakindex[j]] <- ct$V6[iii]
  }
}
#  prior <- sigmoid(prior-1)
# change the prior weights for SNPs within ATAC peaks


z_scores <- summaryfinal$BETA/summaryfinal$SE   
# prepare the z_scores


correlationgenet <- cor(Xmatrix[,commonlist1])   # prepare the LD matrix
print(dim(correlationgenet))

lambda = estimate_s_rss(z_scores, correlationgenet, n=186689)
#  lambda = 0.95
lambda0 = lambda

I <- diag(1,nrow = nrow(correlationgenet), ncol = ncol(correlationgenet))
corfinal <- (1-lambda)*correlationgenet + lambda*I   
# solve the problem of the inconsistency between summary statistics and LD matrix via regularization.


lambda = estimate_s_rss(z_scores, corfinal, n=186689)
print(lambda)
# can check the balance between the new LD matrix corfinal and summary statistics


# conduct the SuSiE analysis:
fitgenet <- susie_rss(bhat = summaryfinal$BETA, shat = summaryfinal$SE, min_abs_corr = 0.5,
                      R = corfinal, n = 186689, L=5)

fitgenet2 <- susie_rss(bhat = summaryfinal$BETA, shat = summaryfinal$SE, min_abs_corr = 0.5,
                       prior_weights = prior, R = corfinal, n = 186689, L=5)


r1 <- colSums(fitgenet$alpha)
for(i in 1:ncol(fitgenet$alpha)){
  r1[i] <- 1-prod(1-fitgenet$alpha[,i])
}

r2 <- colSums(fitgenet2$alpha)
for(i in 1:ncol(fitgenet2$alpha)){
  r2[i] <- 1-prod(1-fitgenet2$alpha[,i])
}
# get the estimated PIPs for the original SuSiE (r1) and ATAC adjusted SuSiE (r2).


# make the scatter plot:
datapip <- data.frame(pip = c(r1,r2), position = rep(startposition,2), 
                      method = c(rep("original",length(r1)),rep("atacprior",length(r1))))
#  datapoint <- datapip[c(whichid, whichid+length(r1)),]
datacs <- datapip[c(peakindex,peakindex+length(r1)),]
# r1: estimated PIP from SuSiE without using ATAC peak information
# r2: estimated PIP from SuSiE using ATAC peak information
# startposition: start position for SNPs
# whichid: the index of leading SNP in GWAS
# peakindex: the index for SNPs within ATAC peaks

path <- glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/finemap_results_0.1/{gene}.{gwas}.{celltype}.{weight}.plot.pdf')
pdf(path, width=5, height=4)
p <- ggplot(datapip, aes(position, pip)) + 
  geom_point(aes(colour = method))+
  labs(y= "PIP", x = "position of SNP")+
  #    geom_point(data = datapoint, col = 'purple', shape=21, size=3)+
  geom_point(data = datacs, col = 'orange', shape=2, size=3)
p
#  print(p)
dev.off()

# save the result into one dataframe
result <- data.frame(SNP_ID = commonlist1, start = startposition, end = endposition,
                     SuSiE_PIP_O = r1, SuSiE_PIP_A = r2,
                     peak_indicator = prior)
path <- glue('/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/variantFineMap/finemap_results_0.1/{gene}.{gwas}.{celltype}.{weight}.result.txt')
write.table(result, path, quote = F, row.names = F)
