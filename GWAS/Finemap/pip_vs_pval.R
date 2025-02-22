library(glue)
library(ggplot2)
library(dplyr)


### Credit: Xiangyu Zhang ###


data <- read.table("/home/ah2428/ShareZhangLab/PTSD/variantFineMap/Finemap_SNPs.csv", sep="\t", header = T)

totalPCL <- read.table(glue("/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/GWAS/GWAS_PTSD/dbGAP_totalPCL_eur_N"), header=T)
reexperiencing <- read.table(glue("/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/GWAS/GWAS_PTSD/dbGAP_reexperiencing_eur_N"), header=T)

data <- data[order(data$GWAS),]

SNP1 <- data$SNP[which(data$GWAS=='Reexperiencing')]
SNP2 <- data$SNP[which(data$GWAS=='TotalPCL')]
SNP <- append(SNP1, SNP2)

pvalue1 <- reexperiencing[match(SNP1, reexperiencing$SNP_ID), c(1,11)]
pvalue2 <- totalPCL[match(SNP2, totalPCL$SNP_ID), c(1,11)]

pvalue <- rbind(pvalue1, pvalue2)
final_data <- cbind(data, pvalue)

plot_data <- final_data[, c(2,3,4,5,7,8,14,16)]
plot_data$log_p <- -log(plot_data$P) ### log-transformation

write.csv(plot_data,'pip_vs_pval_data.csv',quote=F,row.names=F)