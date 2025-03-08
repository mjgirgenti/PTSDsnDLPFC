library(enrichR)
library(ggplot2)



### Credit: Siwei Xu ### 

setEnrichrSite("Enrichr") # Human genes

dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023")


# 502 PTSD only DEG GO

degs <- read.table('/home/ah2428/ShareZhangLab/PTSD/RNA/502_PTSD_only_DEGs.csv', sep = '\t', header = 1)$X0
enriched <- enrichr(degs, dbs)

mf <- enriched$GO_Molecular_Function_2023
Name <- 'GO Molecular Function 2023'
mf <- cbind(Name,mf)

bp <- enriched$GO_Biological_Process_2023
Name <- 'GO Biological Process 2023'
bp <- cbind(Name,bp)

df <- rbind(mf,bp)
df <- df[order(df$Adjusted.P.value),]

write.table(df,'/home/ah2428/tables/502_PTSD_only_GO.txt',quote=F,row.names=F,sep='\t')



# 1184 PTSD snDEG GO 

genes <- degs <- read.table('/home/ah2428/ShareZhangLab/PTSD/RNA/1184_PTSD_snDEGs.csv', sep = '\t', header = 1)$X0
enriched <- enrichr(degs, dbs)

mf <- enriched$GO_Molecular_Function_2023
Name <- 'GO Molecular Function 2023'
mf <- cbind(Name,mf)

bp <- enriched$GO_Biological_Process_2023
Name <- 'GO Biological Process 2023'
bp <- cbind(Name,bp)

df <- rbind(mf,bp)
df <- df[order(df$Adjusted.P.value),]

write.table(df,'/home/ah2428/tables/1184_PTSD_snDEG_GO.txt',quote=F,row.names=F,sep='\t')
