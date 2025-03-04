library(enrichR)
library(ggplot2)



### Credit: Siwei Xu ### 

setEnrichrSite("Enrichr") # Human genes

dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023")


# Fig. 2e 

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

df <- df[1:13,]

plotEnrich(df, showTerms = 13, numChar = 60, y = "Count", orderBy = "Adjusted.P.value")
ggsave('/home/ah2428/502_GO.pdf',dpi=300)


# Extended Data Fig. 3b

genes <- degs <- read.table('/home/ah2428/ShareZhangLab/PTSD/RNA/1184_PTSD_snDEGs.csv', sep = '\t', header = 1)$X0
enriched <- enrichr(degs, dbs)

mf <- enriched$GO_Molecular_Function_2023
mf <- mf[order(mf$Adjusted.P.value), ]
mf$name <- 'MF'
mf <- mf[1:10,]

bp <- enriched$GO_Biological_Process_2023
bp <- bp[order(bp$Adjusted.P.value), ]
bp$name <- 'BP'
bp <- bp[1:10,]

df <- rbind(mf,bp)
df <- df[order(df$Adjusted.P.value),]

df <- df[df$Term !='Ubiquitin-Like Protein Ligase Binding (GO:0044389)',]

plotEnrich(df, showTerms = 20, numChar = 60, y = "Count", orderBy = "Adjusted.P.value")
ggsave('/home/ah2428/1184_GO.pdf',dpi=300)

