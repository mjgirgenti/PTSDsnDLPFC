library(ggplot2)
library(grDevices)
library(ggthemes)
library(RColorBrewer)


### Credit: Yingxin Lin ###

# Figure 6a: LDSC enrichment

dir <- "/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/unique_peaks/results"

traits <- c('TotalPCL','CaseControl','Reexperiencing','Avoidance','Hyperarousal','Schizophrenia','MVP_Depression',
            'Alzheimer','Alcoholism','IBD','Diabetes','Height','Education')

reordered <- c()
for (t in traits) {
    files <- list.files(dir,pattern=t)
    files <- c(files[3],files[4],files[6],files[7],files[2],files[1],files[5])
    reordered <- append(reordered,files)
}

res <- lapply(reordered, function(x) read.delim(file.path(dir, x)))
names(res) <- gsub(".unique.peaks.results", "", reordered)

res_L2 <- do.call(rbind, lapply(res, function(x) x[1, ]))
res_L2$trait <- do.call(rbind, strsplit(rownames(res_L2), "\\."))[, 1]
res_L2$celltype <- do.call(rbind, strsplit(rownames(res_L2), "\\."))[, 2]  
                                
res_L2$trait <- factor(res_L2$trait,levels=traits)
res_L2$celltype <- factor(res_L2$celltype,levels=rev(unique(res_L2$celltype)))

rdbu <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
ggplot(res_L2, aes(x = trait, y = celltype, color = Coefficient_z.score, size = -log10(Enrichment_p))) +
  geom_point() +
  scale_colour_gradientn(colors = rdbu, limits = c(-3, 3), oob = scales::squish) +
  theme_bw() +
  theme(aspect.ratio = 0.5) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle=90,vjust=0.5)) 

ggsave('/home/ah2428/LDSC.pdf')
dev.off()



