library(RRHO2)
library(glue)
library(ggplot2)
library(dplyr)
library(data.table)
library(argparse)
library(gridExtra)


parser <- ArgumentParser()
parser$add_argument('--celltype', type='character')
args <- parser$parse_args()

celltype <- args$celltype

savepath <- glue('/home/ah2428/ShareZhangLab/ah2428/PTSD/RRHO_results_step_20/')
if (!dir.exists(savepath)) {
  dir.create(savepath, recursive=TRUE)
  cat("Directory created:", savepath, "\n")
} else {
  cat("Directory already exists:", savepath, "\n")
}

ptsd = read.table('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_results_df/snRNA/PTSD_vs_CON/MAST_WILCOX_intersect_all.csv',sep='\t',header=1)
mdd = read.table('/home/ah2428/ShareZhangLab/PTSD/RNA/DEG_results_df/snRNA/MDD_vs_CON/MAST_WILCOX_intersect_all.csv',sep='\t',header=1)


RRHO2_heatmap <- function(RRHO_obj, maximum=NULL, minimum=NULL, colorGradient=NULL, labels=NULL, ...)
{
  
  hypermat <- RRHO_obj$hypermat
  method <- RRHO_obj$method
  
  if(is.null(labels)){
    labels <- RRHO_obj$labels	
  }
	
  if(!is.null(maximum)){
    hypermat[hypermat>maximum] <- maximum
  } else {
    maximum <- max(hypermat,na.rm=TRUE)
  }
    
  if(!is.null(minimum)){
    hypermat[hypermat<minimum] <- minimum
  } else {
    minimum <- min(hypermat,na.rm=TRUE)
  }
    
  if(minimum > maximum){
	  stop("minimum > maximum, please check these function arguments!")
  }
  
  color.bar <- function(lut, min, max=-min, 
                        nticks=11, 
                        ticks=seq(min, max, len=nticks), 
                        title='') {
    scale  <- (length(lut)-1)/(max-min)
    plot(c(0,10), c(min,max), type='n', bty='n', 
         xaxt='n', xlab='', yaxt='n', ylab='')
    mtext(title,2,2.3, cex=2)
    axis(2, round(ticks,0), las=1,cex.lab=1)
    for (i in 1:(length(lut)-1)) {
      y  <- (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
  }
  
	if(is.null(colorGradient)){
	  jet.colors  <- colorRampPalette(
	    c("#00007F", "blue", "#007FFF", "cyan", 
	      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		colorGradient <- jet.colors(101)
	}
  layout(matrix(c(rep(1, 6), 2), 1, 7, byrow = TRUE))
  
  breaks <- seq(minimum,maximum,length.out = length(colorGradient) + 1)
  image(hypermat, col = colorGradient,breaks=breaks,
        axes = FALSE, ...)
  
  if(!is.null(labels)){
    mtext(labels[2],2,2,cex=2)
    mtext(labels[1],1,2,cex=2)
  }
  
  if(method == "hyper"){
    atitle <- ifelse(RRHO_obj$log10.ind, "-log10(p-value)", "-log(pval)")
    color.bar(colorGradient, min = minimum, max = maximum, nticks = 6, title = atitle)
  } else if (method == "fisher"){
    atitle <- "-log10(odds ratio)"
    color.bar(colorGradient, min = minimum, max = maximum, nticks = 6, title = atitle)
  } else {
    stop("internal error (1), please report this error to https://github.com/RRHO2/RRHO2/issues")
  }
  invisible(hypermat)
}

replace_zeros_in_column <- function(df, col_name) {
  if (col_name %in% names(df) && is.numeric(df[[col_name]])) {
    min_nonzero <- min(df[[col_name]][df[[col_name]] > 0], na.rm = TRUE)  # Get min nonzero value
    df[[col_name]][df[[col_name]] == 0] <- min_nonzero  # Replace 0s
  }
  return(df)
}

mdd_deg = mdd[mdd$Celltype==celltype,]
ptsd_deg = ptsd[ptsd$Celltype==celltype,]
mdd_new <- replace_zeros_in_column(mdd_deg, "FDR")
ptsd_new <- replace_zeros_in_column(ptsd_deg, 'FDR') 

# Calculate DDE (Differentially Defined Expression)
ptsd_new$DDE <- -log10(ptsd_new$FDR) * ptsd_new$log2FC
mdd_new$DDE <- -log10(mdd_new$FDR) * mdd_new$log2FC

# Prepare the lists of gene names and DDE values for MDD and PTSD
ptsd_list <- ptsd_new[,c('Genename','DDE')]
mdd_list <- mdd_new[,c('Genename','DDE')]

# Find intersection of genes between PTSD and MDD
inter <- intersect(ptsd_list$Genename, mdd_list$Genename)
ptsd_list <- ptsd_list[ptsd_list$Genename %in% inter,]
mdd_list <- mdd_list[mdd_list$Genename %in% inter,]

# Initialize RRHO object
RRHO_obj <- RRHO2_initialize(ptsd_list, mdd_list, labels = c(glue("PTSD {celltype}"), glue("MDD {celltype}")), log10.ind=TRUE, boundary=0.05, method='fisher', stepsize=20)
pdf(file=glue('{savepath}/{celltype}_oddsratio_heatmap.pdf'),width=6,height=6,bg='white')
RRHO2_heatmap(RRHO_obj)
dev.off()

# Initialize RRHO object
RRHO_obj <- RRHO2_initialize(ptsd_list, mdd_list, labels = c(glue("PTSD {celltype}"), glue("MDD {celltype}")), log10.ind=TRUE, boundary=0.05, method='hyper', stepsize=20)
pdf(file=glue('{savepath}/{celltype}_pval_heatmap.pdf'),width=6,height=6,bg='white')
RRHO2_heatmap(RRHO_obj)
dev.off()

uu_genes <- RRHO_obj$genelist_uu$gene_list_overlap_uu
dd_genes <- RRHO_obj$genelist_dd$gene_list_overlap_dd
ud_genes <- RRHO_obj$genelist_ud$gene_list_overlap_ud
du_genes <- RRHO_obj$genelist_du$gene_list_overlap_du

write.table(ptsd_list,glue('{savepath}/{celltype}_ptsd_list.txt',sep='\t',quote=F,row.names=F))
write.table(mdd_list,glue('{savepath}/{celltype}_mdd_list.txt',sep='\t',quote=F,row.names=F))   
write.table(uu_genes,glue('{savepath}/{celltype}_uu_genes.txt',sep='\t',quote=F,row.names=F))
write.table(dd_genes,glue('{savepath}/{celltype}_dd_genes.txt',sep='\t',quote=F,row.names=F))
write.table(ud_genes,glue('{savepath}/{celltype}_ud_genes.txt',sep='\t',quote=F,row.names=F))
write.table(du_genes,glue('{savepath}/{celltype}_du_genes.txt',sep='\t',quote=F,row.names=F))



