library(MAST)
library(ggplot2)
library(EnhancedVolcano)
library(glue)
library(data.table)
library(stringr)
library(dplyr)
library(stats)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--celltype', type='character', help='celltype')
parser$add_argument('--plotname', type='character', help='plotname')
parser$add_argument('--path', type='character', help='path')
args <- parser$parse_args()

celltype <- args$celltype
plotname <- args$plotname
path <- args$path

deg = read.table(glue('{path}/{celltype}_DESEQ2_DEG.csv'),sep=',',header=1)
sig_deg = deg[ which(deg$pvalue < 0.05 & abs(deg$log2FoldChange) > log2(1.2)), ]
num_up <- nrow(sig_deg[sig_deg$log2FoldChange > 0,])
num_down <- nrow(sig_deg[sig_deg$log2FoldChange < 0,])

EnhancedVolcano(deg,subtitle=NULL,title=glue("{celltype} {plotname} DESEQ2"),
                lab=deg$featurekey,x='log2FoldChange',y='pvalue',pCutoff=0.05,labSize=4,FCcutoff=log2(1.2),
                caption=glue("genes in at least 5% of cells = ", nrow(deg), "\n",
                             "SIG UP = ", num_up, "\n", "SIG DOWN = ", num_down, "\n"),
                legendLabels = c("NS", expression(log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))


ggsave(glue('{path}/{celltype}_DESEQ2_volcano.png'))

