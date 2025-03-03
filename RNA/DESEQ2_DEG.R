library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(tidyverse)
library(reshape2)
library(edgeR)
library(glue)

parser <- ArgumentParser()
parser$add_argument('--condition', type='character')
parser$add_argument('--celltype', type='character')
args <- parser$parse_args()

condition <- args$condition
celltype <- args$celltype

# Load data and metadata
seurat_object <- LoadH5Seurat(glue("../PTSD/RNA/DEG_data/{condition}_vs_CON/{condition}_CON_{celltype}.h5seurat"), assays = "RNA", 
                              meta.data = FALSE, misc = FALSE) 
meta <- read.csv(file = glue("../PTSD/RNA/DEG_data/{condition}_vs_CON/{condition}_CON_{celltype}_meta.txt"),
                 sep = '\t',row.names='barcodekey')
seurat_object <- AddMetaData(object = seurat_object, metadata = meta)
seurat_object$samples <- paste0(seurat_object$Condition, '.', seurat_object$Channel)
seurat_object <- NormalizeData(seurat_object)
celltype_to_use <- unique(seurat_object$class)

# Assign sample names to object
Idents(seurat_object) <- paste(seurat_object$Condition, seurat_object$class, sep = "_")
cts <- AggregateExpression(seurat_object, 
                           group.by = c("class", "samples"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)
cts <- cts$RNA

# Split into a list of PseudoBulk Matrix per cell type
cts.t <- t(cts)  # transpose
cts.t <- as.data.frame(cts.t)  # convert to data.frame
splitRows <- gsub('_.*', '', rownames(cts.t))  # get values where to split
cts.split <- split.data.frame(cts.t, f = factor(splitRows))  # split data.frame

# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
    rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
    t(x)
})

counts_cellType <- cts.split.modified

samples_table <- table(seurat_object@meta.data$class, seurat_object@meta.data$Channel)
samples_df <- as.data.frame(samples_table)
samples_df <- dcast(samples_df, Var1 ~ Var2, value.var = "Freq")
rownames(samples_df) <- samples_df$Var1
samples_df$Var1 <- NULL

samples_cellType <- as.vector(samples_df)

remove_samples <- c()
for (i in 1:length(samples_cellType)) {
    sample = samples_cellType[i]
    if (samples_cellType[[i]] < 30) {
        remove_samples <- c(remove_samples, names(samples_cellType)[[i]])
    }
}

if ('MDD' %in% seurat_object$Condition) {
    sums <- sort(rowSums(cts.t))
    sums <- data.frame(sums)
    sums$cond <- substr(rownames(sums), 0, 3)
    sums <- sums[sums$cond == 'MDD',]
    samples <- rownames(sums)[1:7]

    sums <- sort(rowSums(cts.t), decreasing = TRUE)
    sums <- data.frame(sums)
    sums$cond <- substr(rownames(sums), 0, 3)
    sums <- sums[sums$cond == 'CON',]
    samples_con <- rownames(sums)[1:5]

    remove_samples <- unique(c(remove_samples, samples, samples_con))
}

if (!is.null(remove_samples)) {
    print("removing samples...")
    columns_to_remove <- grep(paste(remove_samples, collapse = "|"), names(counts_cellType))
    print(columns_to_remove)
    counts_cellType_filtered <- counts_cellType[-columns_to_remove]
} else {
    counts_cellType_filtered <- counts_cellType
}

counts_cellType_filtered <- as.data.frame(counts_cellType_filtered)

dge <- DGEList(counts = counts_cellType_filtered)
cpm_matrix <- cpm(dge)
cpm_percent <- rowSums(cpm_matrix<1) / ncol(cpm_matrix)  # % of samples that have <0.5 CPM value for each gene
remove_genes <- names(cpm_percent[cpm_percent > 1])
counts_cellType_filtered <- counts_cellType_filtered[!(rownames(counts_cellType_filtered) %in% remove_genes), ]

# filter
keep <- rowSums(counts_cellType_filtered) >= 20  ## PTSD: keeping only genes that have 20 reads and more
counts_cellType_filtered <- counts_cellType_filtered[keep,]

colData <- data.frame(samples = colnames(counts_cellType_filtered))
colData <- colData %>%
mutate(condition = ifelse(grepl(condition, samples), condition, 'CON')) %>%
column_to_rownames(var = 'samples')
colData$condition <- relevel(as.factor(colData$condition), 'CON')

cov_meta <- unique(meta[, c('Channel','AgeDeath','PMI','RIN','Sex','Race','Condition')])
rownames(cov_meta) <- paste0(cov_meta$Condition,'.',cov_meta$Channel)
cov_meta <- cov_meta[ , !(names(cov_meta) %in% c('Channel','Condition'))]
cov_meta$AgeDeath <- as.numeric(cov_meta$AgeDeath)
cov_meta$PMI <- as.numeric(cov_meta$PMI)
cov_meta$RIN <- as.numeric(cov_meta$RIN)
cov_meta$Sex <- relevel(as.factor(cov_meta$Sex), 'M')
cov_meta$Race <- relevel(as.factor(cov_meta$Race), 'W')

colData <- merge(colData, cov_meta, by = 'row.names', all.y=FALSE)
rownames(colData) <- colData$Row.names
colData <- colData[ , !(names(colData) %in% c('Row.names'))]

dds <- DESeqDataSetFromMatrix(countData = counts_cellType_filtered,
                              colData = colData,
                              design = ~AgeDeath + PMI + RIN + Sex + Race + condition)

# run DESeq2
dds <- DESeq(dds, 
             test = "LRT",
             reduced = ~AgeDeath + PMI + RIN + Sex + Race)
res <- results(dds, name = glue("condition_{condition}_vs_CON"))
sorted_res <- as.data.frame(res[order(res$pvalue, decreasing = FALSE), ])
sorted_res_sig <- sorted_res[(!is.na(sorted_res$padj)) & (sorted_res$padj < 0.05), ]
sorted_res_sig <- sorted_res_sig[(sorted_res_sig$log2FoldChange > 0.263) | (sorted_res_sig$log2FoldChange < -0.263), ]
sorted_res_sig$geneName <- row.names(sorted_res_sig)
write.table(sorted_res, glue('./DEG/{condition}-{celltype}_DEG.csv'),sep=',')
