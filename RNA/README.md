# snRNA-seq Analysis Pipeline

Here is the tutorial for snRNA-seq analysis.

## Preprocessing

We first run cellranger and cellbender to generate the initial raw count matrices. To run the pipeline, use: 

```
./preprocess.sh <sample_id> <fastq_location> <ref_genome_location>
```

## Quality Control

Then, we conduct initial quality control and doublet removal:

```
python filter_cells_genes_doublets.py
```

Please change the file path, output path, and sample ID accordingly in the file before running it. It should give the filtered count matrices as a `.zarr` file which can be loaded in pegasus.

## Clustering

Perform the first round of clustering using:
```
python ahyeon_snRNA_pipeline.py
```
Clean the UMAP and perform clustering again using:
```
python ahyeon_snRNA_pipeline_second.py
```

## Differentially Expressed Genes

Conduct DEG analysis using MAST, Wilcox, and DESeq2:  

```
Rscript MAST_DEG.R
Rscript WILCOX_DEG.R
Rscript DESEQ2_DEG.R
```

For all three scripts, please specify the cell type and condition. 

## Cell Aggregation

We also conduct meta-cell construction and WGCNA analysis on the metacells:

```
Rscript create_metacells.R
Rscript wgcna_metacells.R
```

Please specify the dataset name and path (the h5seurat and the corresponding metadata csv files) in the script. 
