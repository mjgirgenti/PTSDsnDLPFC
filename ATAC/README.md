# snATAC-seq Analysis Pipeline

This is the tutorial for snATAC-seq analysis.

For the snATAC-seq analysis, we followed the same protocol as shown in the psychENCODE analysis published [here](https://www.science.org/doi/10.1126/science.adi5199). You may refer to the tutorial from the consortium [here] (https://github.com/gersteinlab/PsychENCODE_SingleCell_Integrative).

Below, we list the major steps and associated scripts in this directory. 

## Preprocessing and initial QC

- Create arrow files from fragments: `create_arrow.R`
- Create bigwig files: `create_bigwig.R`
- Cell clustering: `clustering.R`
- Plot marker genes: `marker_genes_heatmap.R`

## Peak Calling 

- Celltype peaks: `callpeaks.R`
- Celltype-condition peaks: `callpeaks_cluster_cond.R`

## RNA ATAC Integration

- Integrate using high-quality RNA samples: `rna_atac_integration.R`

## Peak-to-Gene Links

- Add peak-to-gene links: `peak2genelinks.R`

## Transcription factor motif analysis

- Add TF deviations matrix: `addDeviationsMatrix.R`
- Get ranked TFs per celltype: `celltype_motifs_rank.R`
- TF Footprinting: `TF_umap_footprint.R`
- TF-CRE-Gene links: `TF_CRE_Gene.R`

## Connect to downstream analysis

- Convert ArchR project into Signac project: `archrtosignac.R`
- Subset cells based on QC metrics: `subset_cells.R`
