# Multiome Analysis Pipeline

This is the instruction and tutorial for the Multiome (RNA+ATAC) part of the analysis.

For the Multiome analysis, we followed the same protocol as shown in the psychENCODE analysis published [here](https://www.science.org/doi/10.1126/science.adi5199). You may refer to the tutorial from the consortium [here] (https://github.com/gersteinlab/PsychENCODE_SingleCell_Integrative).

Below, we will list the major steps and associating script files in this directory. Please refer to the comments in the code as they should be self-explanatory. 

## Preprocessing and initial QC

- Create arrow files from fragments: `create_archr.R`
- Integration RNA with ATAC: `rna_atac_integration.R`

## Peak Calling and CRE set creation

- Peak calling: `callpeaks_archr.R`
- Cell clustering: `multiome_clustering.R`
- Plot marker genes: `multiome_marker_genes.R`

## Peak-to-Gene linkage

- Create p2g linkage: `multiome_p2g.R`
- Find cell type-specific marker peaks: `get_celltype_marker_peaks.R`
- Create TF-CRE-Gene networks: `TF_Peak_Gene.R`
