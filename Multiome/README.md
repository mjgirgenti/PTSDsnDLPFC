# Multiome Analysis Pipeline

This is the instruction and tutorial for the Multiome (RNA+ATAC) part of the analysis.

For the Multiome analysis, we followed the same protocol as shown in the psychENCODE analysis published [here](https://www.science.org/doi/10.1126/science.adi5199). You may refer to the tutorial from the consortium [here] (https://github.com/gersteinlab/PsychENCODE_SingleCell_Integrative).

Below, we list the major steps and associated files in this directory. 

## Preprocessing and initial QC
- Clustering with Signac: `multiome_clustering.R`
- Plot marker genes: `multiome_marker_genes.R`
- Create arrow files from fragments: `create_archr_proj.R`

## Peak Calling and CRE set creation
- Peak calling on 7 celltypes: `callpeaks_archr_celltypes.R`
- Peak calling on 14 subtypes: `callpeaks_archr_subtypes.R`

## Peak-to-Gene linkage
- Create p2g linkage: `multiome_p2g.R`
- Cell type-specific marker peaks: `get_celltype_marker_peaks.R`

