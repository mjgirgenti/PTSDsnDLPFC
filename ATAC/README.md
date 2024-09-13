# snATAC-seq Analysis Pipeline

This is the instruction and tutorial for the snATAC-seq part of the analysis.

For the scATAC-seq analysis, we followed the same protocol as shown in the psychENCODE analysis published [here](https://www.science.org/doi/10.1126/science.adi5199). You may refer to the tutorial from the consortium [here] (https://github.com/gersteinlab/PsychENCODE_SingleCell_Integrative).

Below, we will list the major steps and associating script files in this directory. Please refer to the comments in the code as they should be self-explanatory. 

## Preprocessing and initial QC

- Create arrow files from fragments: `create_arrow.R`
- Create bigwig files: `create_bigwig.R`
- Subset cells based on QC metrics: `subset_cells.R`

## Peak Calling and CRE set creation

- Peak calling: `callpeaks.R`
- Create cell by peak matrix using the generated/added peak set: `addPeakMatrix.R`
- Cell clustering: `clustering.R`
- Plot marker genes: `marker_genes_heatmap.R`

## Transcription factor and enrichment analysis

- Differentially accessed peaks: `diffpeaks_tf.R`
- Add group coverage for downstream analysis: `add_groupcoverages.R`
- Add deviation matrices: `addDeviationsMatrix.R`
- Add ChromVAR TF enrichment analysis: `chromvar.R`
- TF Footprinting: `footprinting,R`

## Connect to downstream analysis

- Convert ArchR project into Signac project: `archrtosignac.R`
