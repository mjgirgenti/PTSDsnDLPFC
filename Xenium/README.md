# Xenium Analysis Pipeline

Single-cell resolution spatial transcriptomics Xenium was used to validate our RNA and CCC results. Please see the following steps for Xenium analysis:

- Resegment for nuclei boundaries: `resegment_jobs.txt`
- Single sample QC: `single_sample_qc.py`
- Clustering: `cluster_combined.py`
- Plot spatially variable genes: `plot_svg.py`
- DEG analysis: `../RNA/MAST_DEG.R`
- Plot snXenium DLPFC slide: `plot_snxenium_slide.py`
- Plot scXenium DLPFC slide: `plot_scxenium_slide.py`
