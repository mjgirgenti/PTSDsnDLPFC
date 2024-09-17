# LDSC

Use LDSC to assess heritability of various psychiatric/non-psychiatric traits using cell type-specific ATAC peaks. Please use the following steps:

- Format raw sumstats files: `sumstats.sh`
- Liftover peaks to hg19: `01_liftover.sh`
- Create annotation files: `02_run_ldsc_1.sh`
- Perform LD score regression: `03_run_ldsc_2.sh`
- Plot results as dotplot: `plot_ldsc_dotplot.R`

 
