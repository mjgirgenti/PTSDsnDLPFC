#!/bin/bash

module load miniconda
conda activate ldsc
module load BEDTools
module load dSQ


### Credit: Flynn Chen ###


base_dir="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC"
download_dir="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/jing_download"
outputs_dir="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/unique_peaks/outputs"
results_dir="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/unique_peaks/results"
bed_dir="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/unique_peaks/peaks_lifted"

mkdir -p $results_dir

rm joblist_02.txt

for sumstats in $base_dir/14_snps_subset/*.sumstats.gz
do
    for baselines in $output_dir/*;
    do
        sumstat=${sumstats##*/}
        sumstat_base=${sumstat%.sumstats.gz}
        baseline=${baselines##*/}
    tissue=${baseline#*.}
    echo "cd $results_dir; python $base_dir/ldsc/ldsc.py --h2 $sumstats --w-ld-chr $base_dir/trial/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --ref-ld-chr ${baselines}/${baseline}.,/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/ldct_dir/output/replicate_baseline/baseline_v1.1. --overlap-annot --frqfile-chr $base_dir/trial/1000G_Phase3_frq/1000G.EUR.QC. --out ${sumstat_base}.${baseline} --print-coefficients" >> joblist_02.txt 
    done
done 

tasklist="$PWD/joblist_02.txt"
mkdir -p $PWD/log
cd $PWD/log
dSQ.py --jobfile $tasklist --submit  --constraint=avx512 --max-jobs 100 --partition=pi_gerstein -t 72:00:00 --mem=5G --job-name=02

