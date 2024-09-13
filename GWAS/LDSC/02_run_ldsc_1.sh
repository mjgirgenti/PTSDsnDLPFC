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

mkdir -p $outputs_dir
rm joblist_01.txt

for f in $download_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.*.bim;
do
	for tissue_beds in $bed_dir/*.hg19.lifted.bed;
	do
		filename=${f##*/}
		basename=${filename%.*}
		bfile=${f%.*}
		chr="${basename##*.}"
		tissue_bed=${tissue_beds##*/}
		tissue_name=${tissue_bed%.hg19.lifted.bed}
		echo "cd $outputs_dir; mkdir -p $tissue_name; python $base_dir/annotate_ldscore.py --annot-bed $tissue_beds --bfile $bfile --annot-name $tissue_name --out $tissue_name/$tissue_name.${chr}; python $base_dir/ldsc/ldsc.py --l2 --bfile $download_dir/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} --ld-wind-cm 1 --annot $tissue_name/$tissue_name.${chr}.annot.gz --out $tissue_name/$tissue_name.${chr} --print-snps $download_dir/hapmap3_snps/hm.${chr}.snp" >> joblist_01.txt
	done
done 


tasklist="$PWD/joblist_01.txt"
mkdir -p $PWD/log
cd $PWD/log
dSQ.py --jobfile $tasklist --submit  --constraint=avx512 --max-jobs 300 --partition=day -t 1:00:00 --mem=10G --job-name=01 

