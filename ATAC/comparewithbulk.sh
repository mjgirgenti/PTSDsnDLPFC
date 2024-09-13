#!/bin/bash

#SBATCH --job-name=bulk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G

BULK1='/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/PEC_bulk_ATAC_individual_merge_sort_cat.sorted.bed'

#path='/home/ah2428/girgenti_scratch60/atac/subclass_proj/PeakBeds_Subclass'
path='/gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/final_proj_copy2/PeakBeds'

clusters=()
for file in $path/*
do
  cluster=$(echo $file | cut -d'/' -f12 | cut -d'.' -f1 )
  echo $cluster
  bedtools window -a $file -b $BULK1 -w 500 -v | uniq > /gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/final_proj_copy2/CompareWithBulk_jason/${cluster}_nonoverlap_bulk_window500.bed
  bedtools window -a $file -b $BULK1 -w 1000 -v | uniq > /gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/final_proj_copy2/CompareWithBulk_jason/${cluster}_nonoverlap_bulk_window1000.bed
  
  total=$(wc -l < $file)
  filetotal1_500=$(wc -l < /gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/final_proj_copy2/CompareWithBulk_jason/${cluster}_nonoverlap_bulk_window500.bed)
  filetotal1_1000=$(wc -l < /gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/ATAC/final_proj_copy2/CompareWithBulk_jason/${cluster}_nonoverlap_bulk_window1000.bed)
  echo "scale=5 ; $filetotal1_500 / $total" | bc
  echo "scale=5 ; $filetotal1_1000 / $total" | bc
  echo "scale=5 ; $filetotal1_500 / $total" | bc >> bulk_window500_percents.txt
  echo "scale=5 ; $filetotal1_1000 / $total" | bc >> bulk_window1000_percents.txt
done

