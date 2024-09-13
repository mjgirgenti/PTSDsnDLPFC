#!/bin/bash
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -p pi_gerstein
#SBATCH --job-name=hg19
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu 10G
#SBATCH --cpus-per-task=1


### Credit: Flynn Chen ###


raw_input="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/unique_peaks/peaks"
processed_output="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/unique_peaks/peaks_lifted"

for file in $raw_input/*
do
        file_name="${file##*/}"
        file_basename="${file_name%.bed}"
        echo "lifting over $file to ${file_basename}.temp.bed to ${processed_output}/${file_basename}.hg19.lifted.bed"
        awk '{print $1"\t"$2"\t"$3"\t"$4}' $file | sort -k 1,1 -k2,2n > $processed_output/${file_basename}.temp.bed
        /gpfs/gibbs/pi/girgenti/JZhang/commonData/PTSD/LDSC/liftOver $processed_output/${file_basename}.temp.bed /gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/hg38ToHg19.over.chain ${processed_output}/${file_basename}.hg19.lifted.bed ${processed_output}/${file_basename}.peaks.bed
done
