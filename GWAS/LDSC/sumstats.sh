#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH -p day
#SBATCH --job-name=sumstat
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH --cpus-per-task=1

module load miniconda
conda activate ldsc

folder="/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/GWAS_PGC2"
cd $folder
for snp in $folder/*.txt
do
    n=$(sed -n '2p' | awk '{print $4}')
    /gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/LDSC/ldsc/munge_sumstats.py --sumstats $snp --N $n --snpid SNP --merge-alleles w_hm3.snplist 
done 
