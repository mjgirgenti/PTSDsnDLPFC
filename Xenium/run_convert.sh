#!/bin/bash
#SBATCH --job-name=ome
#SBATCH --out=omeconvert_%j.log
#SBATCH --time=1-
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --partition=pi_gerstein

module load miniconda
conda activate spatial
python /home/ah2428/project/xenpred/scripts/omeconvert.py /home/ah2428/girgenti/ah2428/xenium/data/5173_resegment/outs/morphology_focus_copy.ome.tif

