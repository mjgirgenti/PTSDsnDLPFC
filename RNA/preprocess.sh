#!/bin/bash

# Conduct sample preprocessing and background removal using cellranger and cellbender

# Sample ID
sample=$1

# FASTQ directory
fastq_path=$2

# Path to reference genome
ref_genome=$3

# Run cellranger
cellranger count --id=${sample} \
                 --transcriptome=${ref_genome} \
                 --fastqs=${fastq_path} \
                 --sample=${sample}_HHT,${sample}_HNT \
                 --include-introns \
                 --chemistry=auto \
                 --localcores=25 --localmem=100

# Run cellbender for background removal
cells=$(tail -1 ${sample}/outs/metrics_summary.csv | cut -d'"' -f2 | sed 's/,//g')
cells=$(echo $((cells+0)))
cellbender remove-background \
        --input ${sample}/outs/raw_feature_bc_matrix.h5 \
        --output ${sample}/${sample}_cellbender.h5 \
        --cuda \
        --expected-cells $cells \
        --total-droplets-included $(($cells+20000)) \
        --fpr=0.01 \
        --epochs 150