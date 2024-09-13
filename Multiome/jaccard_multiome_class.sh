module load BEDTools/2.30.0-GCCcore-10.2.0

multiome_peaksdir=/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/Multiome/25_samples_archr_class_peaks/PeakBeds


celltypes=('EXC' 'INH' 'OLI' 'OPC' 'END' 'AST' 'MIC')


for i in ${celltypes[*]}
do
    for j in ${celltypes[*]}
    do
        echo $i $j
        bedtools jaccard -a $multiome_peaksdir/$i.peaks.bed -b $multiome_peaksdir/$j.peaks.bed | sed -n '2 p'  >> jaccard_multiome_class_outs.txt
    done
done
