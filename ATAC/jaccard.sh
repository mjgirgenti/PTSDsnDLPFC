module load BEDTools/2.30.0-GCCcore-10.2.0

peaksdir=/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/peaks

cd $peaksdir

celltypes=('EXC' 'OLI' 'INH' 'OPC' 'MIC' 'END' 'AST')


for i in ${celltypes[*]}
do
    for j in ${celltypes[*]}
    do
        echo $i $j
        bedtools jaccard -a $i.peaks.bed -b $j.peaks.bed | sed -n '2 p'  >> jaccard_outs.txt
    done
done
