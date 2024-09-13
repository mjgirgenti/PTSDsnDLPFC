sortBed -i End.peaks.bed > End.peaks.sort.bed
sortBed -i Exc.peaks.bed > Exc.peaks.sort.bed
sortBed -i Inh.peaks.bed > Inh.peaks.sort.bed
sortBed -i Mic.peaks.bed > Mic.peaks.sort.bed
sortBed -i Oli.peaks.bed > Oli.peaks.sort.bed
sortBed -i OPC.peaks.bed > Opc.peaks.sort.bed

thresh=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
for t in ${thresh[@]};do
    echo $t
    echo End
    bedtools intersect -a End.peaks.sort.bed -b bulk.bed -f $t -sorted -u > out/End.${t}.peaks.bed
    echo Exc 
    bedtools intersect -a Exc.peaks.sort.bed -b bulk.bed -f $t -sorted -u > out/Exc.${t}.peaks.bed
    echo Inh
    bedtools intersect -a Inh.peaks.sort.bed -b bulk.bed -f $t -sorted -u > out/Inh.${t}.peaks.bed
    echo Mic
    bedtools intersect -a Mic.peaks.sort.bed -b bulk.bed -f $t -sorted -u > out/Mic.${t}.peaks.bed
    echo Oli
    bedtools intersect -a Oli.peaks.sort.bed -b bulk.bed -f $t -sorted -u > out/Oli.${t}.peaks.bed
    echo Opc
    bedtools intersect -a Opc.peaks.sort.bed -b bulk.bed -f $t -sorted -u > out/Opc.${t}.peaks.bed
done

