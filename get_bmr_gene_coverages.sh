while read p; do
bedtools genomecov -ibam /Volumes/graid/mTCGA_mm10_WES_bams/$p.bam -bg -g mm10.chromSizes.txt | awk '$4 > 7' | bedtools merge -i stdin | bedtools intersect -a stdin -b mm10_exons_plus_100.bed -wo | cut -f 1-7 | uniq | awk '{if($2<=$5 && $3<=$6) print $1"\t"$5"\t"$3"\t"$7; else if($2>=$5 && $3<=$6) print $1"\t"$2"\t"$3"\t"$7; else if($2>=$5 && $3>=$6) print $1"\t"$2"\t"$6"\t"$7; else if($2<$5 && $3>$6) print $1"\t"$5"\t"$6"\t"$7}' | bedtools intersect -a stdin -b mm10_bmr_ref_intervals.bed -wo | awk '{print $1"\t"$2"\t"$3"\t"$8}' > bmr_covered_bps/$p.covered.bps.bed 
done<$1
