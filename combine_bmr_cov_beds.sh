
while read p; do

cat $p.covered.bps.bed | awk -v var="$p" '{print $1"\t"$2"\t"$3"\t"$4"\t"var}' >> master_bmr.bed
done<$1
