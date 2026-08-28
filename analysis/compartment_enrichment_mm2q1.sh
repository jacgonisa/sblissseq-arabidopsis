#!/usr/bin/env bash
# Compartment DSB enrichment from minimap2 MAPQ>=1 (uniquely-mapped) BLISS breaks.
# Fold = (region DSB events per Mb of DNA) / (genome-wide DSB events per Mb).
set -euo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
P=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
R=$P/results_TAIR12/analysis/regions
OUT=$P/results_TAIR12/analysis/compartment_enrichment_mm2q1.csv
GMB=$(awk '$1~/^Chr[1-5]$/{s+=$2}END{print s/1e6}' $R/genome_chr15.sizes)   # genome Mb (Chr1-5)

echo "sample,region,events,sites,region_Mb,events_per_Mb,sites_per_Mb,fold_events,fold_sites" > "$OUT"
for s in BA1 BA2 old_BA1_BA2; do
  bg=$P/results_TAIR12/breaks/${s}.mm2q1.break.bedgraph
  totE=$(awk '$1~/^Chr[1-5]$/{e+=$4}END{print e}' "$bg"); totS=$(awk '$1~/^Chr[1-5]$/' "$bg" | wc -l)
  meanE=$(awk -v t=$totE -v g=$GMB 'BEGIN{print t/g}')   # genome events/Mb
  meanS=$(awk -v t=$totS -v g=$GMB 'BEGIN{print t/g}')
  for reg in arms pericentromere centromere nor_45s rdna_5s; do
    rb=$(awk '{s+=$3-$2}END{print s}' $R/${reg}.bed); rMb=$(awk -v b=$rb 'BEGIN{print b/1e6}')
    read e n < <(bedtools intersect -a "$bg" -b $R/${reg}.bed -u 2>/dev/null | awk '{e+=$4;n++}END{print e+0,n+0}')
    awk -v s=$s -v r=$reg -v e=$e -v n=$n -v rMb=$rMb -v meanE=$meanE -v meanS=$meanS 'BEGIN{
      epm=e/rMb; spm=n/rMb;
      printf "%s,%s,%d,%d,%.2f,%.0f,%.0f,%.2f,%.2f\n",s,r,e,n,rMb,epm,spm,epm/meanE,spm/meanS}' >> "$OUT"
  done
done
echo "=== minimap2 MAPQ>=1 compartment enrichment (fold over genome-mean DSB/Mb) ==="
column -t -s, "$OUT"
echo ""
echo "=== ranking by fold_events (BA1) ==="
awk -F, 'NR>1 && $1=="BA1"{print $8"x  "$2}' "$OUT" | sort -rn