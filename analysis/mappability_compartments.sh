#!/usr/bin/env bash
# Mappability-aware compartment DSB density (TAIR12).
# Combines GenMap k=127 uniquely-mappable bp with all-read and MAPQ>=10 break tracks.
set -euo pipefail
BASE=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
R=$BASE/results_TAIR12/analysis/regions
UNIQ=$BASE/results_TAIR12/analysis/mappability/uniq_k127.bed
OUT=$BASE/results_TAIR12/analysis/mappability_compartments.csv

echo "sample,region,reg_bp,mappable_bp,mappable_pct,all_sites,all_sites_per_Mb_DNA,q10_sites,q10_events,q10_sites_per_mappable_Mb,q10_reads_per_site" > "$OUT"
for s in BA1 BA2 old_BA1_BA2; do
  all=$BASE/results_TAIR12/breaks/${s}.break.bedgraph
  q10=$BASE/results_TAIR12/breaks/${s}.q10.break.bedgraph
  for reg in arms pericentromere centromere nor_45s rdna_5s; do
    rbed=$R/${reg}.bed
    rb=$(awk '{s+=$3-$2}END{print s}' "$rbed")
    mb=$(bedtools intersect -a "$UNIQ" -b "$rbed" | awk '{s+=$3-$2}END{print s+0}')
    as=$(bedtools intersect -a "$all" -b "$rbed" -u | wc -l)
    read qs qe < <(bedtools intersect -a "$q10" -b "$rbed" -u | awk '{n++;e+=$4}END{print n+0,e+0}')
    awk -v s=$s -v r=$reg -v rb=$rb -v mb=$mb -v as=$as -v qs=$qs -v qe=$qe 'BEGIN{
      printf "%s,%s,%d,%d,%.2f,%d,%.0f,%d,%d,%.0f,%.2f\n",
        s,r,rb,mb,100*mb/rb,as,as/(rb/1e6),qs,qe,(mb>0?qs/(mb/1e6):0),(qs>0?qe/qs:0)}' >> "$OUT"
  done
done
echo "Wrote $OUT"; column -t -s, "$OUT"
