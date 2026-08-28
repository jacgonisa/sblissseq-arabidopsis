#!/usr/bin/env bash
# DSBs per Mb of UNIQUELY MAPPED read sequence (coverage-normalised by aligned
# bases, not by region length). Denominator = samtools bedcov on the HISAT2
# deduplicated (primary-only) BAM = total aligned bases ("mapped read Mb") in the region.
# Numerator = DISTINCT DSB sites (unique 5' break positions). NB: total break
# *events* == #reads (1 break/read in BLISS), so events/bases is degenerate
# (= 1/read-length); distinct sites is the informative numerator.
# NOR 45S rDNA is its own category (its pile-up would otherwise swamp 'arms').
set -euo pipefail

BASE=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
REG=$BASE/results/analysis/regions
OUT=$BASE/results/analysis/dsb_per_mapped_mb.csv

echo "sample,region,dsb_sites,dsb_events,mapped_bases,sites_per_mapped_mb,events_per_mapped_mb" > "$OUT"
for s in BA1 BA2 old_BA1_BA2; do
  bam=$BASE/results/bam/${s}.hisat2.dedup.bam
  bg=$BASE/results/breaks/${s}.break.bedgraph
  for reg in arms_noNOR pericentromere centromere nor; do
    rbed=$REG/${reg}.bed
    bases=$(samtools bedcov "$rbed" "$bam" | awk '{s+=$NF} END{print s+0}')
    read sites events < <(bedtools intersect -a "$bg" -b "$rbed" -u | awk '{n++; e+=$4} END{print n+0, e+0}')
    sp=$(awk -v b="$sites"  -v m="$bases" 'BEGIN{ printf (m>0)? "%.1f":"NA", b/(m/1e6) }')
    ep=$(awk -v b="$events" -v m="$bases" 'BEGIN{ printf (m>0)? "%.1f":"NA", b/(m/1e6) }')
    echo "${s},${reg},${sites},${events},${bases},${sp},${ep}" >> "$OUT"
    printf "%-12s %-15s sites=%-9s events=%-10s mappedMb=%-9.1f sites/Mb=%-8s events/Mb=%s\n" \
           "$s" "$reg" "$sites" "$events" "$(awk -v m=$bases 'BEGIN{print m/1e6}')" "$sp" "$ep"
  done
done
echo ""
echo "Wrote $OUT"
