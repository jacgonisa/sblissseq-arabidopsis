#!/usr/bin/env bash
# Generate BREAK bed/bedgraph for the TAIR12 (Col-CC) bowtie2 dedup BAMs.
# Same procedure as the ragtag run (bamtobed -> make_break_bed.R), Chr1-5 only.
set -euo pipefail
BASE=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
OUT=$BASE/results_TAIR12/breaks
mkdir -p "$OUT"
for s in BA1 BA2 old_BA1_BA2; do
  bam=$BASE/results_TAIR12/bam/${s}.bowtie2.dedup.bam
  echo "### $s : $bam"
  bedtools bamtobed -i "$bam" | grep -E '^Chr[1-5]	' > "$OUT/${s}.flat.bed"
  Rscript "$BASE/pipeline/scripts/make_break_bed.R" \
      "$OUT/${s}.flat.bed" "$OUT/${s}.break.bed" "$OUT/${s}.break.bedgraph"
  rm -f "$OUT/${s}.flat.bed"
  echo "  break events: $(awk '{s+=$5}END{print s}' "$OUT/${s}.break.bed")  unique sites: $(wc -l < "$OUT/${s}.break.bedgraph")"
done
echo "DONE breaks TAIR12"
