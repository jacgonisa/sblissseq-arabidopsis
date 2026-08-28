#!/usr/bin/env bash
# Map deduplicated BLISS reads to the CEN180 2x178bp dimer (minimap2 -ax sr),
# count per-base coverage + 5' DSB positions on the 356bp dimer (spo11 SE counter).
set -euo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
MM2=$HOME/minimap2/minimap2
SPO11=/mnt/ssd-8tb/spo11-mnaseq-project
DIMER=$SPO11/centromere/index/CEN180_2x178bp.fa
COUNT=$SPO11/scripts/15_cen180_count_se.py
P=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
OUT=$P/results_TAIR12/analysis/cen180_dimer
mkdir -p "$OUT"

for s in BA1 BA2 old_BA1_BA2; do
  bam=$P/results_TAIR12/bam/${s}.bowtie2.dedup.bam
  echo "### $s -> dimer"
  samtools fastq -n "$bam" 2>/dev/null \
    | "$MM2" -ax sr -t 16 "$DIMER" - 2>"$OUT/${s}_minimap2.log" \
    | samtools view - \
    | python3 "$COUNT" "$OUT" "$s"
  mapped=$(awk '{s+=$2}END{print s}' "$OUT/${s}_cen180_se_r1_5prime.tsv")
  echo "  $s: total 5' DSB ends mapped to dimer = $mapped"
done
echo "DIMER_DONE"
