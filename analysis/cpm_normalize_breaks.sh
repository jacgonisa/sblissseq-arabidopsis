#!/usr/bin/env bash
# CPM-normalise BREAK bedgraph files and convert to bigWig.
# sBLISS signal = break counts per position. CPM = count * 1e6 / total_breaks
# allows fair cross-sample comparison despite different library sizes.
set -euo pipefail

BREAKS=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/breaks
OUT=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/bigwig_cpm
CHRSIZES=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/reference/Col-0.chrom.sizes

mkdir -p "$OUT"

for sample in BA1 BA2 old_BA1_BA2; do
    bdg="$BREAKS/${sample}.break.bedgraph"   # HISAT2 both-strand counts
    # total breaks (sum of column 4)
    total=$(awk '{s+=$4} END{print s}' "$bdg")
    echo "$sample: total breaks = $total"

    # CPM-scale, keep only Chr1-5, sort, convert
    awk -v t="$total" 'BEGIN{OFS="\t"} $1 ~ /^Chr[1-5]$/ {printf "%s\t%d\t%d\t%.4f\n",$1,$2,$3,$4*1e6/t}' "$bdg" \
        | LC_ALL=C sort -k1,1 -k2,2n > "$OUT/${sample}.cpm.bedgraph"

    bedGraphToBigWig "$OUT/${sample}.cpm.bedgraph" \
        <(grep -P "^Chr[1-5]\t" "$CHRSIZES") \
        "$OUT/${sample}.cpm.bw" 2>/dev/null \
      || bedGraphToBigWig "$OUT/${sample}.cpm.bedgraph" /tmp/chr15_sizes.txt "$OUT/${sample}.cpm.bw"

    echo "  → $OUT/${sample}.cpm.bw"
done
echo "Done."
