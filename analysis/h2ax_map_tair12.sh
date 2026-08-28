#!/usr/bin/env bash
# Map WT gamma-H2AX IP + INPUT (PE ChIP) to TAIR12 (Col-CC) for correlation with sBLISS.
set -euo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
IDX=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/reference/bowtie2_index_TAIR12/Col-CC
RAW=/mnt/ssd-8tb/NOVOGENE/H2AX_June2026/X204SC26054750-Z01-F001/01.RawData
OUT=/mnt/ssd-8tb/NOVOGENE/H2AX_June2026/results_TAIR12
mkdir -p $OUT/bam $OUT/bigwig $OUT/logs
map_merge(){  # name  dir1 dir2
  local NAME=$1; shift; local bams=()
  for d in "$@"; do
    R1=$(ls $RAW/$d/*_1.fq.gz | paste -sd,); R2=$(ls $RAW/$d/*_2.fq.gz | paste -sd,)
    bowtie2 --very-sensitive-local -p 12 -x $IDX -1 $R1 -2 $R2 2>$OUT/logs/${d}.bt2.log \
      | samtools view -b -f 2 -F 0x900 - | samtools sort -@8 -o $OUT/bam/${d}.bam
    bams+=($OUT/bam/${d}.bam)
    echo "  $d: $(grep 'overall' $OUT/logs/${d}.bt2.log)"
  done
  samtools merge -f $OUT/bam/${NAME}.merged.bam "${bams[@]}"
  samtools sort -n -@8 $OUT/bam/${NAME}.merged.bam -o $OUT/bam/${NAME}.nsort.bam
  samtools fixmate -m $OUT/bam/${NAME}.nsort.bam - | samtools sort -@8 - | samtools markdup -r - $OUT/bam/${NAME}.dedup.bam
  samtools index $OUT/bam/${NAME}.dedup.bam
  rm -f $OUT/bam/${NAME}.merged.bam $OUT/bam/${NAME}.nsort.bam
  bamCoverage -b $OUT/bam/${NAME}.dedup.bam -o $OUT/bigwig/${NAME}.cpm.bw --normalizeUsing CPM --binSize 50 -p 8 2>/dev/null
  echo "$NAME dedup reads: $(samtools view -c -F4 $OUT/bam/${NAME}.dedup.bam)"
}
echo "=== gammaH2AX IP (WT) ==="; map_merge IP_WT   IP_WT_1 IP_WT_2
echo "=== INPUT (WT) ===";        map_merge INPUT_WT INPUT_WT_1 INPUT_WT_2
echo "H2AX_TAIR12_DONE"
