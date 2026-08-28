#!/usr/bin/env bash
# Remap an external PE ChIP fastq pair to OUR TAIR12 (Col-CC) index,
# identically to our γH2AX IP (--very-sensitive-local, properly-paired, dedup).
# Usage: map_external_tair12.sh <NAME> <R1> <R2>
set -euo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
IDX=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/reference/bowtie2_index_TAIR12/Col-CC
OUT=/mnt/ssd-8tb/NOVOGENE/H2AX_June2026/results_external_TAIR12
mkdir -p $OUT/bam $OUT/bigwig $OUT/logs
NAME=$1; R1=$2; R2=$3
echo "[$(date +%H:%M)] mapping $NAME -> TAIR12"
bowtie2 --very-sensitive-local -p 12 -x $IDX -1 "$R1" -2 "$R2" 2>$OUT/logs/$NAME.bt2.log \
  | samtools view -b -f 2 -F 0x900 - | samtools sort -@8 -o $OUT/bam/$NAME.psort.bam
samtools sort -n -@8 $OUT/bam/$NAME.psort.bam -o $OUT/bam/$NAME.nsort.bam
samtools fixmate -m $OUT/bam/$NAME.nsort.bam - | samtools sort -@8 - | samtools markdup -r - $OUT/bam/$NAME.dedup.bam
samtools index $OUT/bam/$NAME.dedup.bam
samtools flagstat $OUT/bam/$NAME.dedup.bam > $OUT/logs/$NAME.flagstat
rm -f $OUT/bam/$NAME.psort.bam $OUT/bam/$NAME.nsort.bam
bamCoverage -b $OUT/bam/$NAME.dedup.bam -o $OUT/bigwig/$NAME.all.cpm.bw --normalizeUsing CPM --binSize 50 -p 8 2>/dev/null
bamCoverage -b $OUT/bam/$NAME.dedup.bam -o $OUT/bigwig/$NAME.q10.cpm.bw --minMappingQuality 10 --normalizeUsing CPM --binSize 50 -p 8 2>/dev/null
echo "[$(date +%H:%M)] DONE $NAME  ($(grep 'overall' $OUT/logs/$NAME.bt2.log))"
