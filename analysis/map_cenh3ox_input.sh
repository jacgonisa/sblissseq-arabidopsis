#!/usr/bin/env bash
set -eo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
BASE=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
IDX=$BASE/reference/bowtie2_index_TAIR12/Col-CC
OUT=$BASE/results_TAIR12; FQ=$BASE/input_cenh3ox
mkdir -p $OUT/bam $OUT/analysis/coverage_cpm $OUT/logs
N=cenh3ox_input
nice -n19 ionice -c3 bowtie2 --very-sensitive-local -p 4 -x $IDX -1 $FQ/${N}_R1.fastq.gz -2 $FQ/${N}_R2.fastq.gz 2>$OUT/logs/$N.bowtie2.log \
  | samtools view -b -f 2 -F 0x900 - | nice -n19 samtools sort -@2 -o $OUT/bam/$N.psort.bam
samtools sort -n -@2 $OUT/bam/$N.psort.bam -o $OUT/bam/$N.nsort.bam
samtools fixmate -m $OUT/bam/$N.nsort.bam - | samtools sort -@2 - | samtools markdup -r - $OUT/bam/$N.dedup.bam
samtools index $OUT/bam/$N.dedup.bam; rm -f $OUT/bam/$N.psort.bam $OUT/bam/$N.nsort.bam
bamCoverage -b $OUT/bam/$N.dedup.bam -o $OUT/analysis/coverage_cpm/$N.cpm.bw --normalizeUsing CPM --binSize 50 -p 4 2>/dev/null
echo "[$(date +%H:%M)] DONE $N: $(samtools view -c -F0x904 $OUT/bam/$N.dedup.bam) reads; $(grep overall $OUT/logs/$N.bowtie2.log)"
