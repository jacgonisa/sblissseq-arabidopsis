#!/usr/bin/env bash
set -eo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
MM2=$HOME/minimap2/minimap2
MMI=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/reference/TAIR12.sr.mmi
OUT=/mnt/ssd-8tb/NOVOGENE/H2AX_June2026/results_external_TAIR12
FQ=/mnt/ssd-8tb/fastq_files/G4s
mkdir -p $OUT/bam $OUT/bigwig $OUT/logs
map(){ local N=$1 R1=$2 R2=$3
  echo "[$(date +%H:%M)] $N (no -f2; keep all mapped primary)"
  $MM2 -ax sr -t 24 $MMI "$R1" "$R2" 2>$OUT/logs/$N.mm2.log \
    | samtools view -b -F 0x904 - | samtools sort -@8 -o $OUT/bam/$N.psort.bam
  samtools sort -n -@8 $OUT/bam/$N.psort.bam -o $OUT/bam/$N.nsort.bam
  samtools fixmate -m $OUT/bam/$N.nsort.bam - 2>/dev/null | samtools sort -@8 - | samtools markdup -r - $OUT/bam/$N.q0.dedup.bam
  samtools index $OUT/bam/$N.q0.dedup.bam; rm -f $OUT/bam/$N.psort.bam $OUT/bam/$N.nsort.bam
  bamCoverage -b $OUT/bam/$N.q0.dedup.bam -o $OUT/bigwig/$N.q0.cpm.bw --normalizeUsing CPM --binSize 50 -p 8 2>/dev/null
  echo "[$(date +%H:%M)] DONE $N : $(samtools view -c -F0x904 $OUT/bam/$N.q0.dedup.bam) reads"
}
map Li_K_rep1_q0 $FQ/Li_K_rep1/SRR6724283_1.fastq.gz $FQ/Li_K_rep1/SRR6724283_2.fastq.gz
map Li_K_rep2_q0 $FQ/Li_K_rep2/SRR6724284_1.fastq.gz $FQ/Li_K_rep2/SRR6724284_2.fastq.gz
map Li_KPDS_q0   $FQ/Li_KPDS/SRR6724285_1.fastq.gz   $FQ/Li_KPDS/SRR6724285_2.fastq.gz
echo "G4_FIXED_DONE"
