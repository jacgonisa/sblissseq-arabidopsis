#!/usr/bin/env bash
set -eo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
MM2=$HOME/minimap2/minimap2
GENOME=/mnt/ssd-4tb/crisanto_project/genome/TAIR12/GCA_028009825.2_Col-CC_genomic_withorganelles.fna
MMI=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/reference/TAIR12.sr.mmi
OUT=/mnt/ssd-8tb/NOVOGENE/H2AX_June2026/results_external_TAIR12
FQ=/mnt/ssd-8tb/fastq_files
mkdir -p $OUT/bam $OUT/bigwig $OUT/logs
[ -f $MMI ] || { echo "[idx] building sr index"; $MM2 -x sr -d $MMI $GENOME; }
map(){  # NAME R1 R2
  local N=$1 R1=$2 R2=$3
  echo "[$(date +%H:%M)] $N : minimap2 -ax sr (MAPQ>=0, multireads kept)"
  $MM2 -ax sr -t 24 $MMI "$R1" "$R2" 2>$OUT/logs/$N.mm2.log \
    | samtools view -b -f 2 -F 0x900 - | samtools sort -@8 -o $OUT/bam/$N.psort.bam
  samtools sort -n -@8 $OUT/bam/$N.psort.bam -o $OUT/bam/$N.nsort.bam
  samtools fixmate -m $OUT/bam/$N.nsort.bam - | samtools sort -@8 - | samtools markdup -r - $OUT/bam/$N.q0.dedup.bam
  samtools index $OUT/bam/$N.q0.dedup.bam
  samtools flagstat $OUT/bam/$N.q0.dedup.bam > $OUT/logs/$N.q0.flagstat
  rm -f $OUT/bam/$N.psort.bam $OUT/bam/$N.nsort.bam
  bamCoverage -b $OUT/bam/$N.q0.dedup.bam -o $OUT/bigwig/$N.q0.cpm.bw --normalizeUsing CPM --binSize 50 -p 8 2>/dev/null
  echo "[$(date +%H:%M)] DONE $N"
}
map WT_RAD51_q0 $FQ/RAD51/WT/RAD51/SRR24938940_1.fastq.gz   $FQ/RAD51/WT/RAD51/SRR24938940_2.fastq.gz
map WT_RPA1A_q0 $FQ/RAD51/WT/RPA1A/SRR24938935_1.fastq.gz   $FQ/RAD51/WT/RPA1A/SRR24938935_2.fastq.gz
map WT_INPUT_q0 $FQ/RAD51/WT/INPUT/SRR24938942_1.fastq.gz   $FQ/RAD51/WT/INPUT/SRR24938942_2.fastq.gz
map ssDRIP_q0   $FQ/ssDRIP/wt-ssDRIP_SRR8097477_1.fastq.gz.fastp $FQ/ssDRIP/wt-ssDRIP_SRR8097477_2.fastq.gz.fastp
echo "ALL_REPAIR_RLOOP_DONE"
