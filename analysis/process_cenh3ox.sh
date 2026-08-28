#!/usr/bin/env bash
set -eo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
BASE=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
RAW=$BASE/X204SC26053809-Z01-F001/01.RawData
OUT=$BASE/results_TAIR12
IDX=$BASE/reference/bowtie2_index_TAIR12/Col-CC
DEMUX=$BASE/pipeline/scripts/bliss_demux.py
MKBREAK=$BASE/pipeline/scripts/make_break_bed.R
SIZES=$OUT/analysis/regions/genome_chr15.sizes
mkdir -p $OUT/merged $OUT/umi $OUT/trimmed $OUT/bam $OUT/breaks $OUT/logs $OUT/analysis/bigwig_cpm $OUT/analysis/bigwig_raw
proc(){ local S=$1 BC=$2; shift 2; local LANES=("$@")
  echo "[$(date +%H:%M)] === $S (barcode $BC), $(echo ${LANES[@]}|wc -w) lane(s) ==="
  cat "${LANES[@]}" > $OUT/merged/${S}_R1.fq.gz
  python $DEMUX $OUT/merged/${S}_R1.fq.gz $OUT/umi/${S}_R1_demux.fq.gz --barcode $BC --max-mm 1 --stats $OUT/logs/${S}.demux.stats.txt
  trim_galore -q 20 --length 20 --output_dir $OUT/trimmed --gzip $OUT/umi/${S}_R1_demux.fq.gz 2> $OUT/logs/${S}.trimgalore.log
  nice -n19 ionice -c3 bowtie2 --very-sensitive -N 1 -p 4 -x $IDX -U $OUT/trimmed/${S}_R1_demux_trimmed.fq.gz 2> $OUT/logs/${S}.bowtie2.log \
    | samtools view -bS | nice -n19 samtools sort -@2 -o $OUT/bam/${S}.bowtie2.bam
  samtools index $OUT/bam/${S}.bowtie2.bam
  umi_tools dedup -I $OUT/bam/${S}.bowtie2.bam -S $OUT/bam/${S}.bowtie2.dedup.bam --umi-separator="+" --log=$OUT/logs/${S}.dedup.log
  samtools index $OUT/bam/${S}.bowtie2.dedup.bam
  bedtools bamtobed -i $OUT/bam/${S}.bowtie2.dedup.bam | grep -E '^Chr[1-5]	' > $OUT/breaks/${S}.flat.bed
  Rscript $MKBREAK $OUT/breaks/${S}.flat.bed $OUT/breaks/${S}.break.bed $OUT/breaks/${S}.break.bedgraph
  # raw + CPM bigwigs (break density)
  LC_ALL=C sort -k1,1 -k2,2n $OUT/breaks/${S}.break.bedgraph > /tmp/${S}.s.bg
  bedGraphToBigWig /tmp/${S}.s.bg $SIZES $OUT/analysis/bigwig_raw/${S}.raw.bw
  tot=$(awk '{s+=$4}END{print s}' /tmp/${S}.s.bg)
  awk -v t=$tot 'BEGIN{OFS="\t"}{print $1,$2,$3,$4*1e6/t}' /tmp/${S}.s.bg > /tmp/${S}.cpm.bg
  bedGraphToBigWig /tmp/${S}.cpm.bg $SIZES $OUT/analysis/bigwig_cpm/${S}.cpm.bw
  rm -f $OUT/breaks/${S}.flat.bed $OUT/merged/${S}_R1.fq.gz /tmp/${S}.s.bg /tmp/${S}.cpm.bg
  echo "[$(date +%H:%M)] DONE $S: dedup reads=$(samtools view -c $OUT/bam/${S}.bowtie2.dedup.bam) break events=$tot sites=$(wc -l <$OUT/breaks/${S}.break.bedgraph)"
}
proc BA1_WT      CATCACGC $RAW/BA1_WT/*_L8_1.fq.gz $RAW/BA1_WT/*_L7_1.fq.gz
proc BA2_OX_40   GTCGTCGC $RAW/BA2_OX_40/*_L6_1.fq.gz $RAW/BA2_OX_40/*_L5_1.fq.gz
proc BA2_OX_80_1 GTCGTCGC $RAW/BA2_OX_80_1/*_L7_1.fq.gz
proc BA2_OX_80_2 GTCGTCGC $RAW/BA2_OX_80_2/*_L6_1.fq.gz $RAW/BA2_OX_80_2/*_L5_1.fq.gz
echo "CENH3OX_BLISS_DONE"
