#!/usr/bin/env bash
# BLISS DSB enrichment vs WGA gDNA, per compartment (minimap2, multimappers kept).
# enrichment = (BLISS reads in region / BLISS total) / (gDNA reads in region / gDNA total)
# gDNA carries the same mappability + copy-number (rDNA collapse) bias, so >1 = DSB-specific enrichment.
set -euo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
P=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026
R=$P/results_TAIR12/analysis/regions
GBAM=${BLISS_GDNABAM:-$P/results_TAIR12/bam/wga_genome.bam}        # gDNA BAM
BSUF=${BLISS_BAMSUFFIX:-minimap2.dedup}                            # BLISS bam suffix
TAG=${BLISS_TAG:-minimap2}
OUT=$P/results_TAIR12/analysis/bliss_vs_gdna_enrichment_${TAG}.csv
regions="arms pericentromere centromere nor_45s rdna_5s"
echo "gDNA=$GBAM  BLISS=*.${BSUF}.bam  tag=$TAG"

# gDNA totals + per-region (Chr1-5)
gTot=$(samtools view -c -F4 "$GBAM" $(echo "Chr1 Chr2 Chr3 Chr4 Chr5"))
declare -A gReg
for reg in $regions; do gReg[$reg]=$(samtools view -c -F4 -L $R/${reg}.bed "$GBAM"); done

GMB=$(awk '$1~/^Chr[1-5]$/{x+=$2}END{print x/1e6}' $R/genome_chr15.sizes)
echo "sample,region,bliss_reads,gdna_reads,region_Mb,fold_nogdna,gdna_fold,enrichment_vs_gdna,log2" > "$OUT"
for s in BA1 BA2 old_BA1_BA2; do
  bam=$P/results_TAIR12/bam/${s}.${BSUF}.bam
  bTot=$(samtools view -c -F4 "$bam" Chr1 Chr2 Chr3 Chr4 Chr5)
  for reg in $regions; do
    b=$(samtools view -c -F4 -L $R/${reg}.bed "$bam"); g=${gReg[$reg]}
    rb=$(awk '{x+=$3-$2}END{print x/1e6}' $R/${reg}.bed)
    awk -v s=$s -v r=$reg -v b=$b -v g=$g -v bT=$bTot -v gT=$gTot -v rb=$rb -v gmb=$GMB 'BEGIN{
      foldB=(b/rb)/(bT/gmb);          # BLISS over uniform (no gDNA)
      foldG=(g/rb)/(gT/gmb);          # gDNA over uniform
      e=(foldG>0)?foldB/foldG:0;       # BLISS/gDNA = enrichment over gDNA expectation
      printf "%s,%s,%d,%d,%.2f,%.2f,%.2f,%.3f,%.3f\n",s,r,b,g,rb,foldB,foldG,e,(e>0?log(e)/log(2):0)}' >> "$OUT"
  done
done
echo "=== BLISS DSB enrichment over WGA gDNA (per compartment) ==="
column -t -s, "$OUT"
echo ""; echo "=== BA1: fold_nogdna vs gdna_fold vs BLISS/gDNA enrichment ==="
awk -F, 'NR>1 && $1=="BA1"{printf "%-15s  no-gDNA=%5.2fx  gDNA=%5.2fx  BLISS/gDNA=%5.2fx\n",$2,$6,$7,$8}' "$OUT" | sort -t= -k4 -rn