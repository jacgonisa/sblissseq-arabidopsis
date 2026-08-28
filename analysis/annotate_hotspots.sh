#!/usr/bin/env bash
# Annotate DSB hotspots with overlapping genes/TEs and assess cross-sample
# reproducibility.
set -euo pipefail

HS=${BLISS_HSDIR:-/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/hotspots}
REG=${BLISS_REGDIR:-/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/regions}
OUT=$HS

# ── 1. Reproducibility: how many hotspots are shared across samples? ──────────
echo "=== Hotspot reproducibility (1kb tiles) ==="
for s in BA1 BA2 old_BA1_BA2; do
  # keep BED6 (cols 1-6); bedtools rejects the float fold in col7
  cut -f1-6 $HS/${s}_1kb.hotspots.bed | sort -k1,1 -k2,2n > $HS/${s}.sorted.bed
done

ba1_ba2=$(bedtools intersect -a $HS/BA1.sorted.bed -b $HS/BA2.sorted.bed -u | wc -l)
ba1_n=$(wc -l < $HS/BA1.sorted.bed)
ba2_n=$(wc -l < $HS/BA2.sorted.bed)
echo "BA1 hotspots: $ba1_n | BA2 hotspots: $ba2_n | shared: $ba1_ba2"

# 3-way intersect → robust hotspot set
bedtools intersect -a $HS/BA1.sorted.bed -b $HS/BA2.sorted.bed -u | \
  bedtools intersect -a - -b $HS/old_BA1_BA2.sorted.bed -u \
  > $HS/robust_hotspots.bed
echo "Robust hotspots (in all 3 samples): $(wc -l < $HS/robust_hotspots.bed)"

# ── 2. Annotate robust hotspots with genes ───────────────────────────────────
echo ""
echo "=== Annotating robust hotspots ==="
bedtools intersect -a $HS/robust_hotspots.bed -b $REG/genes.bed6 -wa -wb 2>/dev/null | \
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,$(NF-2)}' | sort -k4,4nr > $OUT/robust_hotspots.genes.tsv

genic=$(bedtools intersect -a $HS/robust_hotspots.bed -b $REG/genes.bed6 -u | wc -l)
te=$(bedtools intersect -a $HS/robust_hotspots.bed -b $REG/TEs.bed6 -u | wc -l)
nor=$(bedtools intersect -a $HS/robust_hotspots.bed -b $REG/nor.bed -u 2>/dev/null | wc -l)
total=$(wc -l < $HS/robust_hotspots.bed)
echo "Robust hotspots overlapping genes:        $genic / $total"
echo "Robust hotspots overlapping TEs:          $te / $total"
echo "Robust hotspots overlapping NOR 45S rDNA: $nor / $total  (repetitive-pileup artefacts)"

# ── 3. Where do robust hotspots fall (arms/peri/cen)? ────────────────────────
# NOTE: no blacklist is applied — centromeres are analysed, not removed; NOR
# 45S rDNA tiles are kept but flagged above as artefacts.
echo ""
echo "=== Robust hotspot region distribution ==="
cen=$(bedtools intersect -a $HS/robust_hotspots.bed -b $REG/centromere.bed -u 2>/dev/null | wc -l)
peri=$(bedtools intersect -a $HS/robust_hotspots.bed -b $REG/pericentromere.bed -u 2>/dev/null | wc -l)
# arms = robust hotspots NOT in peri and NOT in cen
arm=$(bedtools intersect -a $HS/robust_hotspots.bed -b $REG/pericentromere.bed $REG/centromere.bed -v 2>/dev/null | wc -l)
echo "arm:            $arm"
echo "pericentromere: $peri"
echo "centromere:     $cen"

# ── 4. Region-annotated robust hotspot table ─────────────────────────────────
# chr start end score region is_NOR top_gene  (region: arm/pericentromere/centromere)
echo ""
echo "=== Writing region-annotated robust hotspot table ==="
# per-hotspot top overlapping gene (first = highest score, table is score-sorted)
awk 'BEGIN{OFS="\t"} !seen[$1":"$2":"$3]++ {print $1":"$2":"$3, $5}' \
    $OUT/robust_hotspots.genes.tsv > $HS/.gene_lookup.tsv
# region + NOR flags via counts, then resolve label
bedtools intersect -a $HS/robust_hotspots.bed -b $REG/centromere.bed     -c 2>/dev/null | \
  bedtools intersect -a - -b $REG/pericentromere.bed -c 2>/dev/null | \
  bedtools intersect -a - -b $REG/nor.bed            -c 2>/dev/null | \
  awk 'BEGIN{OFS="\t"}
       {key=$1":"$2":"$3;
        region=($7>0)?"centromere":(($8>0)?"pericentromere":"arm");
        isnor=($9>0)?"TRUE":"FALSE";
        print key, $1, $2, $3, $5, region, isnor}' | \
  sort -k5,5nr > $HS/.robust_tmp.tsv
echo -e "chr\tstart\tend\tscore\tregion\tis_NOR\ttop_gene" > $OUT/robust_hotspots.annotated.tsv
awk 'BEGIN{OFS="\t"}
     NR==FNR{g[$1]=$2; next}
     {gene=($1 in g)?g[$1]:".";
      print $2,$3,$4,$5,$6,$7,gene}' \
     $HS/.gene_lookup.tsv $HS/.robust_tmp.tsv >> $OUT/robust_hotspots.annotated.tsv
rm -f $HS/.gene_lookup.tsv $HS/.robust_tmp.tsv
echo "Wrote $OUT/robust_hotspots.annotated.tsv ($(($(wc -l < $OUT/robust_hotspots.annotated.tsv)-1)) hotspots)"
echo "  centromere rows: $(awk -F'\t' 'NR>1 && $5=="centromere"' $OUT/robust_hotspots.annotated.tsv | wc -l)"
echo "  is_NOR=TRUE rows: $(awk -F'\t' 'NR>1 && $6=="TRUE"' $OUT/robust_hotspots.annotated.tsv | wc -l)"

# ── 5. Top 20 robust hotspot genes ───────────────────────────────────────────
echo ""
echo "=== Top 20 genes at robust hotspots (by break count) ==="
head -20 $OUT/robust_hotspots.genes.tsv

# cleanup
rm -f $HS/*.sorted.bed
echo ""
echo "Outputs: robust_hotspots.bed, robust_hotspots.genes.tsv, robust_hotspots.annotated.tsv"
