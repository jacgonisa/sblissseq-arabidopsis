#!/usr/bin/env bash
# deepTools metaplots of sBLISS DSB signal at genomic features.
# Uses CPM-normalised bigWigs (Chr1-5) for cross-sample comparison.
set -euo pipefail

REG=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/regions
BW=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/bigwig_cpm
MAT=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/metaplots
FIG=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/figures
THREADS=8

mkdir -p "$MAT" "$FIG"

BWS="$BW/BA1.cpm.bw $BW/BA2.cpm.bw $BW/old_BA1_BA2.cpm.bw"
LABELS="BA1 BA2 old_BA1_BA2"

# ── 1. TSS reference-point (±2 kb) ────────────────────────────────────────────
computeMatrix reference-point --referencePoint TSS \
    -S $BWS -R "$REG/genes.bed6" \
    -b 2000 -a 2000 --binSize 50 \
    --missingDataAsZero -p $THREADS \
    -o "$MAT/TSS.matrix.gz" 2>/dev/null
plotProfile -m "$MAT/TSS.matrix.gz" -o "$FIG/profile_TSS.png" \
    --samplesLabel $LABELS --refPointLabel TSS \
    --plotTitle "DSB density around TSS" --perGroup 2>/dev/null

# ── 2. TES reference-point (±2 kb) ────────────────────────────────────────────
computeMatrix reference-point --referencePoint TES \
    -S $BWS -R "$REG/genes.bed6" \
    -b 2000 -a 2000 --binSize 50 \
    --missingDataAsZero -p $THREADS \
    -o "$MAT/TES.matrix.gz" 2>/dev/null
plotProfile -m "$MAT/TES.matrix.gz" -o "$FIG/profile_TES.png" \
    --samplesLabel $LABELS --refPointLabel TES \
    --plotTitle "DSB density around TES" --perGroup 2>/dev/null

# ── 3. Gene body scale-regions (TSS..TES + 2 kb flanks) ───────────────────────
computeMatrix scale-regions \
    -S $BWS -R "$REG/genes.bed6" \
    -b 2000 -a 2000 --regionBodyLength 4000 --binSize 50 \
    --missingDataAsZero -p $THREADS \
    -o "$MAT/genebody.matrix.gz" 2>/dev/null
plotProfile -m "$MAT/genebody.matrix.gz" -o "$FIG/profile_genebody.png" \
    --samplesLabel $LABELS --startLabel TSS --endLabel TES \
    --plotTitle "DSB density over gene bodies" --perGroup 2>/dev/null
plotHeatmap -m "$MAT/genebody.matrix.gz" -o "$FIG/heatmap_genebody.png" \
    --samplesLabel $LABELS --startLabel TSS --endLabel TES \
    --colorMap viridis 2>/dev/null

# ── 4. TE scale-regions (all repeat_regions) ──────────────────────────────────
computeMatrix scale-regions \
    -S $BWS -R "$REG/TEs.bed6" \
    -b 1000 -a 1000 --regionBodyLength 2000 --binSize 50 \
    --missingDataAsZero -p $THREADS \
    -o "$MAT/TE.matrix.gz" 2>/dev/null
plotProfile -m "$MAT/TE.matrix.gz" -o "$FIG/profile_TE.png" \
    --samplesLabel $LABELS --startLabel "TE start" --endLabel "TE end" \
    --plotTitle "DSB density over transposons" --perGroup 2>/dev/null

# ── 5. Per-superfamily TE profile (BA1 only, to compare TE classes) ──────────
TE_BEDS=""
TE_LABELS=""
for te in Gypsy_LTR_retrotransposon Copia_LTR_retrotransposon Mutator_TIR_transposon helitron L1_LINE_retrotransposon; do
    TE_BEDS="$TE_BEDS $REG/TE_${te}.bed6"
    TE_LABELS="$TE_LABELS ${te%%_*}"
done
computeMatrix scale-regions \
    -S "$BW/BA1.cpm.bw" -R $TE_BEDS \
    -b 1000 -a 1000 --regionBodyLength 2000 --binSize 50 \
    --missingDataAsZero -p $THREADS \
    -o "$MAT/TE_superfamily.matrix.gz" 2>/dev/null
plotProfile -m "$MAT/TE_superfamily.matrix.gz" -o "$FIG/profile_TE_superfamily.png" \
    --regionsLabel $TE_LABELS --startLabel start --endLabel end \
    --plotTitle "BA1 DSB density by TE superfamily" --perGroup 2>/dev/null

echo "Metaplots complete. Figures in $FIG/"
