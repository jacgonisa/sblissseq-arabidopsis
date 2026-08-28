# sBLISS-seq Pipeline — Arabidopsis thaliana DSB Mapping

Snakemake pipeline for mapping physiological DNA double-strand breaks (DSBs)
from **in-suspension BLISS sequencing (sBLISS)** data, adapted from
[Hidmi et al. 2024 (STAR Protocols)](https://doi.org/10.1016/j.xpro.2024.103059)
for *Arabidopsis thaliana*.

---

## Overview

sBLISS labels DSB ends *in situ* with barcoded adapters containing a T7
promoter, followed by *in vitro* transcription, reverse transcription, and
Illumina sequencing. This pipeline processes the resulting paired-end reads
into single-nucleotide DSB maps.

### Pipeline steps

```
Raw FASTQ (R1 only — R2 contains adapter only)
    │
    ├─ 1. FastQC + MultiQC          raw read QC (duplication, adapter content)
    │
    ├─ 2. bliss_demux.py            extract UMI + verify Internal_Index (barcode)
    │       UMI            = R1[:8]
    │       Internal_Index = R1[8:16]
    │       header → @readname+Internal_Index+UMI
    │       sequence → R1[16:]  (16 bp stripped)
    │       Hamming distance ≤ 1 mismatch allowed on barcode
    │
    ├─ 3. Trim Galore -q 20 --length 20
    │       quality trim, adapter removal (RA3), poly-G tail removal
    │
    ├─ 4. FastQC + MultiQC          post-trim QC
    │
    ├─ 5. HISAT2 --no-spliced-alignment --no-unal
    │       DNA-mode alignment to Col-0 genome
    │       samtools view | sort | index
    │
    ├─ 6. umi_tools dedup --umi-separator="+"
    │       PCR duplicate removal using UMI + 5' alignment position
    │       --output-stats for per-UMI dedup statistics
    │
    ├─ 7. DSB site extraction
    │       5' end of each deduplicated R1 read = exact DSB coordinate
    │       output: tabix-indexed BED + per-position counts
    │
    ├─ 8. bamCoverage (deepTools)
    │       CPM-normalised bigWig (1 bp resolution) for IGV/UCSC
    │
    └─ 9. MultiQC full              all logs summarised
```

---

## Read structure

```
R1:  [UMI: 8 bp][Internal_Index: 8 bp][genomic sequence from DSB][RA3 adapter if short insert]
R2:  [barcode: 8 bp][UMI: 8 bp][GATCGTCG linker][T7 backbone]  ← not used after demux
```

### Internal_Index (barcode) sequences

| Adapter | Internal_Index in R1 | Notes |
|---------|---------------------|-------|
| BA1 | `CATCACGC` | BA1 top strand 3' blunt end |
| BA2 | `GTCGTCGC` | BA2 top strand 3' blunt end |

---

## Samples (June 2026 run — X204SC26044732-Z01-F001)

| Sample | Library ID | Internal_Index | Notes |
|--------|-----------|----------------|-------|
| BA1 | EKDL260003554 | CATCACGC | BA1 adapter |
| BA2 | EKDL260003555 | GTCGTCGC | BA2 adapter — 2 lanes merged |
| old_BA1_BA2 | EKDL260003553 | GTCGTCGC | Named "oldBA1" but BA2 adapter used |

---

## Reference genome

**Col-0** (`Col-0.ragtag_chrs.mito.chloro.renamed.fa`)
- Chromosomes: Chr1–Chr5
- Mitochondria: NC_037304.1
- Chloroplast: NC_000932.1

HISAT2 index: `reference/hisat2_index/Col-0`

---

## Requirements

```bash
conda env create -f envs/bliss.yaml
conda activate bliss
```

Key tools: `hisat2`, `samtools`, `umi_tools`, `trim-galore`, `fastqc`,
`multiqc`, `deeptools`, `bedtools`, `snakemake`

---

## Usage

```bash
# 1. QC only (raw + post-trim)
snakemake --configfile config.yaml --cores 8 \
    results/multiqc_raw/multiqc_report.html \
    results/multiqc_trimmed/multiqc_report.html

# 2. Full pipeline
snakemake --configfile config.yaml --cores 8

# 3. Specific sample
snakemake --configfile config.yaml --cores 8 results/bam/BA1.dedup.bam
```

---

## Configuration (`config.yaml`)

| Key | Description |
|-----|-------------|
| `samples.<name>.r1` | Path to R1 FASTQ |
| `samples.<name>.internal_index` | Expected 8 bp barcode |
| `adapters.max_barcode_mismatches` | Hamming distance tolerance (default 1) |
| `reference.hisat2_index` | HISAT2 index prefix |
| `params.threads` | Threads per job |

---

## Outputs

```
results/
├── multiqc_raw/         raw read QC report
├── multiqc_trimmed/     post-trim QC report (FastQC + Trim Galore)
├── logs/
│   ├── *.demux.stats.txt       barcode pass/fail counts
│   ├── *.trimgalore.log        Trim Galore summary
│   ├── *.hisat2.summary        alignment rates
│   └── *.dedup.log / .stats    UMI dedup statistics
├── bam/
│   ├── <sample>.bam            aligned reads (HISAT2)
│   └── <sample>.dedup.bam      PCR-deduplicated reads
├── dsb/
│   ├── <sample>.dsb.bed.gz     tabix-indexed DSB sites (BED6)
│   └── <sample>.dsb.stats.txt  total reads / unique sites / mean depth
└── bigwig/
    └── <sample>.bw             CPM-normalised coverage (1 bp bins)
```

---

## Key QC metrics to check

| Metric | Expected range |
|--------|---------------|
| Barcode pass rate (`bliss_demux`) | > 80% |
| Reads with adapter (Trim Galore) | 30–50% (normal for sBLISS short inserts) |
| HISAT2 overall alignment rate | > 70% |
| UMI dedup retention | 20–80% (depends on library complexity) |
| Mean reads per DSB site | ~1.0–1.5 (low = good complexity) |

---

---

## Downstream analysis

### Annotation files used

| File | Source | Purpose |
|------|--------|---------|
| `/mnt/ssd-4tb/HIFI_NAMIL/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds.fa_liftoff.edta.gff3` | Liftoff (TAIR10→ragtag) + EDTA | 27,953 genes + 46,707 TEs — **use this for IGV and all analyses** |
| `/home/jg2070/Desktop/PhD/PoreC/experiment_results/2025September/centromeres.bed` | PoreC project | Centromere core coordinates (Col-0 ragtag) |
| `/home/jg2070/Desktop/PhD/PoreC/experiment_results/2025September/data/Col-0.ragtag_scaffolds_ARMS.bed` | PoreC project | p-arm / centromere / q-arm boundaries |

Pericentromere defined as centromere ±2 Mb (clamped to chromosome bounds), consistent with TAIR12 convention.

---

### DSB enrichment results (HISAT2, BA1)

#### Genomic features (ragtag annotation)

| Feature | Genome bg | DSB % | Fold enrichment |
|---------|----------|-------|----------------|
| Genes | 44.6% | 40.0% | **0.90×** (depleted) |
| TEs | 24.5% | 30.1% | **1.23×** (enriched) |

DSBs are consistently depleted in genic regions and enriched in TE regions across all three samples. BA1 > BA2 > old_BA1_BA2 in terms of TE enrichment magnitude.

#### Chromosomal compartments (reconciled — `reconcile_density.R`)

Three normalisations, same numerator (distinct break sites), **NOR rDNA always its own category**:

| Compartment | (a) sites / Mb DNA | (b) sites / 1k reads | (c) sites / mapped-Mb |
|-------------|-------------------:|---------------------:|----------------------:|
| Arm (NOR-excl.) | 33,285 | **949** | 7,686 |
| Pericentromere  | 53,531 | 859 | 6,908 |
| Centromere core | 25,470 | 939 | 7,447 |
| rDNA (NOR)      | 489,120 | **226** | 1,735 |

> **(b) and (c) are the same metric** — they differ only by mean read length (~125 bp), so
> (c) ≈ (b) × 7.8; across all region×sample points r = **0.998**. They give the **same ranking**:
> **Arm ≈ Centromere > Pericentromere**, with rDNA the lone outlier.

> **Correction:** an earlier version reported "centromere > pericentromere > arms (939 vs 549)".
> The 549 was an artefact — the "arms" region had **included the NOR rDNA** (a low-complexity
> pile-up: 4.3M reads → ~1.0M sites = 226/1k). With NOR separated, **arms = 949 ≈ centromere**.
> Per mappable read, genuine compartments are uniform within ~10%; the centromere is **not**
> specially enriched (nor depleted). Raw per-Mb (a) only looks heterochromatin-skewed because
> more sequence maps there. See `figures/density_reconciliation.png`, `figures/cenh3_zoom.png`.

> **rDNA paradox:** NOR looks *enriched* in raw counts but *depleted* per read. Both are artefacts
> of the collapsed-repeat assembly + 8 bp UMI saturation (see "rDNA / NOR" note below) — rDNA
> fragility is real but **unquantifiable** with this assay/assembly.

#### Plot

```
results/dsb_enrichment.pdf   — 3-panel figure (raw fold, coverage-normalised, gene/TE)
results/dsb_enrichment.png
results/dsb_enrichment.csv   — raw data table
```

Generated by: `pipeline/scripts/plot_dsb_enrichment.R`

---

---

## Comprehensive downstream analysis

Full report: **`results/analysis/sBLISS_report.html`** (self-contained, all figures embedded).

### Scripts (`pipeline/scripts/`)

| Script | Purpose |
|--------|---------|
| `cpm_normalize_breaks.sh` | CPM-normalise BREAK bedgraphs → bigWig (cross-sample comparable) |
| `genome_density.R` | Chromosome-scale DSB density landscape; takes a tile-width arg (run at 100 kb / 10 kb / 1 kb), regions overlaid, no masking |
| `sample_correlation.R` | Pairwise sample correlation (10 kb bins, CPM; Pearson/Spearman; NOR bins flagged) |
| `hotspots_venn_classification.R` | 3-way hotspot Venn (all + rDNA-excluded) and compartment classification (arm/peri/cen/rDNA) |
| `dsb_per_mapped_mb.sh` + `dsb_per_mapped_mb_plot.R` | DSB sites per Mb of uniquely-mapped reads (bedcov aligned bases); NOR separate |
| `reconcile_density.R` | Reconciles the 3 density normalisations with consistent NOR-separated regions |
| `cenh3_zoom.R` | Centromere zoom overlaying sBLISS DSB density with CENH3 ChIP-seq (ragtag, 10 kb) |
| `metaplots.sh` | deepTools profiles/heatmaps at TSS, TES, gene bodies, TEs, TE superfamilies |
| `hotspots_brgenomics.R` | BRGenomics tiling + Poisson hotspot calling (raw counts; optional blacklist; region + NOR flags) |
| `annotate_hotspots.sh` | Cross-sample reproducibility + gene/TE/region annotation; emits `robust_hotspots.annotated.tsv` |
| `plot_dsb_enrichment.R` | Region & feature enrichment bar plots (raw + coverage-normalised) |
| `build_report.py` | Assemble the self-contained HTML report |

### No blacklist — regions are flagged, not masked
Hotspot calling and the density landscape apply **no blacklist**. Centromeres and the
NOR2/NOR4 45S rDNA arrays (Chr2/Chr4 0–1 Mb) are analysed like any other region; every
hotspot is **annotated** with its compartment (arm/pericentromere/centromere) and an
`is_NOR` flag instead (`hotspots/robust_hotspots.annotated.tsv`). The NOR arrays are
collapsed-rDNA repeat pile-ups and dominate the strongest raw hotspots — they are kept
visible (flagged) so they can be excluded at interpretation rather than hidden.
`regions/nor.bed` and `regions/centromere.bed` hold these annotation intervals;
`regions/blacklist_full.bed` is retired (kept only for provenance). Centromeres are
analysed region-based (BED), **not** via gene/TE GFF.

#### Why rDNA / NOR is enriched (raw) but depleted (per-read) — and why it's a blind spot
The 45S rDNA is genuinely one of the most fragile loci, but **neither metric quantifies it**:
- **Assembly collapse:** the hundreds–thousands of real rDNA copies are represented by a ~1 Mb
  stub. Reads from *all* copies pile onto that stub → huge raw break counts (NOR dominates raw
  hotspots; 489k sites/Mb). This **over-states** rDNA (it's copy-number pile-up, not a per-copy rate).
- **Positional + UMI saturation:** with only ~2 Mb of assembled rDNA, nearly every position is
  already a break site, and the 8 bp UMI (65,536 combos) saturates at ultra-high-coverage positions
  → `umi_tools` collapses distinct real breaks into one. So *new* unique sites stop accumulating:
  NOR has 4.4 reads per unique site vs 1.05 on the arms. That makes per-read density look low
  (226 vs 949 sites/1k reads) — an artefact, **not** evidence rDNA is unbroken.
- **Bottom line:** rDNA fragility is real and consistent with the raw signal, but it is
  **unquantifiable** here. Proper measurement needs an rDNA-resolved assembly (or per-copy-unit
  normalisation) and a longer UMI. This is exactly why rDNA is flagged separately.

### Key biological findings
- **Compartment density (reconciled, NOR separated):** per mappable read/Mb,
  **Arm ≈ Centromere > Pericentromere** — uniform within ~10%; centromere is not specially
  enriched. (b) per-1k-reads and (c) per-mapped-Mb are proportional (r=0.998). Raw per-Mb (a)
  over-weights heterochromatin via mappability. **rDNA (NOR) is a blind spot** — enriched raw,
  depleted per-read; both artefacts of collapsed-repeat assembly + UMI saturation.
- **Reproducibility:** BA1↔BA2 Pearson r=0.90 (10 kb, log CPM, NOR excluded); old_BA1_BA2 weaker (0.69–0.82).
  Hotspot 3-way overlap: 1,886 all / 581 genuine (rDNA-excluded). See `figures/hotspots_venn.png`,
  `figures/hotspots_classification.png`, `figures/sample_correlation.png`, `figures/genome_density_{100kb,10kb,1kb}.png`.
- **Features:** DSBs depleted in genes (0.81–0.90×), enriched in TEs (1.23–1.62×).
- **Hotspots (no blacklist, raw counts):** 1,886 robust (3-way reproducible) 1-kb hotspots.
  **1,305 are NOR 45S rDNA artefacts** (`is_NOR=TRUE`); the remaining ~581 are genuine —
  490 pericentromeric, 57 **centromeric** (now recovered), 34 on arms. Strongest genuine
  hotspots cluster in Chr2 pericentromere (~4.5 Mb); top arm *gene* hotspot = AT5G10250
  (Chr5:3.26 Mb, 2,137 breaks); top centromere hotspot = Chr4:7.03 Mb (593 breaks).
- **BA1 ≈ BA2** (concordant); **old_BA1_BA2** is low-complexity (90% dup) — qualitative only.

### Hotspot detection (BRGenomics, paper snippet — blacklist optional)
```r
smp      <- import(break_track, format="bedGraph")   # RAW per-base break counts
# blacklist is OPTIONAL: pass "none" to skip filtering (centromeres + NORs analysed)
tiles    <- tileGenome(chrom_sizes, tilewidth=1000)
counts   <- getCountsByRegions(smp, tiles, field="score")
# → Poisson test per tile vs genome mean, BH-adjust, call padj<0.001 & fold>=3
# → annotate each hotspot with region (arm/peri/cen) + is_NOR flag
```
Run with: `Rscript hotspots_brgenomics.R <break.bedgraph> Col-0.chrom.sizes none 1000 <out_prefix>`

---

## Reference

Hidmi O, Oster S, Shatleh D, Monin J, Aqeilan RI.
*Protocol for mapping physiological DSBs using in-suspension break labeling
in situ and sequencing.*
STAR Protocols 5, 103059 (2024).
https://doi.org/10.1016/j.xpro.2024.103059
