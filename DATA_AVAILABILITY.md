# Data availability

This repository contains **code, documentation, self-contained HTML reports, and
small derived text results** (region BEDs, hotspot tables, summary CSVs, figures).
The large primary and intermediate data (raw FASTQ, BAM, bigWig, genome FASTA and
aligner indices) are **not** in git — they are regenerable from the pipeline or
available from the sources below.

## Raw sequencing reads (sBLISS-seq, this study)

| Delivery (Novogene) | Libraries | Notes |
|---------------------|-----------|-------|
| `X204SC26044732-Z01-F001` | BA1, BA2, old_BA1_BA2 | primary DSB-mapping run |
| `X204SC26053809-Z01-F001` | BA1_WT, BA2_OX_40, BA2_OX_80_1/2 | CENH3-overexpression comparison |

Reads are PE150 (R1 + R2); only **R1** is used (R2 is the T7/backbone side and is
uninformative for the break position). To be deposited to ENA/SRA on publication.

## Reference genomes

| Genome | Accession / source | Used for |
|--------|--------------------|----------|
| TAIR12 / Col-CC (T2T, rDNA-resolved) | `GCA_028009825.2` (`..._Col-CC_genomic_withorganelles.fna`) | primary analysis (`results_TAIR12`) |
| Col-0 ragtag (Chr1–5 + organelles) | project assembly `Col-0.ragtag_chrs.mito.chloro.renamed.fa` | earlier run (`results`) |

Build aligner indices locally: `bowtie2-build genome.fna <prefix>` /
`hisat2-build genome.fa <prefix>`, then point `reference.*_index` in the config at them.

## Annotation

- `TAIR12_1Feb26.gff3` — merged Gnomon + TAIR12 genes, EDTA TEs, TRASH satellites,
  RepeatMasker rRNA (used to derive `results/regions/*.bed`).
- Region BEDs (centromere, pericentromere, arms, 45S/5S rDNA, CEN178, TEs, genes,
  lncRNA, TRASH satellites, G4, methylation, mappability) are included under
  `results/regions/` so the downstream analyses can be reproduced without the GFF.

## External datasets used for cross-validation (public)

| Signal | Accession |
|--------|-----------|
| RAD51 (WT) | `SRR24938940` |
| RPA1A (WT) | `SRR24938935` |
| Input (for RAD51/RPA1A) | `SRR24938942` |
| ssDRIP-seq (R-loops) | `SRR8097477` |
| γH2A.X (WT) | `SRR8434381` |
| G4-seq (OQS peaks) | GEO `GSE110582` (Marsico et al. 2019) |
| WGA gDNA (naked-DNA control) | `E-MTAB-6257` (ERR2215864 / ERR2215865, random/WGA) |

The WGA gDNA control is central to the analysis: it distinguishes genuine DSB
enrichment from copy-number / mappability artifacts (see the reports, §4–§5, §7).
