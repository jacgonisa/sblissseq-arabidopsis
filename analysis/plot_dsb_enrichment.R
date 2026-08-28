#!/usr/bin/env Rscript
# Plot DSB enrichment across genomic features and chromosomal regions
# Two panels:
#   Left:  raw fold-enrichment (DSBs per bp, relative to genome background)
#   Right: coverage-normalised (DSBs per 1000 mapped reads)
# Input: CSV produced by the shell collection script

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
csv_in  <- ifelse(length(args) >= 1, args[1],
  "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/dsb_enrichment.csv")
out_pdf <- ifelse(length(args) >= 2, args[2],
  "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/dsb_enrichment.pdf")
out_png <- sub("\\.pdf$", ".png", out_pdf)

df <- read.csv(csv_in)

# ── tidy labels ──────────────────────────────────────────────────────────────
df$region <- factor(df$region,
  levels = c("centromere", "pericentromere", "arms"),
  labels = c("Centromere", "Pericentromere", "Arms"))

df$sample <- factor(df$sample,
  levels = c("BA1", "BA2", "old_BA1_BA2"),
  labels = c("BA1", "BA2", "old_BA1_BA2"))

# ── compute fold-enrichment over genome background ───────────────────────────
genome_bp <- 134556983
df <- df %>%
  mutate(
    genomic_frac  = genomic_bp / genome_bp,
    dsb_frac      = dsb_sites  / sum(dsb_sites[region == region]),
    dsb_per_bp    = dsb_sites  / genomic_bp * 1e6,  # per Mb
    fold_raw      = (dsb_sites / sum(dsb_sites)) / (genomic_bp / genome_bp)
  ) %>%
  group_by(sample) %>%
  mutate(
    total_dsb     = sum(dsb_sites),
    fold_raw      = (dsb_sites / total_dsb) / (genomic_bp / genome_bp)
  ) %>%
  ungroup()

colours <- c("Centromere" = "#d62728", "Pericentromere" = "#ff7f0e", "Arms" = "#1f77b4")

# ── Panel A: raw fold enrichment ─────────────────────────────────────────────
pA <- ggplot(df, aes(x = sample, y = fold_raw, fill = region)) +
  geom_col(position = "dodge", width = 0.7, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = colours) +
  labs(
    title    = "A  Raw fold-enrichment (by genomic bp)",
    subtitle = "Dashed line = expected if uniform",
    x = NULL, y = "Fold enrichment over background",
    fill = "Region"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 11)
  )

# ── Panel B: coverage-normalised density ─────────────────────────────────────
pB <- ggplot(df, aes(x = sample, y = dsb_per_1k_mapped, fill = region)) +
  geom_col(position = "dodge", width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = colours) +
  labs(
    title    = "B  Coverage-normalised DSB density",
    subtitle = "DSBs per 1,000 mapped reads (corrects for mapping bias)",
    x = NULL, y = "DSBs per 1,000 mapped reads",
    fill = "Region"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 11)
  )

# ── Panel C: gene + TE enrichment ────────────────────────────────────────────
# Load gene/TE numbers if available
feat_data <- data.frame(
  sample  = rep(c("BA1", "BA2", "old_BA1_BA2"), each = 2),
  feature = rep(c("Genes", "TEs"), 3),
  fold    = c(0.90, 1.23,   # BA1
              0.86, 1.38,   # BA2
              0.81, 1.62)   # old_BA1_BA2
)
feat_data$sample  <- factor(feat_data$sample, levels = c("BA1","BA2","old_BA1_BA2"))
feat_data$feature <- factor(feat_data$feature, levels = c("Genes","TEs"))

pC <- ggplot(feat_data, aes(x = sample, y = fold, fill = feature)) +
  geom_col(position = "dodge", width = 0.7, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = c("Genes" = "#2ca02c", "TEs" = "#9467bd")) +
  labs(
    title    = "C  DSB enrichment at genomic features",
    subtitle = "Background: genes=44.6%, TEs=24.5% of genome",
    x = NULL, y = "Fold enrichment over background",
    fill = "Feature"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 11)
  )

# ── Combine and save ──────────────────────────────────────────────────────────
combined <- (pA | pB) / pC +
  plot_annotation(
    title   = "sBLISS-seq DSB Enrichment Analysis — Arabidopsis thaliana Col-0",
    subtitle = "HISAT2 alignment, UMI-deduplicated. Annotation: Liftoff+EDTA (Col-0 ragtag).",
    theme   = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(out_pdf, combined, width = 12, height = 9)
ggsave(out_png, combined, width = 12, height = 9, dpi = 150)
cat("Saved:", out_pdf, "\n")
cat("Saved:", out_png, "\n")
