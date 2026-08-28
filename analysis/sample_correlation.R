#!/usr/bin/env Rscript
# Pairwise sample correlation of genome-wide DSB density (CPM, 10 kb bins).
# NOR 45S rDNA bins are flagged so the rDNA pile-ups don't silently inflate r.
# Usage: Rscript sample_correlation.R [out_png]

options(scipen = 999)
suppressMessages({
  library(rtracklayer); library(BRGenomics); library(GenomicRanges)
  library(ggplot2); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
# path overrides (default ragtag; set env vars for TAIR12)
regdir  <- Sys.getenv("BLISS_REGDIR",   file.path(base, "results/analysis/regions"))
bwdir   <- Sys.getenv("BLISS_BWDIR",    file.path(base, "results/analysis/bigwig_cpm"))
chrfile <- Sys.getenv("BLISS_CHRSIZES", file.path(base, "reference/Col-0.chrom.sizes"))
out_png <- ifelse(length(args) >= 1, args[1],
  file.path(base, "results/analysis/figures/sample_correlation.png"))

chrsz <- read.table(chrfile)
chrsz <- chrsz[grepl("^Chr[1-5]$", chrsz$V1), ]
seqlens <- setNames(chrsz$V2, chrsz$V1)
tiles <- tileGenome(seqlens, tilewidth = 10000, cut.last.tile.in.chrom = TRUE)

samples <- c("BA1", "BA2", "old_BA1_BA2")
mat <- sapply(samples, function(s) {
  bw <- import(file.path(bwdir, sprintf("%s.cpm.bw", s)), format = "bigwig")
  getCountsByRegions(bw, tiles, field = "score")
})
df <- as.data.frame(mat)

# flag NOR bins
nor <- import(file.path(regdir, "nor.bed"), format = "bed")
df$is_NOR <- FALSE
df$is_NOR[queryHits(findOverlaps(tiles, nor))] <- TRUE

# keep informative bins
df <- df[rowSums(df[samples]) > 0, ]

# ── correlation values (genuine bins only = NOR excluded) ─────────────────────
g <- df[!df$is_NOR, samples]
cor_tab <- function(method) {
  m <- cor(log10(g + 1), method = method)
  round(m, 3)
}
cat("Pearson (log10, NOR-excluded):\n");  print(cor_tab("pearson"))
cat("\nSpearman (NOR-excluded):\n");       print(cor_tab("spearman"))

# ── pairwise scatter panels ───────────────────────────────────────────────────
pairs_df <- bind_rows(
  data.frame(x = df$BA1, y = df$BA2,         pair = "BA1 vs BA2",         is_NOR = df$is_NOR),
  data.frame(x = df$BA1, y = df$old_BA1_BA2, pair = "BA1 vs old_BA1_BA2", is_NOR = df$is_NOR),
  data.frame(x = df$BA2, y = df$old_BA1_BA2, pair = "BA2 vs old_BA1_BA2", is_NOR = df$is_NOR)
)
# annotate Pearson r (NOR-excluded) per pair
rlab <- function(a, b) sprintf("r = %.3f", cor(log10(g[[a]] + 1), log10(g[[b]] + 1)))
labs_df <- data.frame(
  pair = c("BA1 vs BA2", "BA1 vs old_BA1_BA2", "BA2 vs old_BA1_BA2"),
  lab  = c(rlab("BA1","BA2"), rlab("BA1","old_BA1_BA2"), rlab("BA2","old_BA1_BA2")))

p <- ggplot(pairs_df, aes(x + 1, y + 1)) +
  geom_point(aes(colour = is_NOR), size = 0.4, alpha = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("FALSE" = "#1f77b4", "TRUE" = "#d62728"),
                      labels = c("FALSE" = "genuine bin", "TRUE" = "NOR 45S rDNA"),
                      name = NULL) +
  facet_wrap(~pair) +
  geom_text(data = labs_df, aes(x = 1.5, y = max(pairs_df$y + 1), label = lab),
            hjust = 0, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE) +
  labs(title = "sBLISS sample correlation (DSB density, 10 kb bins, CPM)",
       subtitle = "Pearson r computed on log10(CPM+1), NOR 45S rDNA bins excluded (shown red)",
       x = "CPM + 1", y = "CPM + 1") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

ggsave(out_png, p, width = 12, height = 5, dpi = 150)
ggsave(sub("\\.png$", ".pdf", out_png), p, width = 12, height = 5)
cat("\nSaved:", out_png, "\n")
