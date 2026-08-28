#!/usr/bin/env Rscript
# Chromosome-scale DSB density landscape (CPM-normalised, coarse tiles).
# Overlays centromere, pericentromere and NOR regions.
# Usage: Rscript genome_density.R <out_png>

suppressMessages({
  library(rtracklayer); library(BRGenomics); library(GenomicRanges)
  library(ggplot2); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
out_png <- ifelse(length(args) >= 1, args[1],
  "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/figures/genome_density.png")
tilew <- ifelse(length(args) >= 2, as.integer(args[2]), 100000L)
mask_nor <- length(args) >= 3 && tolower(args[3]) %in% c("masknor","mask","nor","true")
ymax <- suppressWarnings(as.numeric(Sys.getenv("BLISS_YMAX","")))   # fixed y-axis cap (shared across chr) if set
tile_lab <- if (tilew >= 1e6) sprintf("%g Mb", tilew/1e6) else
            if (tilew >= 1e3) sprintf("%g kb", tilew/1e3) else sprintf("%d bp", tilew)

base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
# Path overrides (default = ragtag; set these env vars for the TAIR12 run)
regdir  <- Sys.getenv("BLISS_REGDIR",   file.path(base, "results/analysis/regions"))
bwdir   <- Sys.getenv("BLISS_BWDIR",    file.path(base, "results/analysis/bigwig_cpm"))
chrfile <- Sys.getenv("BLISS_CHRSIZES", file.path(base, "reference/Col-0.chrom.sizes"))
chrsz <- read.table(chrfile)
chrsz <- chrsz[grepl("^Chr[1-5]$", chrsz$V1), ]
seqlens <- setNames(chrsz$V2, chrsz$V1)
tiles <- tileGenome(seqlens, tilewidth = tilew, cut.last.tile.in.chrom = TRUE)

samples <- c("BA1", "BA2", "old_BA1_BA2")
dens <- lapply(samples, function(s) {
  bw <- import(file.path(bwdir, sprintf("%s.cpm.bw", s)),
               format = "bigwig")
  cnt <- getCountsByRegions(bw, tiles, field = "score")
  data.frame(sample = s,
             chr = as.character(seqnames(tiles)),
             pos = (start(tiles) + end(tiles)) / 2 / 1e6,
             cpm = cnt / (tilew / 1000))   # CPM per kb
})
df <- bind_rows(dens)
df$sample <- factor(df$sample, levels = samples)

# By default NO masking — centromeres and NORs are shown (NOR spikes labelled).
# With mask_nor=TRUE, ONLY the NOR 45S rDNA bins are blanked (set NA) so the
# rest of the genome is readable; centromeres are still shown.
df$pos_bp <- df$pos * 1e6
if (mask_nor) {
  norbl <- read.table(file.path(regdir, "nor.bed"))
  names(norbl)[1:3] <- c("chr","start","end")
  for (i in seq_len(nrow(norbl))) {
    df$cpm[df$chr == norbl$chr[i] &
           df$pos_bp >= norbl$start[i] &
           df$pos_bp <= norbl$end[i]] <- NA
  }
}

# region overlays
read_reg <- function(f) {
  r <- read.table(file.path(regdir, f))
  r[grepl("^Chr[1-5]$", r$V1), 1:3] |> setNames(c("chr", "start", "end"))
}
cen  <- read_reg("centromere.bed");     cen$type  <- "Centromere"
peri <- read_reg("pericentromere.bed"); peri$type <- "Pericentromere"
nor_lab <- if (mask_nor) "NOR 45S rDNA (masked)" else "NOR 45S rDNA"
nor  <- read_reg("nor.bed");            nor$type  <- nor_lab
regions <- bind_rows(cen, peri, nor)
regions$xmin <- regions$start / 1e6
regions$xmax <- regions$end   / 1e6

p <- ggplot(df, aes(pos, cpm, colour = sample)) +
  geom_rect(data = regions, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = type),
            alpha = 0.15) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~chr, scales = if (!is.na(ymax)) "free_x" else "free", ncol = 1, strip.position = "right") +
  scale_fill_manual(values = setNames(c("red","orange","purple"),
                       c("Centromere","Pericentromere", nor_lab))) +
  scale_colour_manual(values = c("BA1" = "#1f77b4", "BA2" = "#2ca02c",
                                 "old_BA1_BA2" = "#d62728")) +
  labs(title = sprintf("sBLISS DSB density landscape (CPM per kb, %s windows)", tile_lab),
       subtitle = if (mask_nor)
         "NOR 45S rDNA bins MASKED (Chr2/Chr4 0-1 Mb) so arms/pericentromeres are visible. Shaded: centromere (red), pericentromere (orange), masked NOR (purple)."
       else
         "Shaded: centromere (red), pericentromere (orange), NOR 45S rDNA (purple). No blacklist - all regions shown; per-chromosome y-scale.",
       x = "Position (Mb)", y = "DSB CPM per kb",
       colour = "Sample", fill = "Region") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

if (!is.na(ymax)) p <- p + coord_cartesian(ylim = c(0, ymax))
ggsave(out_png, p, width = 11, height = 12, dpi = 150)
ggsave(sub("\\.png$", ".pdf", out_png), p, width = 11, height = 12)
cat("Saved:", out_png, "\n")
