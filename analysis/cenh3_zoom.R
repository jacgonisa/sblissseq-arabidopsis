#!/usr/bin/env Rscript
# Centromere zoom: sBLISS DSB density vs CENH3 ChIP-seq, per chromosome.
# Both at 10 kb on the same Col-0 ragtag assembly. Two stacked tracks per chr.
# Usage: Rscript cenh3_zoom.R [out_png]
options(scipen = 999)
suppressMessages({
  library(rtracklayer); library(BRGenomics); library(GenomicRanges)
  library(ggplot2); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
out_png <- ifelse(length(args) >= 1, args[1],
  file.path(base, "results/analysis/figures/cenh3_zoom.png"))
cenh3_file <- "/home/jg2070/Desktop/PhD/crossover/only_kmer_pipeline/10-plotting_crossover/genome_and_annotation/CENH3.Col.ChIP.10kb.Jan2023.trimmed.bedgraph"
FLANK <- 1.5e6   # zoom = centromere core +/- 1.5 Mb

chrs <- paste0("Chr", 1:5)
chrsz <- read.table(file.path(base, "reference/Col-0.chrom.sizes"))
chrsz <- chrsz[chrsz$V1 %in% chrs, ]
seqlens <- setNames(chrsz$V2, chrsz$V1)

cen <- read.table(file.path(base, "results/analysis/regions/centromere.bed"),
                  col.names = c("chr","start","end"))
cen <- cen[cen$chr %in% chrs, ]
zoom <- cen %>% mutate(zlo = pmax(0, start - FLANK),
                       zhi = pmin(seqlens[chr], end + FLANK))

# ── sBLISS: 10 kb CPM tiles ───────────────────────────────────────────────────
tiles <- tileGenome(seqlens, tilewidth = 10000, cut.last.tile.in.chrom = TRUE)
samples <- c("BA1","BA2","old_BA1_BA2")
sb <- lapply(samples, function(s) {
  bw <- import(file.path(base, sprintf("results/analysis/bigwig_cpm/%s.cpm.bw", s)),
               format = "bigwig")
  data.frame(chr = as.character(seqnames(tiles)),
             mid = (start(tiles)+end(tiles))/2,
             value = getCountsByRegions(bw, tiles, field = "score")/10,  # CPM per kb
             series = s, track = "sBLISS DSB density (CPM/kb)")
}) %>% bind_rows()

# ── CENH3 ChIP (already 10 kb, ragtag) ────────────────────────────────────────
ce <- read.table(cenh3_file, header = FALSE)
ce$chr <- sub(":.*", "", ce$V1)
ce <- ce %>% filter(chr %in% chrs) %>%
  transmute(chr, mid = (V2+V3)/2, value = V4,
            series = "CENH3", track = "CENH3 ChIP (norm.)")

df <- bind_rows(sb, ce)
df$track <- factor(df$track, levels = c("sBLISS DSB density (CPM/kb)", "CENH3 ChIP (norm.)"))

# clip to zoom window per chromosome
df <- df %>% inner_join(zoom[,c("chr","zlo","zhi")], by = "chr") %>%
  filter(mid >= zlo & mid <= zhi)
df$mid_Mb <- df$mid/1e6
cen$xmin <- cen$start/1e6; cen$xmax <- cen$end/1e6

cols <- c(BA1="#1f77b4", BA2="#2ca02c", old_BA1_BA2="#d62728", CENH3="black")

p <- ggplot(df, aes(mid_Mb, value, colour = series)) +
  geom_rect(data = cen, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "purple", alpha = 0.10) +
  geom_line(linewidth = 0.5) +
  facet_grid(track ~ chr, scales = "free", switch = "y") +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Centromere zoom: sBLISS DSB density vs CENH3 ChIP-seq",
       subtitle = "Per chromosome, centromere core +/- 1.5 Mb. Shaded = CENH3 centromere core (region BED).",
       x = "Position (Mb)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        strip.placement = "outside", strip.background.y = element_blank())

ggsave(out_png, p, width = 14, height = 6, dpi = 150)
ggsave(sub("\\.png$",".pdf",out_png), p, width = 14, height = 6)
cat("Saved:", out_png, "\n")

# correlation sBLISS(BA1) vs CENH3 within zoom (10 kb bins)
# bin key (mid rounded to 10 kb) to bridge 1-based tile vs 0-based bedgraph offset
m <- df %>% filter(series %in% c("BA1","CENH3")) %>%
  mutate(bin = paste(chr, floor(mid/1e4))) %>%
  select(bin, series, value) %>% tidyr::pivot_wider(names_from=series, values_from=value)
cat(sprintf("Pearson r (BA1 sBLISS vs CENH3, zoom 10kb bins): %.3f\n",
            cor(m$BA1, m$CENH3, use = "complete.obs")))
cat(sprintf("Spearman r: %.3f\n",
            cor(m$BA1, m$CENH3, method = "spearman", use = "complete.obs")))
