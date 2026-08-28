#!/usr/bin/env Rscript
# DSB hotspot detection with BRGenomics (Hidmi et al. snippet, extended).
#
#   smp        <- import(bigwig)
#   smp        <- remove blacklist overlaps  (findOverlaps)   [optional]
#   tiles      <- tileGenome(chrom_sizes, tilewidth)
#   counts     <- getCountsByRegions(smp, tiles, field="score")
#   hotspots   <- tiles where counts exceed a Poisson significance threshold
#
# The blacklist is OPTIONAL. Pass "none" (or an empty file) to disable all
# blacklist filtering — centromeres and NORs are then analysed like any other
# region. Hotspots are instead ANNOTATED with their region (arm/peri/cen) and a
# NOR flag so repetitive-pileup artefacts remain identifiable downstream.
#
# Input track: pass the RAW per-base break track (counts), not a CPM track —
# the Poisson test below assumes integer counts. A .bedgraph or .bigwig is
# accepted (format auto-detected from the extension).
#
# Usage: Rscript hotspots_brgenomics.R <track.bedgraph|bigwig> <chrom.sizes> <blacklist|none> <tilewidth> <out_prefix>

options(scipen = 999)   # no scientific notation in coordinates
suppressMessages({
  library(rtracklayer)
  library(BRGenomics)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
bw_file    <- args[1]
chrsz_file <- args[2]
bl_file    <- args[3]
tilewidth  <- as.integer(args[4])
out_prefix <- args[5]

in_fmt <- if (grepl("\\.(bw|bigwig)$", bw_file, ignore.case = TRUE)) "bigwig" else "bedGraph"
cat("Importing", in_fmt, "track:", bw_file, "\n")
smp <- import(bw_file, format = in_fmt)
seqlevels(smp, pruning.mode = "coarse") <- grep("^Chr[1-5]$", seqlevels(smp), value = TRUE)

# ── blacklist filtering (paper snippet) — OPTIONAL ────────────────────────────
# A blacklist of "none"/"NA" (or an empty/missing file) disables all filtering.
use_blacklist <- !(is.na(bl_file) || tolower(bl_file) %in% c("none", "na", "")) &&
                 file.exists(bl_file) &&
                 file.info(bl_file)$size > 0
if (use_blacklist) {
  bls <- import(bl_file, format = "bed")
  ov  <- data.frame(findOverlaps(smp, bls))
  if (nrow(ov) > 0) {
    smp <- smp[-unique(ov$queryHits)]
    cat("Removed", length(unique(ov$queryHits)), "blacklisted intervals\n")
  }
} else {
  cat("No blacklist applied — centromeres and NORs analysed (flagged, not removed)\n")
}

# ── tile genome ───────────────────────────────────────────────────────────────
chrsz <- read.table(chrsz_file, col.names = c("chr", "len"))
seqlens <- setNames(chrsz$len, chrsz$chr)
tiles <- tileGenome(seqlens, tilewidth = tilewidth,
                    cut.last.tile.in.chrom = TRUE)

# ── count breaks per tile ─────────────────────────────────────────────────────
counts <- getCountsByRegions(smp, tiles, field = "score")
tiles$score <- counts

# drop blacklisted tiles entirely (only when a blacklist is in use)
if (use_blacklist) {
  bl_tiles <- unique(queryHits(findOverlaps(tiles, bls)))
  if (length(bl_tiles) > 0) tiles <- tiles[-bl_tiles]
}

# ── region / NOR annotation (so artefacts stay identifiable) ──────────────────
reg_dir <- Sys.getenv("BLISS_REGDIR", "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results/analysis/regions")
load_reg <- function(f) {
  p <- file.path(reg_dir, f)
  if (!file.exists(p) || file.info(p)$size == 0) return(GRanges())
  import(p, format = "bed")
}
cen  <- load_reg("centromere.bed")
peri <- load_reg("pericentromere.bed")
nor  <- load_reg("nor.bed")

region <- rep("arm", length(tiles))
if (length(peri)) region[queryHits(findOverlaps(tiles, peri))] <- "pericentromere"
if (length(cen))  region[queryHits(findOverlaps(tiles, cen))]  <- "centromere"
tiles$region <- region
tiles$is_NOR <- FALSE
if (length(nor)) tiles$is_NOR[queryHits(findOverlaps(tiles, nor))] <- TRUE

# ── Poisson hotspot calling ───────────────────────────────────────────────────
# H0: breaks distributed uniformly at the genome-wide mean rate per tile.
lambda <- mean(tiles$score)
tiles$pval  <- ppois(tiles$score - 1, lambda, lower.tail = FALSE)
tiles$padj  <- p.adjust(tiles$pval, method = "BH")
tiles$fold  <- tiles$score / lambda

hotspots <- tiles[tiles$padj < 0.001 & tiles$fold >= 3]
hotspots <- hotspots[order(-hotspots$score)]

cat(sprintf("Tilewidth: %d bp | mean breaks/tile (lambda): %.1f\n", tilewidth, lambda))
cat(sprintf("Total tiles tested: %d\n", length(tiles)))
cat(sprintf("Hotspots (padj<0.001 & fold>=3): %d\n", length(hotspots)))

# ── write outputs ─────────────────────────────────────────────────────────────
export_bed <- function(gr, file) {
  df <- data.frame(chr = as.character(seqnames(gr)),
                   start = start(gr) - 1, end = end(gr),
                   name = sprintf("hotspot_%d", seq_along(gr)),
                   score = round(gr$score, 2),
                   strand = ".",
                   fold = round(gr$fold, 2),
                   padj = signif(gr$padj, 3),
                   region = gr$region,
                   is_NOR = gr$is_NOR)
  write.table(df, file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}
export_bed(hotspots, paste0(out_prefix, ".hotspots.bed"))

# report how many hotspots fall in each region / are NOR artefacts
cat("Hotspots by region:\n"); print(table(hotspots$region))
cat(sprintf("Hotspots overlapping NOR 45S rDNA: %d\n", sum(hotspots$is_NOR)))

# full tile table (for genome-density plotting)
tile_df <- data.frame(chr = as.character(seqnames(tiles)),
                      start = start(tiles) - 1, end = end(tiles),
                      score = tiles$score, fold = tiles$fold,
                      padj = tiles$padj,
                      region = tiles$region, is_NOR = tiles$is_NOR)
write.table(tile_df, paste0(out_prefix, ".tiles.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Written:", paste0(out_prefix, ".hotspots.bed"), "\n")
cat("Written:", paste0(out_prefix, ".tiles.tsv"), "\n")
