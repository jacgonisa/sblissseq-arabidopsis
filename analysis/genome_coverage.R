#!/usr/bin/env Rscript
# Simple genome-wide READ COVERAGE (CPM) landscape — actual read depth, not 5' breaks.
# Usage: Rscript genome_coverage.R <out_png> [tilewidth]
suppressMessages({ library(rtracklayer); library(BRGenomics); library(GenomicRanges)
                   library(ggplot2); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE)
out_png <- ifelse(length(args)>=1, args[1],
  "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis/figures/genome_coverage_100kb.png")
tilew <- ifelse(length(args)>=2, as.integer(args[2]), 100000L)
tlab  <- if(tilew>=1e6) sprintf("%g Mb",tilew/1e6) else sprintf("%g kb",tilew/1e3)

base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
regdir  <- Sys.getenv("BLISS_REGDIR",   file.path(base,"results_TAIR12/analysis/regions"))
bwdir   <- Sys.getenv("BLISS_BWDIR",    file.path(base,"results_TAIR12/analysis/coverage_cpm"))
chrfile <- Sys.getenv("BLISS_CHRSIZES", file.path(base,"results_TAIR12/analysis/regions/genome_chr15.sizes"))
suffix  <- Sys.getenv("BLISS_BWSUFFIX", ".cov.cpm.bw")
qtag    <- if (grepl("q10", suffix)) " — MAPQ>=10 (uniquely-mapped reads only)" else " — all reads"

chrsz <- read.table(chrfile); chrsz <- chrsz[grepl("^Chr[1-5]$",chrsz$V1),]
seqlens <- setNames(chrsz$V2, chrsz$V1)
tiles <- tileGenome(seqlens, tilewidth=tilew, cut.last.tile.in.chrom=TRUE)
samples <- c("BA1","BA2","old_BA1_BA2")
df <- bind_rows(lapply(samples, function(s){
  bw <- import(file.path(bwdir, paste0(s,suffix)), format="bigwig")
  cnt <- getCountsByRegions(bw, tiles, field="score")
  data.frame(sample=s, chr=as.character(seqnames(tiles)),
             pos=(start(tiles)+end(tiles))/2/1e6, cov=cnt/(tilew/1000))
}))
df$sample <- factor(df$sample, levels=samples)

read_reg <- function(f){ r<-read.table(file.path(regdir,f)); r[grepl("^Chr[1-5]$",r$V1),1:3] |> setNames(c("chr","start","end")) }
cen<-read_reg("centromere.bed"); cen$type<-"Centromere"
peri<-read_reg("pericentromere.bed"); peri$type<-"Pericentromere"
nor<-read_reg("nor.bed"); nor$type<-"45S rDNA (NOR)"
regions<-bind_rows(cen,peri,nor); regions$xmin<-regions$start/1e6; regions$xmax<-regions$end/1e6

p <- ggplot(df, aes(pos, cov, colour=sample)) +
  geom_rect(data=regions, inherit.aes=FALSE,
            aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=type), alpha=0.15) +
  geom_line(linewidth=0.4) +
  facet_wrap(~chr, scales="free", ncol=1, strip.position="right") +
  scale_fill_manual(values=c("Centromere"="red","Pericentromere"="orange","45S rDNA (NOR)"="purple")) +
  scale_colour_manual(values=c(BA1="#1f77b4",BA2="#2ca02c",old_BA1_BA2="#d62728")) +
  labs(title=sprintf("TAIR12: genome-wide READ COVERAGE (CPM per kb, %s tiles)%s", tlab, qtag),
       subtitle="Read depth (dedup BAM, CPM). Spikes = read pile-ups (rDNA/repeats); dips = low/unmappable coverage. Per-chromosome y-scale.",
       x="Position (Mb)", y="Read coverage (CPM/kb)", colour="Sample", fill="Region") +
  theme_bw(base_size=10)+theme(legend.position="bottom", plot.title=element_text(face="bold"))
ggsave(out_png, p, width=11, height=12, dpi=150)
cat("Saved:", out_png, "\n")
