#!/usr/bin/env Rscript
# BiCroLab blissNPanalysis-style plots, adapted for Arabidopsis (TAIR12).
#  Plot 1: DSB-location count vs UMI threshold (1..10)  [library complexity]
#  Plot 4: circular genome-wide DSB density (1 Mb bins, z-capped at mean+3SD)
suppressMessages({ library(data.table); library(ggplot2); library(circlize) })

base   <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
BRK    <- file.path(base, "results_TAIR12/breaks")
CHR    <- file.path(base, "results_TAIR12/analysis/regions/genome_chr15.sizes")
FIG    <- file.path(base, "results_TAIR12/analysis/figures")
samples<- c("BA1","BA2","old_BA1_BA2")
scol   <- c(BA1="#1f77b4", BA2="#2ca02c", old_BA1_BA2="#d62728")

chrsz  <- fread(CHR, header=FALSE); setnames(chrsz, c("chr","len"))
chrsz  <- chrsz[grepl("^Chr[1-5]$", chr)]

# load break bedgraphs (chr,start,end,score=UMIs per location), Chr1-5
bg <- lapply(samples, function(s){
  d <- fread(file.path(BRK, paste0(s,".break.bedgraph")), header=FALSE)
  setnames(d, c("chr","start","end","score")); d[grepl("^Chr[1-5]$", chr)]
}); names(bg) <- samples

# ── Plot 1: DSB locations with >= t UMIs (t=1..10) ───────────────────────────
thr <- 1:10
d1 <- rbindlist(lapply(samples, function(s)
  data.table(sample=s, threshold=thr,
             n=sapply(thr, function(t) sum(bg[[s]]$score >= t)))))
d1$sample <- factor(d1$sample, levels=samples)
p1 <- ggplot(d1, aes(threshold, n, colour=sample)) +
  geom_line(linewidth=1) + geom_point(size=1.8) +
  scale_y_log10(labels=scales::comma) +
  scale_x_continuous(breaks=thr) +
  scale_colour_manual(values=scol, name=NULL) +
  labs(title="DSB locations vs UMI threshold (BiCroLab Plot 1)",
       subtitle="Number of distinct break sites with at least N UMIs — library complexity / saturation",
       x="UMIs per location (>=)", y="DSB locations (log scale)") +
  theme_bw(base_size=12) + theme(legend.position="bottom",
       plot.title=element_text(face="bold"), panel.grid.minor=element_blank())
ggsave(file.path(FIG,"bicrolab_umi_threshold.png"), p1, width=8, height=5.5, dpi=150)
fwrite(dcast(d1, threshold~sample, value.var="n"),
       file.path(base,"results_TAIR12/analysis/bicrolab_umi_threshold.csv"))
cat("Saved bicrolab_umi_threshold.png + .csv\n"); print(dcast(d1, threshold~sample, value.var="n"))

# ── 2 kb DSB-density distribution (genome-wide) ──────────────────────────────
seqlens <- setNames(chrsz$len, chrsz$chr)
binw2 <- 2000
bins2 <- rbindlist(lapply(chrsz$chr, function(c){
  L<-seqlens[c]; s<-seq(0,L-1,by=binw2); data.table(chr=c,start=s,end=pmin(s+binw2,L))}))
setkey(bins2, chr, start, end)
dist2 <- rbindlist(lapply(samples, function(s){
  d<-copy(bg[[s]]); setkey(d,chr,start,end)
  ov<-foverlaps(d,bins2,type="within",nomatch=0L)
  agg<-ov[,.(v=sum(score)),by=.(chr,start,end)]
  m<-merge(bins2,agg,by=c("chr","start","end"),all.x=TRUE); m[is.na(v),v:=0]
  data.table(sample=s, v=m$v)}))
dist2$sample<-factor(dist2$sample,levels=samples)
pd <- ggplot(dist2[v>0], aes(v, colour=sample)) +
  stat_ecdf(linewidth=1) + scale_x_log10(labels=scales::comma) +
  scale_colour_manual(values=scol,name=NULL) +
  labs(title="DSBs per 2 kb window — distribution (mouse-tutorial window)",
       subtitle="ECDF of DSB UMIs per non-empty 2 kb window; long right tail = hot windows (centromere/NOR)",
       x="DSB UMIs per 2 kb window (log)", y="cumulative fraction of windows") +
  theme_bw(base_size=12)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
ggsave(file.path(FIG,"bicrolab_2kb_distribution.png"), pd, width=8, height=5.5, dpi=150)
cat("Saved bicrolab_2kb_distribution.png\n")
cat("2kb windows: median/95th/99th/max DSB UMIs per sample:\n")
print(dist2[v>0, .(median=median(v), p95=quantile(v,.95), p99=quantile(v,.99), max=max(v)), by=sample])

# ── Plot 4: circular genome-wide DSB density (1 Mb, z-capped) ────────────────
binw <- 1e6
seqlens <- setNames(chrsz$len, chrsz$chr)
# build 1Mb bins per chr
bins <- rbindlist(lapply(chrsz$chr, function(c){
  L <- seqlens[c]; starts <- seq(0, L-1, by=binw)
  data.table(chr=c, start=starts, end=pmin(starts+binw, L))
}))
# count score per bin per sample (foverlaps)
countbin <- function(s){
  d <- bg[[s]]; setkey(d, chr, start, end)
  b <- copy(bins); setkey(b, chr, start, end)
  ov <- foverlaps(d, b, type="within", nomatch=0L)
  agg <- ov[, .(v=sum(score)), by=.(chr,start,end)]
  m <- merge(bins, agg, by=c("chr","start","end"), all.x=TRUE); m[is.na(v), v:=0]
  # z-cap at mean+3SD (BiCroLab)
  M <- mean(m$v); SD <- sd(m$v); cap <- M+3*SD
  m[v>cap, v:=cap]; m
}
cdata <- lapply(samples, countbin); names(cdata) <- samples

png(file.path(FIG,"bicrolab_circos.png"), width=2000, height=2000, res=240)
circos.clear()
circos.par(start.degree=90, gap.degree=c(rep(2,nrow(chrsz)-1),8))
gdf <- data.frame(chr=chrsz$chr, start=0, end=chrsz$len)
circos.genomicInitialize(gdf, plotType=c("axis","labels"))
for (s in samples){
  m <- cdata[[s]]
  circos.genomicTrack(as.data.frame(m[,.(chr,start,end,v)]), numeric.column=4,
    track.height=0.13, panel.fun=function(region,value,...){
      circos.genomicLines(region, value, type="h", col=scol[[s]], lwd=1.2)
    })
  circos.text(0, 0.5, s, sector.index=chrsz$chr[1], facing="bending.inside",
              cex=0.6, col=scol[[s]], adj=c(1.2,0.5))
}
title("Genome-wide DSB density (BiCroLab Plot 4) — 1 Mb bins, z-capped mean+3SD")
circos.clear(); dev.off()
cat("Saved bicrolab_circos.png\n")
