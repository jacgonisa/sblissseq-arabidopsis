#!/usr/bin/env Rscript
# BA1 minimap2 genome-density landscape at 3 read-processing levels:
#  raw (all primary, incl MAPQ0, no dedup) | dedup (incl MAPQ0) | dedup+MAPQ>=1.
# Shows how multimapper handling changes the centromere/NOR signal.
suppressMessages({ library(data.table); library(ggplot2) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
BRK  <- file.path(base,"results_TAIR12/breaks")
CHR  <- file.path(base,"results_TAIR12/analysis/regions/genome_chr15.sizes")
R    <- file.path(base,"results_TAIR12/analysis/regions")
FIG  <- file.path(base,"results_TAIR12/analysis/figures")
tilew <- 1e5
levels_ <- c(raw="BA1.mm2raw.break.bedgraph", dedup="BA1.mm2dedup.break.bedgraph",
             "dedup+MAPQ>=1"="BA1.mm2q1.break.bedgraph")
lcol <- c(raw="#888888", dedup="#1f77b4", "dedup+MAPQ>=1"="#d62728")

chrsz <- fread(CHR,header=FALSE); setnames(chrsz,c("chr","len")); chrsz <- chrsz[grepl("^Chr[1-5]$",chr)]
df <- rbindlist(lapply(names(levels_), function(L){
  d <- fread(file.path(BRK, levels_[[L]]), header=FALSE); setnames(d,c("chr","start","end","score"))
  d <- d[grepl("^Chr[1-5]$",chr)]
  tot <- sum(d$score)
  d[, tile := floor(start/tilew)]
  agg <- d[, .(cpmkb=sum(score)*1e6/tot/(tilew/1000)), by=.(chr,tile)]
  agg[, pos := (tile*tilew + tilew/2)/1e6][, level := L][]
}))
df$level <- factor(df$level, levels=names(levels_))

read_reg <- function(f){ r<-fread(file.path(R,f),header=FALSE); r<-r[grepl("^Chr[1-5]$",V1)]; data.table(chr=r$V1,xmin=r$V2/1e6,xmax=r$V3/1e6)}
cen<-read_reg("centromere.bed");cen$type<-"Centromere"
peri<-read_reg("pericentromere.bed");peri$type<-"Pericentromere"
nor<-read_reg("nor_45s.bed");nor$type<-"45S NOR"
reg<-rbindlist(list(cen,peri,nor))

p <- ggplot(df, aes(pos,cpmkb,colour=level)) +
  geom_rect(data=reg, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=type), alpha=0.13) +
  geom_line(linewidth=0.45) +
  facet_wrap(~chr, scales="free", ncol=1, strip.position="right") +
  scale_fill_manual(values=c("Centromere"="red","Pericentromere"="orange","45S NOR"="purple"),name="Region") +
  scale_colour_manual(values=lcol, name="Read processing") +
  labs(title="BA1 minimap2 DSB density — effect of multimapper handling (100 kb, CPM/kb)",
       subtitle="raw (all primary, incl MAPQ0) vs dedup (incl MAPQ0) vs dedup+MAPQ>=1 (multimappers dropped). Per-chr y-scale.",
       x="Position (Mb)", y="DSB CPM per kb") +
  theme_bw(base_size=10)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
ggsave(file.path(FIG,"genome_density_mm2_processing.png"), p, width=11, height=12, dpi=150)
cat("Saved genome_density_mm2_processing.png\n")
# summary: centromere+NOR fraction of total signal at each level
for(L in names(levels_)){
  d<-fread(file.path(BRK,levels_[[L]]),header=FALSE); setnames(d,c("chr","start","end","score")); d<-d[grepl("^Chr[1-5]$",chr)]
  tot<-sum(d$score)
  cat(sprintf("%-16s total events=%d\n", L, tot))
}
