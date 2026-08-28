#!/usr/bin/env Rscript
# Genome-wide BLISS DSB density vs WGA gDNA coverage (CPM/kb, 100 kb bins).
# Shows gDNA is ~uniform while BLISS spikes at the NOR -> why BLISS/gDNA reveals real rDNA enrichment.
suppressMessages({ library(rtracklayer); library(BRGenomics); library(GenomicRanges)
                   library(ggplot2); library(dplyr) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
A <- file.path(base,"results_TAIR12/analysis"); FIG <- file.path(A,"figures")
chrsz <- read.table(file.path(A,"regions/genome_chr15.sizes")); chrsz <- chrsz[grepl("^Chr[1-5]$",chrsz$V1),]
seqlens <- setNames(chrsz$V2,chrsz$V1); tilew <- 1e5
tiles <- tileGenome(seqlens, tilewidth=tilew, cut.last.tile.in.chrom=TRUE)
tracks <- list("sBLISS DSB (BA1)"=file.path(A,"bigwig_cpm/BA1.cpm.bw"),
               "WGA gDNA coverage"=file.path(A,"coverage_cpm/wga.cov.cpm.bw"))
tcol <- c("sBLISS DSB (BA1)"="#d62728","WGA gDNA coverage"="#555555")
df <- bind_rows(lapply(names(tracks), function(n){
  bw <- import(tracks[[n]], format="bigwig")
  data.frame(track=n, chr=as.character(seqnames(tiles)),
             pos=(start(tiles)+end(tiles))/2/1e6,
             cpmkb=getCountsByRegions(bw,tiles,field="score")/(tilew/1000))}))
df$track <- factor(df$track, levels=names(tracks))
rr <- function(f){r<-read.table(file.path(A,"regions",f));r<-r[grepl("^Chr[1-5]$",r$V1),];data.frame(chr=r$V1,xmin=r$V2/1e6,xmax=r$V3/1e6)}
cen<-rr("centromere.bed");cen$type<-"Centromere";peri<-rr("pericentromere.bed");peri$type<-"Pericentromere";nor<-rr("nor_45s.bed");nor$type<-"45S NOR"
reg<-bind_rows(cen,peri,nor)
p <- ggplot(df, aes(pos,cpmkb,colour=track)) +
  geom_rect(data=reg, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=type), alpha=0.13) +
  geom_line(linewidth=0.45) +
  facet_wrap(~chr, scales="free", ncol=1, strip.position="right") +
  scale_fill_manual(values=c("Centromere"="red","Pericentromere"="orange","45S NOR"="purple"),name="Region") +
  scale_colour_manual(values=tcol, name=NULL) +
  labs(title="Genome-wide: sBLISS DSB density vs WGA gDNA coverage (100 kb, CPM/kb)",
       subtitle="WGA gDNA is ~uniform; sBLISS spikes at the 45S NOR -> the rDNA enrichment is DSB-specific, not DNA-content. Per-chr y-scale.",
       x="Position (Mb)", y="CPM per kb") +
  theme_bw(base_size=10)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
ggsave(file.path(FIG,"genome_bliss_vs_wga.png"), p, width=11, height=12, dpi=150)
cat("Saved genome_bliss_vs_wga.png\n")
