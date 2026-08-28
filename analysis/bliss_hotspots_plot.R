#!/usr/bin/env Rscript
# Genome-wide map of robust DSB hotspots (TAIR12), coloured by category (rDNA / cen / peri / arm).
suppressMessages({ library(data.table); library(ggplot2) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"; A <- file.path(base,"results_TAIR12/analysis")
h <- fread(file.path(A,"hotspots/robust_hotspots.annotated.tsv"))
h[, cat := ifelse(is_NOR=="TRUE","45S rDNA (NOR)",
            ifelse(region=="centromere","Centromere", ifelse(region=="pericentromere","Pericentromere","Arm")))]
# tag 5S overlaps
fives <- fread(file.path(A,"regions/rdna_5s.bed"),header=FALSE)
for(i in 1:nrow(fives)) h[chr==fives$V1[i] & start>=fives$V2[i] & end<=fives$V3[i], cat:="5S rDNA"]
h$cat <- factor(h$cat, levels=c("45S rDNA (NOR)","5S rDNA","Pericentromere","Centromere","Arm"))
h$mb <- (h$start+h$end)/2/1e6
chrsz <- fread(file.path(A,"regions/genome_chr15.sizes"),header=FALSE)[grepl("^Chr[1-5]$",V1)]
rdr <- function(fn){r<-fread(file.path(A,"regions",fn),header=FALSE)[grepl("^Chr[1-5]$",V1)];data.frame(chr=r$V1,xmin=r$V2/1e6,xmax=r$V3/1e6)}
cen<-rdr("centromere.bed");cen$type<-"Centromere";peri<-rdr("pericentromere.bed");peri$type<-"Pericentromere";nor<-rdr("nor_45s.bed");nor$type<-"45S NOR"
reg<-rbind(cen,peri,nor)
cols <- c("45S rDNA (NOR)"="#d62728","5S rDNA"="#17becf","Pericentromere"="#ff7f0e","Centromere"="#9467bd","Arm"="#1f77b4")
# counts for subtitle (non-rDNA)
n_tot<-nrow(h); n_nor<-sum(h$cat=="45S rDNA (NOR)"); n_non<-n_tot-n_nor
sub<-sprintf("%d robust 1 kb hotspots: %d are 45S rDNA; of the %d non-rDNA, %d%% pericentromeric, %d centromere, %d arm",
  n_tot, n_nor, n_non, round(100*sum(h$cat=="Pericentromere")/n_non), sum(h$cat=="Centromere"), sum(h$cat=="Arm"))
p <- ggplot(h, aes(mb, score, colour=cat)) +
  geom_rect(data=reg, inherit.aes=FALSE, aes(xmin=xmin,xmax=xmax,ymin=1,ymax=Inf,fill=type), alpha=0.12) +
  geom_point(size=0.8, alpha=0.7) +
  facet_wrap(~chr, scales="free_x", ncol=1, strip.position="right") +
  scale_y_log10() +
  scale_colour_manual(values=cols, name="hotspot class") +
  scale_fill_manual(values=c("Centromere"="#9467bd","Pericentromere"="#ff7f0e","45S NOR"="#d62728"), name="region", guide="none") +
  labs(title="Genome-wide robust DSB hotspots (TAIR12, 1 kb, 3-way reproducible)",
       subtitle=sub, x="Position (Mb)", y="hotspot DSB score (log10)") +
  theme_bw(base_size=10)+theme(legend.position="bottom", plot.title=element_text(face="bold"))
ggsave(file.path(A,"figures/hotspots_genomewide.png"), p, width=12, height=11, dpi=150)
cat("Saved hotspots_genomewide.png\n")
print(h[, .N, by=cat][order(-N)])
