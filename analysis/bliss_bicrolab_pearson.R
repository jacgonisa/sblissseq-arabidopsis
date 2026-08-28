#!/usr/bin/env Rscript
# BiCroLab Plot 2 — inter-sample Pearson correlation of DSB counts in 100 kb bins.
# Score-weighted bin counts (UMIs), all non-empty bins; Pearson; lower-triangle heatmap.
suppressMessages({ library(data.table); library(ggplot2) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
BRK  <- file.path(base,"results_TAIR12/breaks")
CHR  <- file.path(base,"results_TAIR12/analysis/regions/genome_chr15.sizes")
FIG  <- file.path(base,"results_TAIR12/analysis/figures")
samples <- c("BA1","BA2","old_BA1_BA2")
BINW <- as.numeric(Sys.getenv("BLISS_BINW", "1e5"))
TAG  <- Sys.getenv("BLISS_BINTAG", "100kb")

chrsz <- fread(CHR,header=FALSE); setnames(chrsz,c("chr","len")); chrsz <- chrsz[grepl("^Chr[1-5]$",chr)]
bins <- rbindlist(lapply(1:nrow(chrsz), function(i){
  L<-chrsz$len[i]; s<-seq(0,L-1,by=BINW); data.table(chr=chrsz$chr[i],start=s,end=pmin(s+BINW,L))}))
bins[, bin:=.I]; setkey(bins, chr, start, end)

mat <- sapply(samples, function(s){
  d <- fread(file.path(BRK,paste0(s,".break.bedgraph")),header=FALSE)
  setnames(d,c("chr","start","end","score")); d <- d[grepl("^Chr[1-5]$",chr)]; setkey(d,chr,start,end)
  ov <- foverlaps(d,bins,type="within",nomatch=0L)
  agg <- ov[, .(v=sum(score)), by=bin]
  out <- rep(0, nrow(bins)); out[agg$bin] <- agg$v; out
})
mat <- mat[rowSums(mat)>0, ]           # non-empty bins (BiCroLab)
cat("non-empty 100kb bins:", nrow(mat), "\n")
P <- cor(mat, method="pearson"); S <- cor(mat, method="spearman")
cat("Pearson:\n"); print(round(P,3)); cat("Spearman:\n"); print(round(S,3))

# heatmap (Pearson, full matrix with values)
m <- as.data.table(as.table(P)); setnames(m,c("x","y","r"))
m$x <- factor(m$x,levels=samples); m$y <- factor(m$y,levels=rev(samples))
p <- ggplot(m, aes(x,y,fill=r)) + geom_tile(colour="white") +
  geom_text(aes(label=sprintf("%.3f",r)), size=4.2) +
  scale_fill_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B",midpoint=0.9,limits=c(0.6,1),name="Pearson r") +
  labs(title=sprintf("Inter-sample DSB correlation (BiCroLab Plot 2) — %s bins", TAG),
       subtitle=sprintf("Score-weighted DSB UMIs in %s bins (%s non-empty), Pearson", TAG, format(nrow(mat),big.mark=",")),
       x=NULL,y=NULL) +
  coord_equal() + theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"), panel.grid=element_blank())
ggsave(file.path(FIG, sprintf("bicrolab_pearson_heatmap_%s.png",TAG)), p, width=6.5, height=5.5, dpi=150)
cat("Saved bicrolab_pearson_heatmap_",TAG,".png\n")
