#!/usr/bin/env Rscript
# Hotspot calling against the WGA gDNA baseline (Poisson, lambda from gDNA per tile,
# not flat). Corrects mappability + copy-number -> should shrink rDNA dominance.
suppressMessages({ library(data.table) })
base<-"/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"; A<-file.path(base,"results_TAIR12/analysis")
tiles<-fread(file.path(A,"hotspots/gdna_tiles.bedcov"),header=FALSE); setnames(tiles,c("chr","start","end","gbases"))
tiles<-tiles[grepl("^Chr[1-5]$",chr)]; setkey(tiles,chr,start,end)
gtot<-sum(as.numeric(tiles$gbases))
rd<-function(fn){r<-fread(file.path(A,"regions",fn),header=FALSE)[grepl("^Chr[1-5]$",V1)];setkey(r,V1,V2,V3);r}
cen<-rd("centromere.bed");peri<-rd("pericentromere.bed");nor<-rd("nor_45s.bed")
inreg<-function(dt,reg){ o<-foverlaps(dt,reg,by.x=c("chr","start","end"),type="any",nomatch=0L,which=TRUE); idx<-rep(FALSE,nrow(dt)); idx[unique(o$xid)]<-TRUE; idx }
samples<-c("BA1","BA2","old_BA1_BA2"); hot<-list()
for(s in samples){
  bg<-fread(file.path(base,sprintf("results_TAIR12/breaks/%s.break.bedgraph",s)),header=FALSE); setnames(bg,c("chr","start","end","score")); bg<-bg[grepl("^Chr[1-5]$",chr)]; setkey(bg,chr,start,end)
  ov<-foverlaps(bg,tiles,type="within",nomatch=0L); agg<-ov[,.(b=sum(score)),by=.(chr,start,end)]
  t<-merge(tiles,agg,by=c("chr","start","end"),all.x=TRUE); t[is.na(b),b:=0]
  btot<-sum(as.numeric(t$b))
  t[, exp := pmax(gbases*(btot/gtot), 1)]          # expected from gDNA (>=1 pseudocount)
  t[, fold := b/exp]
  t[, p := ppois(b-1, exp, lower.tail=FALSE)]
  t[, padj := p.adjust(p, "BH")]
  h<-t[padj<0.001 & fold>=3 & b>=10]               # min 10 breaks to avoid low-count noise
  h[, key:=paste(chr,start)]
  hot[[s]]<-h
  cat(sprintf("%-12s gDNA-baseline hotspots: %d\n", s, nrow(h)))
}
robust<-Reduce(intersect, lapply(hot,function(h)h$key))
rb<-hot[["BA1"]][key %in% robust]
rb[, is_NOR:=inreg(rb,nor)][, is_cen:=inreg(rb,cen)][, is_peri:=inreg(rb,peri)]
rb[, cat:=ifelse(is_NOR,"45S rDNA",ifelse(is_cen,"Centromere",ifelse(is_peri,"Pericentromere","Arm")))]
cat("\n=== gDNA-baseline ROBUST hotspots (3-way):", length(robust), "===\n")
print(rb[,.N,by=cat][order(-N)])
cat("\nvs flat-Poisson robust (6322: 5856 rDNA / 417 peri / 34 cen / 15 arm)\n")
fwrite(rb[,.(chr,start,end,b,fold,padj,cat)], file.path(A,"hotspots/robust_hotspots_vs_gdna.csv"))
cat("Wrote robust_hotspots_vs_gdna.csv\n")
