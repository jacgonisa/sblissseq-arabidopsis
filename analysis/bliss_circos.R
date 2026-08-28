#!/usr/bin/env Rscript
# Polished circular genome-wide DSB density (BiCroLab Plot 4 style), TAIR12.
# Region ring (centromere/pericentromere/NOR) + filled-area DSB tracks per sample, y-capped.
suppressMessages({ library(data.table); library(circlize) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
BRK  <- file.path(base,"results_TAIR12/breaks"); R <- file.path(base,"results_TAIR12/analysis/regions")
FIG  <- file.path(base,"results_TAIR12/analysis/figures")
samples <- c("BA1","BA2","old_BA1_BA2"); scol <- c(BA1="#1f77b4",BA2="#2ca02c",old_BA1_BA2="#d62728")
binw <- 5e5
chrsz <- fread(file.path(R,"genome_chr15.sizes"),header=FALSE); setnames(chrsz,c("chr","len")); chrsz<-chrsz[grepl("^Chr[1-5]$",chr)]
seqlens <- setNames(chrsz$len,chrsz$chr)
bins <- rbindlist(lapply(chrsz$chr,function(c){s<-seq(0,seqlens[c]-1,by=binw); data.table(chr=c,start=s,end=pmin(s+binw,seqlens[c]))}))
setkey(bins,chr,start,end)
dens <- lapply(samples,function(s){
  d<-fread(file.path(BRK,paste0(s,".break.bedgraph")),header=FALSE); setnames(d,c("chr","start","end","score")); d<-d[grepl("^Chr[1-5]$",chr)]
  tot<-sum(d$score); setkey(d,chr,start,end)
  ov<-foverlaps(d,bins,type="within",nomatch=0L); agg<-ov[,.(v=sum(score)),by=.(chr,start,end)]
  m<-merge(bins,agg,by=c("chr","start","end"),all.x=TRUE); m[is.na(v),v:=0][, v:=v*1e6/tot/(binw/1000)]; m })
names(dens)<-samples
# common y-cap = 95th pct across samples (so NOR clips, arms visible)
ycap <- quantile(unlist(lapply(dens,function(m)m$v)), 0.97)
rdr <- function(fn){r<-fread(file.path(R,fn),header=FALSE); r<-r[grepl("^Chr[1-5]$",V1)]; data.frame(chr=r$V1,start=r$V2,end=r$V3)}
cen<-rdr("centromere.bed"); peri<-rdr("pericentromere.bed"); nor<-rdr("nor_45s.bed")

png(file.path(FIG,"bicrolab_circos.png"), width=2100, height=2100, res=250)
circos.clear()
circos.par(start.degree=90, gap.degree=c(rep(1.5,nrow(chrsz)-1),6), track.margin=c(0.004,0.004))
circos.genomicInitialize(data.frame(chr=chrsz$chr,start=0,end=chrsz$len), plotType=c("axis","labels"), axis.labels.cex=0.5, labels.cex=0.9)
# region ring (robust: circos.rect per region, explicit colours)
regall <- rbind(data.frame(cen,col="#d62728",stringsAsFactors=FALSE),
                data.frame(peri,col="#ff7f0e",stringsAsFactors=FALSE),
                data.frame(nor,col="#9467bd",stringsAsFactors=FALSE))
circos.track(ylim=c(0,1), track.height=0.05, bg.border=NA, panel.fun=function(x,y){
  ch<-CELL_META$sector.index; sub<-regall[regall$chr==ch,]
  if(nrow(sub)) for(i in seq_len(nrow(sub))) circos.rect(sub$start[i],0,sub$end[i],1,col=sub$col[i],border=NA)
})
# per-sample filled-area density tracks
for(s in samples){
  circos.genomicTrack(as.data.frame(dens[[s]][,.(chr,start,end,v)]), numeric.column=4,
    ylim=c(0,ycap), track.height=0.17, bg.border="#eeeeee",
    panel.fun=function(region,value,...){
      v<-pmin(value[[1]],ycap)
      circos.genomicLines(region, data.frame(v), area=TRUE, col=scol[[s]], border=scol[[s]], lwd=0.6)
    })
  circos.yaxis(side="left", at=c(0,round(ycap)), labels.cex=0.4, sector.index="Chr1", track.index=get.current.track.index())
}
title("Genome-wide sBLISS DSB density (CPM/kb, 500 kb bins)", cex.main=1)
text(0,0.18,"ring: centromere (red) /\npericentromere (orange) / NOR (purple)",cex=0.55,col="#555555")
text(0,0,paste(c("inner→outer:",samples),collapse="\n"),cex=0.6)
text(0,-0.22,sprintf("(y capped at %.0f CPM/kb;\nNOR spikes clipped)",ycap),cex=0.5,col="#888888")
legend("bottomright", legend=samples, fill=scol[samples], cex=0.6, bty="n")
circos.clear(); dev.off()
cat("Saved bicrolab_circos.png  (ycap=",round(ycap,1)," CPM/kb)\n")
