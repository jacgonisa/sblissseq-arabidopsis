suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
M<-"/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis/metaplots"
FIG<-"/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis/figures"
read_prof<-function(tab){
  tb<-readLines(tab); rows<-tb[grepl("^BA1|^wga",tb)]
  out<-rbindlist(lapply(rows,function(l){v<-strsplit(l,"\t")[[1]]
    data.table(track=ifelse(grepl("BA1",v[1]),"sBLISS DSB","WGA gDNA"),
               i=seq_along(v[3:length(v)]), val=as.numeric(v[3:length(v)]))}))
  out[, v:=val/mean(val), by=track]; out }
panel<-function(tab,title,kind){
  d<-read_prof(tab); n<-max(d$i)
  if(kind=="tss"){ d[, x:=(i-1)/(n-1)*4000-2000]; vl<-0; xlab<-"bp from TSS"; brk<-c(-2000,0,2000); lbl<-c("-2kb","TSS","+2kb") }
  else { d[, x:=i]; up<-round(n*2000/7000); bo<-round(n*5000/7000); vl<-c(up,bo); xlab<-NULL; brk<-c(up,bo); lbl<-c("start","end") }
  pk<-d[track=="sBLISS DSB"][which.max(v)]; pkg<-d[track=="WGA gDNA"][i==pk$i]
  ggplot(d,aes(x,v,colour=track))+geom_line(linewidth=1.1)+
    {if(kind=="tss") geom_vline(xintercept=vl,linetype=2,colour="grey50") else geom_vline(xintercept=vl,linetype=2,colour="grey50")}+
    scale_colour_manual(values=c("sBLISS DSB"="#d62728","WGA gDNA"="#555555"),name=NULL)+
    scale_x_continuous(breaks=brk,labels=lbl)+
    labs(title=sprintf("%s (sBLISS peak %.2f× vs gDNA %.2f×)",title,pk$v,pkg$v),x=xlab,y="signal / mean")+
    theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold",size=10))
}
pA<-panel(file.path(M,"TSS_prof.tab"),"TSS","tss")
pB<-panel(file.path(M,"genebody_prof.tab"),"Gene body","scale")
pC<-panel(file.path(M,"TE_prof.tab"),"Transposons","scale")
ggsave(file.path(FIG,"metaplot_shape_bliss_vs_gdna.png"), pA|pB|pC, width=13, height=4.6, dpi=140)
cat("Saved metaplot_shape_bliss_vs_gdna.png\n")
