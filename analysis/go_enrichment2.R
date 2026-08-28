suppressMessages({library(data.table); library(clusterProfiler); library(org.At.tair.db)})
A<-"/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis"
g<-fread(file.path(A,"go/gene_dsb_gdna.tsv"),header=FALSE)
setnames(g,c("chr","start","end","gene","bliss","gdna"))
g[, len:=(end-start)/1000][, bliss_kb:=bliss/len]
g<-g[bliss>=5 & gdna>=20]
g[, enr := (bliss/sum(bliss)) / (gdna/sum(gdna)) ]
allg<-unique(g$gene)
run<-function(metric, ont, frac=0.05){
  top<-unique(g[order(-get(metric))][1:round(frac*.N)]$gene)
  e<-tryCatch(enrichGO(gene=top, universe=allg, OrgDb=org.At.tair.db, keyType="TAIR",
     ont=ont, pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=FALSE),error=function(x)NULL)
  df<-if(is.null(e))data.frame() else as.data.frame(e)
  cat(sprintf("\n## metric=%s ont=%s top=%d -> %d sig terms\n",metric,ont,length(top),nrow(df)))
  if(nrow(df)) print(head(df[,c("Description","GeneRatio","BgRatio","p.adjust","Count")],10))
  if(nrow(df)) fwrite(df, file.path(A,sprintf("go/go_%s_%s.csv",metric,ont)))
  invisible(df)
}
for(m in c("enr","bliss_kb")) for(o in c("BP","MF","CC")) run(m,o)
