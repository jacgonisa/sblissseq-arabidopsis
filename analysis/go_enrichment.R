suppressMessages({library(data.table); library(clusterProfiler); library(org.At.tair.db); library(ggplot2)})
A<-"/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis"
g<-fread(file.path(A,"go/gene_dsb_gdna.tsv"),header=FALSE)
setnames(g,c("chr","start","end","gene","bliss","gdna"))
g[, len:=(end-start)/1000][, bliss_kb:=bliss/len]
g<-g[bliss>=5 & gdna>=20]                      # require min coverage
g[, enr := (bliss/sum(bliss)) / (gdna/sum(gdna)) ]   # DSB-specific enrichment per gene
allg<-unique(g$gene)
top<-unique(g[order(-enr)][1:round(.05*.N)]$gene)     # top 5% by DSB/gDNA enrichment
cat("background genes:",length(allg)," top set:",length(top),"\n")
ego<-enrichGO(gene=top, universe=allg, OrgDb=org.At.tair.db, keyType="TAIR",
              ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1, readable=FALSE)
if(is.null(ego) || nrow(as.data.frame(ego))==0){ cat("No significant GO BP terms.\n") } else {
  df<-as.data.frame(ego); fwrite(df, file.path(A,"go/go_BP_topDSBgenes.csv"))
  cat("\n=== Top GO:BP terms (top-DSB protein-coding genes) ===\n")
  print(head(df[,c("ID","Description","GeneRatio","BgRatio","p.adjust","Count")],15))
   p<-barplot(ego, showCategory=15) + ggtitle("GO:BP — top-5% DSB-enriched protein-coding genes (gDNA-normalised)")
  ggsave(file.path(A,"figures/go_BP_topDSBgenes.png"), p, width=10, height=7, dpi=150)
  cat("Saved go_BP_topDSBgenes.png + csv\n")
}
