#!/usr/bin/env Rscript
# MASTER comparison (BA1, minimap2): compartment x method.
# Methods cross MAPQ handling (>=0 keep multimappers vs >=1 drop) and normalisation
# (per Mb of DNA / per mapped read / divided by WGA gDNA).
suppressMessages({ library(data.table); library(ggplot2) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"; FIG <- file.path(base,"results_TAIR12/analysis/figures")
q1 <- fread(file.path(base,"results_TAIR12/analysis/compartment_enrichment_mm2q1.csv"))   # MAPQ>=1
gd <- fread(file.path(base,"results_TAIR12/analysis/bliss_vs_gdna_enrichment_minimap2.csv")) # MAPQ>=0 + gDNA
q1 <- q1[sample=="BA1"]; gd <- gd[sample=="BA1"]
# MAPQ>=1 per-read fold (sites per read, relative to genome mean)
q1[, per1k := sites/events*1000]; gper1k <- sum(q1$sites)/sum(q1$events)*1000
q1[, perread_fold := per1k/gper1k]
m <- merge(q1[,.(region, mapq1_perMb=fold_events, mapq1_perRead=perread_fold)],
           gd[,.(region, mapq0_perMb=fold_nogdna, gdna=enrichment_vs_gdna)], by="region")
long <- melt(m, id.vars="region", variable.name="method", value.name="fold")
lab <- c(arms="Arm",pericentromere="Pericentromere",centromere="Centromere",nor_45s="45S rDNA (NOR)",rdna_5s="5S rDNA")
long$region <- factor(lab[long$region], levels=rev(c("45S rDNA (NOR)","5S rDNA","Pericentromere","Centromere","Arm")))
mlab <- c(mapq0_perMb="MAPQ≥0\nper Mb of DNA\n(multimappers kept)",
          mapq1_perMb="MAPQ≥1\nper Mb of DNA\n(multimappers dropped)",
          mapq1_perRead="MAPQ≥1\nper mapped read",
          gdna="MAPQ≥0\n÷ WGA gDNA\n(unbiased)")
long$method <- factor(mlab[as.character(long$method)], levels=mlab)
long[, l2 := log2(pmax(fold,1e-3))]
p <- ggplot(long, aes(method, region, fill=l2)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.2f",fold)), size=3.6) +
  scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B", midpoint=0,
                       breaks=log2(c(0.25,0.5,1,2,4,8)), labels=c("0.25","0.5","1","2","4","8"),
                       name="fold\n(log2)") +
  labs(title="The same sBLISS data, four ways: enrichment depends on MAPQ handling × normalisation",
       subtitle="BA1, minimap2. Fold vs genome expectation. Note the 45S NOR: depleted (drop multimappers) → enriched (keep) → 12× (÷gDNA, the unbiased value).",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold"), panel.grid=element_blank(),
        axis.text.x=element_text(size=8))
ggsave(file.path(FIG,"master_comparison.png"), p, width=11, height=5.5, dpi=150)
cat("Saved master_comparison.png\n"); print(dcast(long, region~method, value.var="fold"))
