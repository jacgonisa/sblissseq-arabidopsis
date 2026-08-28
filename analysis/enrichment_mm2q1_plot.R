#!/usr/bin/env Rscript
# Compartment DSB fold-enrichment (minimap2 MAPQ>=1, uniquely-mapped), with mappability context.
suppressMessages({ library(data.table); library(ggplot2) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
FIG  <- file.path(base,"results_TAIR12/analysis/figures")
d <- fread(file.path(base,"results_TAIR12/analysis/compartment_enrichment_mm2q1.csv"))
lab <- c(arms="Arm", pericentromere="Pericentromere", centromere="Centromere",
         nor_45s="45S rDNA (NOR)", rdna_5s="5S rDNA")
mappable <- c("Arm"=97.5,"Pericentromere"=80.5,"Centromere"=3.9,"5S rDNA"=5.3,"45S rDNA (NOR)"=2.0)
d$region <- factor(lab[d$region], levels=c("Pericentromere","Arm","5S rDNA","Centromere","45S rDNA (NOR)"))
d$sample <- factor(d$sample, levels=c("BA1","BA2","old_BA1_BA2"))
cols <- c(BA1="#1f77b4",BA2="#2ca02c",old_BA1_BA2="#d62728")

p <- ggplot(d, aes(region, fold_events, fill=sample)) +
  geom_col(position=position_dodge(0.8), width=0.75) +
  geom_hline(yintercept=1, linetype=2, colour="grey40") +
  geom_text(aes(label=sprintf("%.2f",fold_events)), position=position_dodge(0.8), vjust=-0.3, size=2.7) +
  annotate("text", x=1:5, y=-0.07,
           label=paste0(mappable[levels(d$region)],"%"), size=2.6, colour="grey35") +
  scale_fill_manual(values=cols, name=NULL) +
  labs(title="Compartment DSB enrichment — minimap2, MAPQ>=1 (uniquely-mapped)",
       subtitle="Fold over genome-mean DSB density per Mb of DNA. Grey % = uniquely-mappable fraction (k=127).\nPericentromere > Arms; Centromere & 45S NOR appear depleted = mappability blind spot, NOT true absence.",
       x=NULL, y="Fold enrichment (DSB/Mb vs genome mean)") +
  theme_bw(base_size=11) + theme(legend.position="bottom", plot.title=element_text(face="bold"),
       panel.grid.major.x=element_blank())
ggsave(file.path(FIG,"compartment_enrichment_mm2q1.png"), p, width=9, height=6, dpi=150)
cat("Saved compartment_enrichment_mm2q1.png\n")
