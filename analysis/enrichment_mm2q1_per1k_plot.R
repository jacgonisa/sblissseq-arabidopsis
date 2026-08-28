#!/usr/bin/env Rscript
# Depth-normalised compartment density: DSB sites per 1,000 mapped reads (minimap2 MAPQ>=1).
# Contrast with the per-Mb-of-DNA fold (mappability-confounded).
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"; FIG <- file.path(base,"results_TAIR12/analysis/figures")
d <- fread(file.path(base,"results_TAIR12/analysis/compartment_enrichment_mm2q1.csv"))
d[, per1k := sites/events*1000]
lab <- c(arms="Arm", pericentromere="Pericentromere", centromere="Centromere",
         nor_45s="45S rDNA (NOR)", rdna_5s="5S rDNA")
d$region <- factor(lab[d$region], levels=c("Centromere","Arm","Pericentromere","5S rDNA","45S rDNA (NOR)"))
d$sample <- factor(d$sample, levels=c("BA1","BA2","old_BA1_BA2"))
cols <- c(BA1="#1f77b4",BA2="#2ca02c",old_BA1_BA2="#d62728")

pL <- ggplot(d, aes(region, fold_events, fill=sample)) +
  geom_col(position=position_dodge(0.8), width=0.75) + geom_hline(yintercept=1, linetype=2, colour="grey40") +
  geom_text(aes(label=sprintf("%.2f",fold_events)), position=position_dodge(0.8), vjust=-0.3, size=2.5) +
  scale_fill_manual(values=cols, name=NULL) +
  labs(title="A  Per Mb of DNA (spatial) — mappability-confounded",
       subtitle="Fold vs genome mean. Pericentromere>Arms; cen/NOR low = few reads map there.",
       x=NULL, y="fold (DSB/Mb vs genome mean)") +
  theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold"),
       axis.text.x=element_text(angle=20,hjust=1),panel.grid.major.x=element_blank())

pR <- ggplot(d, aes(region, per1k, fill=sample)) +
  geom_col(position=position_dodge(0.8), width=0.75) +
  geom_text(aes(label=sprintf("%.0f",per1k)), position=position_dodge(0.8), vjust=-0.3, size=2.5) +
  scale_fill_manual(values=cols, name=NULL) +
  labs(title="B  Per 1,000 mapped reads (depth-normalised)",
       subtitle="DSBs per mapped read. ~Uniform (Centromere~Arms~Peri); only 45S NOR lower.",
       x=NULL, y="distinct DSB sites / 1,000 mapped reads") +
  theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold"),
       axis.text.x=element_text(angle=20,hjust=1),panel.grid.major.x=element_blank())

ggsave(file.path(FIG,"compartment_enrichment_mm2q1_norm.png"), pL|pR, width=13, height=6, dpi=150)
cat("Saved compartment_enrichment_mm2q1_norm.png\n")
