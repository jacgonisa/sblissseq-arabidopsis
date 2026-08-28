#!/usr/bin/env Rscript
# Mappability-aware compartment DSB density (TAIR12). Three views:
#  A. breaks per Mb of DNA (all reads)         -> rDNA dominates (break production)
#  B. uniquely-mappable % (GenMap k=127)        -> why unique metrics fail in repeats
#  C. q10 sites per mappable-Mb                 -> break density of OBSERVABLE sequence
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
FIG  <- file.path(base, "results_TAIR12/analysis/figures")
d <- read.csv(file.path(base, "results_TAIR12/analysis/mappability_compartments.csv"))

lab <- c(arms="Arm", pericentromere="Pericentromere", centromere="Centromere",
         nor_45s="45S rDNA (NOR)", rdna_5s="5S rDNA")
ord <- c("Arm","Pericentromere","Centromere","5S rDNA","45S rDNA (NOR)")
d$region <- factor(lab[d$region], levels=ord)
d$sample <- factor(d$sample, levels=c("BA1","BA2","old_BA1_BA2"))
cols <- c(BA1="#1f77b4", BA2="#2ca02c", old_BA1_BA2="#d62728")

pA <- ggplot(d, aes(region, all_sites_per_Mb_DNA, fill=sample)) +
  geom_col(position=position_dodge(0.8), width=0.75) +
  scale_fill_manual(values=cols, name=NULL) +
  labs(title="A  Break density per Mb of DNA (all reads) — rDNA is the most break-dense locus",
       subtitle="Distinct break sites per Mb of region. The 45S NOR towers ~8x over the arms.",
       x=NULL, y="Break sites / Mb of DNA") +
  theme_bw(base_size=11)+theme(legend.position="top",plot.title=element_text(face="bold"),
       axis.text.x=element_text(angle=20,hjust=1),panel.grid.major.x=element_blank())

mp <- d %>% filter(sample=="BA1") %>% select(region, mappable_pct)
pB <- ggplot(mp, aes(region, mappable_pct, fill=region)) +
  geom_col(width=0.7) + geom_text(aes(label=sprintf("%.1f%%",mappable_pct)), vjust=-0.3, size=3) +
  scale_fill_manual(values=c(Arm="#1f77b4",Pericentromere="#ff7f0e",Centromere="#9467bd",
       "5S rDNA"="#17becf","45S rDNA (NOR)"="#d62728"), guide="none") +
  labs(title="B  Uniquely-mappable fraction (GenMap, k=127 = read length)",
       subtitle="Centromere 3.9% and 45S NOR 2.0% are near-invisible to unique mapping — this is why per-read metrics hide them.",
       x=NULL, y="% uniquely mappable") + ylim(0,105) +
  theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),
       axis.text.x=element_text(angle=20,hjust=1),panel.grid.major.x=element_blank())

pC <- ggplot(d, aes(region, q10_sites_per_mappable_Mb, fill=sample)) +
  geom_col(position=position_dodge(0.8), width=0.75) +
  scale_fill_manual(values=cols, name=NULL) +
  labs(title="C  Break density of the OBSERVABLE sequence (MAPQ>=10 sites per uniquely-mappable Mb)",
       subtitle="Within mappable DNA: 5S rDNA & centromere exceed the arms. The 45S's visible 2% looks arm-like but cannot represent the array.",
       x=NULL, y="MAPQ>=10 sites / mappable-Mb") +
  theme_bw(base_size=11)+theme(legend.position="none",plot.title=element_text(face="bold"),
       axis.text.x=element_text(angle=20,hjust=1),panel.grid.major.x=element_blank())

ggsave(file.path(FIG,"mappability_compartments.png"), pA/pB/pC, width=10, height=12, dpi=150)
cat("Saved mappability_compartments.png\n")
