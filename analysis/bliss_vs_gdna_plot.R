#!/usr/bin/env Rscript
# BLISS DSB enrichment: without gDNA (uniform expectation) vs ÷gDNA (DNA-content corrected).
# arg1 = tag (minimap2 / bowtie2)
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
args <- commandArgs(trailingOnly=TRUE); TAG <- ifelse(length(args)>=1,args[1],"minimap2")
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"; FIG <- file.path(base,"results_TAIR12/analysis/figures")
d <- fread(file.path(base, sprintf("results_TAIR12/analysis/bliss_vs_gdna_enrichment_%s.csv",TAG)))
lab <- c(arms="Arm",pericentromere="Pericentromere",centromere="Centromere",nor_45s="45S rDNA (NOR)",rdna_5s="5S rDNA")
ord <- c("45S rDNA (NOR)","5S rDNA","Pericentromere","Centromere","Arm")
d$region <- factor(lab[d$region], levels=ord); d$sample <- factor(d$sample, levels=c("BA1","BA2","old_BA1_BA2"))
cols <- c(BA1="#1f77b4",BA2="#2ca02c",old_BA1_BA2="#d62728")
mk <- function(col, ttl, sub) ggplot(d, aes(region, get(col), fill=sample)) +
  geom_col(position=position_dodge(0.8), width=0.75) + geom_hline(yintercept=1, linetype=2, colour="grey40") +
  geom_text(aes(label=sprintf("%.1f",get(col))), position=position_dodge(0.8), vjust=-0.3, size=2.5) +
  scale_fill_manual(values=cols,name=NULL) + scale_y_log10() +
  labs(title=ttl, subtitle=sub, x=NULL, y=paste0(col," (log scale)")) +
  theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold"),
       axis.text.x=element_text(angle=20,hjust=1),panel.grid.major.x=element_blank())
pA <- mk("fold_nogdna","A  Without gDNA (over uniform expectation)","Multimappers kept; NOR dominated by copy-number+breaks.")
pB <- mk("enrichment_vs_gdna","B  Divided by WGA gDNA (DNA-content corrected)","Real DSB enrichment: 45S rDNA genuinely break-rich (~12-18x).")
ggsave(file.path(FIG, sprintf("bliss_vs_gdna_%s.png",TAG)), pA|pB, width=13, height=6, dpi=150)
cat("Saved bliss_vs_gdna_",TAG,".png\n")
