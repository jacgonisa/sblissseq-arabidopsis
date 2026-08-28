#!/usr/bin/env Rscript
# BiCroLab Plot 3 — DSB composition: promoter/gene/intergenic + gene biotype.
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
FIG  <- file.path(base,"results_TAIR12/analysis/figures")
comp <- fread(file.path(base,"results_TAIR12/analysis/bicrolab_composition.csv"))
bt   <- fread(file.path(base,"results_TAIR12/analysis/bicrolab_biotype.csv"))
lv   <- c("BA1","BA2","old_BA1_BA2")

# A: promoter/gene/intergenic %
comp[, pct:=100*dsb_umis/sum(dsb_umis), by=sample]
comp$sample <- factor(comp$sample, levels=lv)
comp$category <- factor(comp$category, levels=c("intergenic","gene","promoter"))
pA <- ggplot(comp, aes(sample, pct, fill=category)) +
  geom_col(width=0.7, colour="white") +
  geom_text(aes(label=sprintf("%.0f%%",pct)), position=position_stack(vjust=0.5), size=3, colour="white") +
  scale_fill_manual(values=c(promoter="#1f77b4", gene="#2ca02c", intergenic="#cccccc"), name=NULL) +
  labs(title="A  DSB composition (promoter / gene / intergenic)",
       subtitle="% of total DSB UMIs. Promoter = TSS −2 kb.", x=NULL, y="% of DSB UMIs") +
  theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold"),
       panel.grid.major.x=element_blank())

# B: gene biotype (collapse small to 'other')
bt[, biotype:=sub("^=","",biotype)]
bt[, pct:=100*dsb_umis/sum(dsb_umis), by=sample]
keep <- bt[, .(m=max(pct)), by=biotype][m>=2, biotype]
bt[, grp:=ifelse(biotype %in% keep, biotype, "other")]
btg <- bt[, .(pct=sum(pct)), by=.(sample,grp)]
btg$sample <- factor(btg$sample, levels=lv)
pB <- ggplot(btg, aes(sample, pct, fill=grp)) +
  geom_col(width=0.7, colour="white") +
  geom_text(aes(label=ifelse(pct>=4,sprintf("%.0f%%",pct),"")), position=position_stack(vjust=0.5), size=3, colour="white") +
  scale_fill_brewer(palette="Set2", name="gene biotype") +
  labs(title="B  Genic DSBs by gene biotype",
       subtitle="% of genic DSB UMIs (rRNA = rDNA genes)", x=NULL, y="% of genic DSB UMIs") +
  theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold"),
       panel.grid.major.x=element_blank())

ggsave(file.path(FIG,"bicrolab_composition.png"), pA|pB, width=11, height=5.5, dpi=150)
cat("Saved bicrolab_composition.png\n")
print(dcast(comp, sample~category, value.var="pct"))
