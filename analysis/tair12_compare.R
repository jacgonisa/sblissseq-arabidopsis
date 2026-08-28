#!/usr/bin/env Rscript
# TAIR12 compartment density + the ragtag(collapsed) vs TAIR12(resolved) rDNA comparison.
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })
base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
FIG  <- file.path(base, "results_TAIR12/analysis/figures")

t12 <- read.csv(file.path(base, "results_TAIR12/analysis/compartment_density.csv"))
rag <- read.csv(file.path(base, "results/analysis/density_reconciliation.csv"))

# ── Panel A: TAIR12 compartment density (sites per 1k reads), all compartments ─
labA <- c(arms="Arm", pericentromere="Pericentromere", centromere="Centromere",
          nor_45s="45S rDNA (NOR)", rdna_5s="5S rDNA")
a <- t12 %>% mutate(region=factor(labA[region],
        levels=c("Arm","Centromere","Pericentromere","5S rDNA","45S rDNA (NOR)")),
        sample=factor(sample, levels=c("BA1","BA2","old_BA1_BA2")))
colA <- c(Arm="#1f77b4",Centromere="#9467bd",Pericentromere="#ff7f0e",
          "5S rDNA"="#17becf","45S rDNA (NOR)"="#d62728")
pA <- ggplot(a, aes(sample, sites_per_1k_reads, fill=region)) +
  geom_col(position=position_dodge(0.8), width=0.75, colour="white", linewidth=0.3) +
  geom_text(aes(label=round(sites_per_1k_reads)), position=position_dodge(0.8),
            vjust=-0.3, size=2.5) +
  scale_fill_manual(values=colA, name=NULL) +
  labs(title="A  TAIR12 compartment DSB density (per-read)",
       subtitle="Distinct sites per 1,000 mapped reads. Arm ~ Centromere > Pericentromere; 5S ~ arms; 45S still partly saturated.",
       x=NULL, y="Sites / 1,000 reads") +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(face="bold"), panel.grid.major.x=element_blank())

# ── Panel B: rDNA paradox resolution — ragtag vs TAIR12 (45S NOR) ──────────────
rag_nor <- rag %>% filter(region=="nor") %>%
  transmute(sample, assembly="ragtag (collapsed, 2 Mb)",
            sites_per_1k=per_1k_reads, reads_per_site=dsb_events/dsb_sites)
t12_nor <- t12 %>% filter(region=="nor_45s") %>%
  transmute(sample, assembly="TAIR12 (resolved, 9.3 Mb)",
            sites_per_1k=sites_per_1k_reads, reads_per_site=reads_per_site)
nor <- bind_rows(rag_nor, t12_nor)
nor$sample <- factor(nor$sample, levels=c("BA1","BA2","old_BA1_BA2"))

pB <- ggplot(nor, aes(sample, sites_per_1k, fill=assembly)) +
  geom_col(position=position_dodge(0.8), width=0.7, colour="white", linewidth=0.3) +
  geom_text(aes(label=round(sites_per_1k)), position=position_dodge(0.8), vjust=-0.3, size=3) +
  scale_fill_manual(values=c("ragtag (collapsed, 2 Mb)"="#bbbbbb",
                             "TAIR12 (resolved, 9.3 Mb)"="#d62728"), name=NULL) +
  labs(title="B  45S rDNA per-read density recovers on the resolved assembly",
       subtitle="Distinct break sites per 1,000 reads in the NOR. Collapse halved the apparent density.",
       x=NULL, y="Sites / 1,000 reads (NOR)") +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(face="bold"), panel.grid.major.x=element_blank())

pC <- ggplot(nor, aes(sample, reads_per_site, fill=assembly)) +
  geom_col(position=position_dodge(0.8), width=0.7, colour="white", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.2f",reads_per_site)), position=position_dodge(0.8), vjust=-0.3, size=3) +
  geom_hline(yintercept=1.05, linetype=2, colour="grey40") +
  annotate("text", x=0.7, y=1.15, label="arms ~1.05", size=2.8, colour="grey30", hjust=0) +
  scale_fill_manual(values=c("ragtag (collapsed, 2 Mb)"="#bbbbbb",
                             "TAIR12 (resolved, 9.3 Mb)"="#d62728"), name=NULL) +
  labs(title="C  Read redundancy in the NOR (saturation proxy)",
       subtitle="Reads per unique site. Resolution roughly halves redundancy; residual = multimapping + 8 bp UMI ceiling.",
       x=NULL, y="Reads / unique site") +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(face="bold"), panel.grid.major.x=element_blank())

ggsave(file.path(FIG,"tair12_compartments.png"), pA, width=9, height=5.5, dpi=150)
ggsave(file.path(FIG,"rdna_resolution_compare.png"), pB | pC, width=12, height=5.5, dpi=150)
cat("Saved tair12_compartments.png and rdna_resolution_compare.png\n")
