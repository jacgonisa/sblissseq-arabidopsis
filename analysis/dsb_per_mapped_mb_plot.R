#!/usr/bin/env Rscript
# Barplot: distinct DSB sites per Mb of uniquely-mapped read sequence, by compartment.
# rDNA (NOR) is its own category. Input: results/analysis/dsb_per_mapped_mb.csv
suppressMessages({ library(ggplot2); library(dplyr) })

base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
df <- read.csv(file.path(base, "results/analysis/dsb_per_mapped_mb.csv"))

df$region <- factor(df$region,
  levels = c("nor","centromere","pericentromere","arms_noNOR"),
  labels = c("rDNA (NOR)","Centromere","Pericentromere","Arm (NOR-excl.)"))
df$sample <- factor(df$sample, levels = c("BA1","BA2","old_BA1_BA2"))

cols <- c("rDNA (NOR)"="#d62728","Centromere"="#9467bd",
          "Pericentromere"="#ff7f0e","Arm (NOR-excl.)"="#1f77b4")

p <- ggplot(df, aes(sample, sites_per_mapped_mb, fill = region)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = format(round(sites_per_mapped_mb), big.mark=",")),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = cols) +
  labs(title = "DSB density per Mb of uniquely-mapped reads",
       subtitle = "Distinct break sites per Mb of aligned bases (mapped read-Mb). NOR rDNA shown separately.",
       x = NULL, y = "Distinct DSB sites per mapped-Mb", fill = "Compartment") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(), legend.position = "bottom")

ggsave(file.path(base, "results/analysis/figures/dsb_per_mapped_mb.png"), p,
       width = 9, height = 5.5, dpi = 150)
ggsave(file.path(base, "results/analysis/figures/dsb_per_mapped_mb.pdf"), p,
       width = 9, height = 5.5)
cat("Saved dsb_per_mapped_mb.png\n")
