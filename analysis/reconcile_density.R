#!/usr/bin/env Rscript
# Reconcile the three DSB-density normalisations with CONSISTENT region
# definitions (NOR 45S rDNA always its own category, never lumped into 'arms').
# Numerator throughout = distinct break SITES.
#   (a) raw spatial      = sites / region Mb
#   (b) per 1k reads     = sites / mapped reads * 1000      (depth-normalised)
#   (c) per mapped-Mb    = sites / aligned bases * 1e6      (depth-normalised)
# (b) and (c) differ only by mean read length, so they must give the same ranking.
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
REG  <- file.path(base, "results/analysis/regions")
df <- read.csv(file.path(base, "results/analysis/dsb_per_mapped_mb.csv"))  # sites,events,mapped_bases

# region sizes (bp) from the BED files
bedlen <- function(f) { b <- read.table(file.path(REG, f)); sum(b$V3 - b$V2) }
size <- c(arms_noNOR = bedlen("arms_noNOR.bed"),
          pericentromere = bedlen("pericentromere.bed"),
          centromere = bedlen("centromere.bed"),
          nor = bedlen("nor.bed"))

df <- df %>% mutate(
  region_mb     = size[region] / 1e6,
  per_region_Mb = dsb_sites / region_mb,          # (a) raw spatial
  per_1k_reads  = dsb_sites / dsb_events * 1000,   # (b) events == #reads
  per_mapped_Mb = sites_per_mapped_mb,            # (c)
  mean_read_len = mapped_bases / dsb_events)

write.csv(df[,c("sample","region","dsb_sites","dsb_events","mapped_bases",
                "per_region_Mb","per_1k_reads","per_mapped_Mb","mean_read_len")],
          file.path(base,"results/analysis/density_reconciliation.csv"), row.names = FALSE)
cat("=== reconciliation table ===\n"); print(df %>%
  select(sample, region, per_region_Mb, per_1k_reads, per_mapped_Mb, mean_read_len) %>%
  mutate(across(where(is.numeric), round, 1)) %>% as.data.frame())
cat(sprintf("\nProportionality check (c)/(b): mean=%.2f sd=%.2f  (~1000/readlen)\n",
            mean(df$per_mapped_Mb/df$per_1k_reads), sd(df$per_mapped_Mb/df$per_1k_reads)))
cat(sprintf("Pearson r between (b) and (c) across all region x sample: %.4f\n",
            cor(df$per_1k_reads, df$per_mapped_Mb)))

# ── figure: 3 metrics side by side, consistent regions ────────────────────────
lab <- c(arms_noNOR="Arm", pericentromere="Pericentromere",
         centromere="Centromere", nor="rDNA (NOR)")
long <- df %>%
  select(sample, region, per_region_Mb, per_1k_reads, per_mapped_Mb) %>%
  pivot_longer(c(per_region_Mb, per_1k_reads, per_mapped_Mb),
               names_to = "metric", values_to = "value")
long$region <- factor(lab[long$region],
  levels = c("Arm","Pericentromere","Centromere","rDNA (NOR)"))
long$sample <- factor(long$sample, levels = c("BA1","BA2","old_BA1_BA2"))
long$metric <- factor(long$metric,
  levels = c("per_region_Mb","per_1k_reads","per_mapped_Mb"),
  labels = c("(a) Raw spatial\n(sites per Mb of DNA)",
             "(b) Per 1,000 mapped reads\n(depth-normalised)",
             "(c) Per mapped-Mb\n(depth-normalised)"))

cols <- c(Arm="#1f77b4", Pericentromere="#ff7f0e", Centromere="#9467bd", "rDNA (NOR)"="#d62728")
p <- ggplot(long, aes(region, value, fill = region)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.3) +
  facet_grid(metric ~ sample, scales = "free_y", switch = "y") +
  scale_fill_manual(values = cols, name = NULL) +
  labs(title = "Reconciling DSB-density normalisations (consistent regions; NOR separated)",
       subtitle = "Numerator = distinct break sites. (b) and (c) are proportional (r=0.998) and agree:\nArm ~ Centromere > Pericentromere; rDNA the outlier. (a) differs only because of mappability.",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.placement = "outside")
ggsave(file.path(base,"results/analysis/figures/density_reconciliation.png"), p,
       width = 11, height = 8, dpi = 150)
ggsave(file.path(base,"results/analysis/figures/density_reconciliation.pdf"), p,
       width = 11, height = 8)
cat("\nSaved density_reconciliation.png\n")
