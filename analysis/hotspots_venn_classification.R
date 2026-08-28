#!/usr/bin/env Rscript
# 3-way Venn of DSB hotspots across samples + region classification.
# Categories: arm / pericentromere / centromere, with NOR 45S rDNA as its own class.
# Hotspot BEDs carry: chr start end name score strand fold padj region is_NOR
# Usage: Rscript hotspots_venn_classification.R

options(scipen = 999)
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

base <- "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
HS   <- Sys.getenv("BLISS_HSDIR",  file.path(base, "results/analysis/hotspots"))
FIG  <- Sys.getenv("BLISS_FIGDIR", file.path(base, "results/analysis/figures"))
samples <- c("BA1", "BA2", "old_BA1_BA2")

read_hs <- function(s) {
  d <- read.table(file.path(HS, sprintf("%s_1kb.hotspots.bed", s)), sep = "\t",
                  stringsAsFactors = FALSE)
  names(d)[c(1,2,3,9,10)] <- c("chr","start","end","region","is_NOR")
  d$key <- paste(d$chr, d$start, sep = ":")
  d$category <- ifelse(d$is_NOR == "TRUE" | d$is_NOR == TRUE, "rDNA (NOR)", d$region)
  d
}
hs <- lapply(samples, read_hs); names(hs) <- samples
sets     <- lapply(hs, function(d) d$key)                       # all hotspots
sets_gen <- lapply(hs, function(d) d$key[d$category != "rDNA (NOR)"])  # non-rDNA

# ── 3-circle Venn (hand-drawn, dependency-free) ───────────────────────────────
circle <- function(cx, cy, r, n = 200) {
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = cx + r*cos(t), y = cy + r*sin(t))
}
venn_counts <- function(S) {
  A <- S[[1]]; B <- S[[2]]; C <- S[[3]]
  c(A     = length(setdiff(A, union(B,C))),
    B     = length(setdiff(B, union(A,C))),
    C     = length(setdiff(C, union(A,B))),
    AB    = length(setdiff(intersect(A,B), C)),
    AC    = length(setdiff(intersect(A,C), B)),
    BC    = length(setdiff(intersect(B,C), A)),
    ABC   = length(Reduce(intersect, list(A,B,C))))
}
draw_venn <- function(S, title) {
  cn <- venn_counts(S)
  cx <- c(-0.45, 0.45, 0);  cy <- c(0.35, 0.35, -0.45); r <- 0.85
  circ <- bind_rows(
    cbind(circle(cx[1],cy[1],r), set="BA1"),
    cbind(circle(cx[2],cy[2],r), set="BA2"),
    cbind(circle(cx[3],cy[3],r), set="old_BA1_BA2"))
  lab <- data.frame(
    x = c(-0.95, 0.95, 0,    0,    -0.62, 0.62, 0),
    y = c(0.75, 0.75, -1.05, 0.78, -0.30, -0.30, 0.08),
    n = c(cn["A"],cn["B"],cn["C"],cn["AB"],cn["AC"],cn["BC"],cn["ABC"]))
  setlab <- data.frame(x=c(-1.05,1.05,0), y=c(1.35,1.35,-1.45),
                       set=c("BA1","BA2","old_BA1_BA2"))
  ggplot() +
    geom_polygon(data=circ, aes(x,y,fill=set,group=set), alpha=0.30, colour="grey30") +
    geom_text(data=lab, aes(x,y,label=format(n, big.mark=",")), size=4.2, fontface="bold") +
    geom_text(data=setlab, aes(x,y,label=set,colour=set), size=4.5, fontface="bold", show.legend=FALSE) +
    scale_fill_manual(values=c(BA1="#1f77b4",BA2="#2ca02c",old_BA1_BA2="#d62728")) +
    scale_colour_manual(values=c(BA1="#1f77b4",BA2="#2ca02c",old_BA1_BA2="#d62728")) +
    coord_equal(xlim=c(-1.7,1.7), ylim=c(-1.7,1.7)) +
    labs(title=title) +
    theme_void(base_size=12) +
    theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5))
}

pv1 <- draw_venn(sets,     sprintf("All hotspots (n: BA1=%d, BA2=%d, old=%d)",
                                   length(sets[[1]]),length(sets[[2]]),length(sets[[3]])))
pv2 <- draw_venn(sets_gen, sprintf("Genuine only, rDNA excluded (n: %d / %d / %d)",
                                   length(sets_gen[[1]]),length(sets_gen[[2]]),length(sets_gen[[3]])))

suppressMessages(library(patchwork))
ggsave(file.path(FIG,"hotspots_venn.png"), pv1 | pv2, width=12, height=6, dpi=150)
ggsave(file.path(FIG,"hotspots_venn.pdf"), pv1 | pv2, width=12, height=6)
cat("Saved hotspots_venn.png\n")
cat("\nVenn counts (all):\n");     print(venn_counts(sets))
cat("\nVenn counts (genuine):\n"); print(venn_counts(sets_gen))

# ── Region classification barplot (per sample + robust set) ───────────────────
cls <- bind_rows(lapply(samples, function(s)
  data.frame(sample=s, category=hs[[s]]$category)))
# robust set
rb <- read.table(file.path(HS,"robust_hotspots.annotated.tsv"), header=TRUE, sep="\t",
                 stringsAsFactors=FALSE)
rb$category <- ifelse(rb$is_NOR=="TRUE"|rb$is_NOR==TRUE, "rDNA (NOR)", rb$region)
cls <- bind_rows(cls, data.frame(sample="Robust (3-way)", category=rb$category))

cls$category <- factor(cls$category,
  levels=c("rDNA (NOR)","centromere","pericentromere","arm"),
  labels=c("rDNA (NOR)","Centromere","Pericentromere","Arm"))
cls$sample <- factor(cls$sample, levels=c("BA1","BA2","old_BA1_BA2","Robust (3-way)"))

cnt <- cls %>% count(sample, category)
cat("\nClassification counts:\n"); print(as.data.frame(cnt))

pcls <- ggplot(cnt, aes(sample, n, fill=category)) +
  geom_col(width=0.7, colour="white", linewidth=0.3) +
  geom_text(aes(label=n), position=position_stack(vjust=0.5), size=3, colour="white") +
  scale_fill_manual(values=c("rDNA (NOR)"="#d62728","Centromere"="#9467bd",
                             "Pericentromere"="#ff7f0e","Arm"="#1f77b4")) +
  labs(title="DSB hotspot classification by chromosomal compartment",
       subtitle="NOR 45S rDNA shown as a separate (artefact) category",
       x=NULL, y="Number of 1 kb hotspots", fill="Category") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold"), panel.grid.major.x=element_blank())

ggsave(file.path(FIG,"hotspots_classification.png"), pcls, width=8, height=6, dpi=150)
ggsave(file.path(FIG,"hotspots_classification.pdf"), pcls, width=8, height=6)
cat("\nSaved hotspots_classification.png\n")
