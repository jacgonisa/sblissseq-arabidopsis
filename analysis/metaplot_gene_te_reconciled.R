suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
# Reconciles report sections 7a (log2 BLISS/gDNA) and 7b (signal/mean overlay).
# Both use the SAME two deeptools profiles (BA1 sBLISS break density + WGA gDNA
# coverage). 7a looked noisy only because it divided sparse-5'-break CPM by dense
# coverage CPM with a +1 pseudocount. The honest, clean version is to (1) scale
# each track to its own mean (shape), overlay them [top row = 7b], then (2) take
# log2(shapeBLISS / shapeGDNA) [bottom row = corrected 7a]. Row 2 is exactly the
# vertical gap between the row-1 curves, in log space: DSB enrichment ABOVE the
# naked-DNA (gDNA) accessibility baseline, per position.
M  <- Sys.getenv("BLISS_METADIR", "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis/metaplots")
FIG<- Sys.getenv("BLISS_FIGDIR",  "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis/figures")

read_prof <- function(tab){
  tb  <- readLines(tab)
  rows<- tb[grepl("^BA1|^wga", tb)]
  rbindlist(lapply(rows, function(l){
    v <- strsplit(l, "\t")[[1]]
    data.table(track = ifelse(grepl("BA1", v[1]), "sBLISS DSB", "WGA gDNA"),
               i = seq_along(v[3:length(v)]), val = as.numeric(v[3:length(v)]))
  }))[, shape := val/mean(val), by = track]
}

xaxis <- function(d, kind){
  n <- max(d$i)
  if (kind == "tss"){
    d[, x := (i-1)/(n-1)*4000 - 2000]; list(d=d, vl=0, brk=c(-2000,0,2000),
      lbl=c("-2kb","TSS","+2kb"), xlab="bp from TSS")
  } else {
    up <- round(n*2000/7000); bo <- round(n*5000/7000)
    d[, x := i]; list(d=d, vl=c(up,bo), brk=c(up,bo), lbl=c("5' start","3' end"),
      xlab=if(kind=="te") "transposon (scaled + 2kb flanks)" else "gene body (scaled + 2kb flanks)")
  }
}

overlay <- function(tab, title, kind){
  a <- xaxis(read_prof(tab), kind); d <- a$d
  pk  <- d[track=="sBLISS DSB"][which.max(shape)]
  pkg <- d[track=="WGA gDNA"][i==pk$i]
  ggplot(d, aes(x, shape, colour=track)) +
    geom_hline(yintercept=1, linetype=3, colour="grey70") +
    geom_vline(xintercept=a$vl, linetype=2, colour="grey60") +
    geom_line(linewidth=1.05) +
    scale_colour_manual(values=c("sBLISS DSB"="#d62728","WGA gDNA"="#555555"), name=NULL) +
    scale_x_continuous(breaks=a$brk, labels=a$lbl) +
    labs(title=sprintf("%s  (sBLISS %.2f× vs gDNA %.2f×)", title, pk$shape, pkg$shape),
         x=NULL, y="signal / mean") +
    theme_bw(base_size=11) + theme(legend.position="bottom",
      plot.title=element_text(face="bold", size=9.5))
}

ratio <- function(tab, kind){
  a <- xaxis(read_prof(tab), kind); d <- a$d
  w <- dcast(d, i + x ~ track, value.var="shape")
  setnames(w, c("sBLISS DSB","WGA gDNA"), c("bliss","gdna"))
  w[, l2 := log2(bliss/gdna)]
  pk <- w[which.max(l2)]
  ggplot(w, aes(x, l2)) +
    geom_hline(yintercept=0, colour="grey40") +
    geom_vline(xintercept=a$vl, linetype=2, colour="grey60") +
    geom_area(aes(fill = l2 > 0), alpha=.35) +
    geom_line(linewidth=1.05, colour="#1a1a1a") +
    scale_fill_manual(values=c(`TRUE`="#d62728", `FALSE`="#1f6feb"), guide="none") +
    scale_x_continuous(breaks=a$brk, labels=a$lbl) +
    annotate("text", x=pk$x, y=pk$l2, label=sprintf("  peak %+.2f (=%.2f×)", pk$l2, 2^pk$l2),
             hjust=0, vjust=-0.4, size=3, fontface="bold") +
    labs(x=a$xlab, y="log2( sBLISS / gDNA )") +
    theme_bw(base_size=11)
}

specs <- list(c("TSS_prof.tab","TSS","tss"),
              c("genebody_prof.tab","Gene body","gene"),
              c("TE_prof.tab","Transposons","te"))
top <- lapply(specs, function(s) overlay(file.path(M,s[1]), s[2], s[3]))
bot <- lapply(specs, function(s) ratio  (file.path(M,s[1]), s[3]))

fig <- (top[[1]] | top[[2]] | top[[3]]) / (bot[[1]] | bot[[2]] | bot[[3]]) +
  plot_annotation(
    title = "DSB density at genes & transposons, gDNA-controlled (sections 7a+7b reconciled)",
    subtitle = "Top: each track scaled to its own mean — sBLISS departs from the flat gDNA baseline.  Bottom: log2(sBLISS/gDNA) on those shapes = DSB enrichment above naked-DNA expectation.",
    theme = theme(plot.title=element_text(face="bold", size=12),
                  plot.subtitle=element_text(size=9.5, colour="grey30")))

ggsave(file.path(FIG,"metaplot_gene_te_reconciled.png"), fig, width=13.5, height=8.2, dpi=140)

# print the numeric summary used in the report caption
for (s in specs){
  a <- xaxis(read_prof(file.path(M,s[1])), s[3])
  w <- dcast(a$d, i + x ~ track, value.var="shape")
  setnames(w, c("sBLISS DSB","WGA gDNA"), c("bliss","gdna")); w[, l2:=log2(bliss/gdna)]
  pk <- w[which.max(l2)]; tr <- w[which.min(l2)]
  cat(sprintf("%-12s  peak log2=%+.2f (=%.2f x) at x=%.0f ; min log2=%+.2f (=%.2f x) at x=%.0f\n",
              s[2], pk$l2, 2^pk$l2, pk$x, tr$l2, 2^tr$l2, tr$x))
}
cat("Saved metaplot_gene_te_reconciled.png\n")
