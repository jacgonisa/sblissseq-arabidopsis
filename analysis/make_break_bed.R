#!/usr/bin/env Rscript
# Creates BREAK BED and BREAK BEDGRAPH files from a flat BED file.
# Implements the R snippet from Hidmi et al. 2024 (STAR Protocols):
#
#   bDF[bDF$STRAND=="+",]$END   <- bDF[bDF$STRAND=="+",]$START + 1
#   bDF[bDF$STRAND=="-",]$START <- bDF[bDF$STRAND=="-",]$END   - 1
#   bDF1 <- aggregate(SCORE ~ ., bDF, sum)    # strand-specific BREAK BED
#   bDF2 <- aggregate(SCORE ~ ., bDF1, sum)   # both-strand BREAK BEDGRAPH

options(scipen = 999)   # prevent scientific notation in coordinates (e.g. 8e+05 → 800000)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript make_break_bed.R <flat.bed> <output.break.bed> <output.break.bedgraph>")
}

flat_bed   <- args[1]
out_bed    <- args[2]
out_bdg    <- args[3]

cat("Reading flat BED:", flat_bed, "\n")
bDF <- read.table(flat_bed, header = FALSE, sep = "\t",
                  col.names = c("CHR", "START", "END", "ID", "SCORE", "STRAND"),
                  colClasses = c("character","integer","integer","character","numeric","character"))
cat("Reads loaded:", nrow(bDF), "\n")

# a. Resize each entry to width 1 at the 5' end (exact DSB coordinate)
#    + strand: 5' end = START → keep START, set END = START + 1
#    - strand: 5' end = END   → keep END,   set START = END - 1
bDF[bDF$STRAND == "+", "END"]   <- bDF[bDF$STRAND == "+", "START"] + 1
bDF[bDF$STRAND == "-", "START"] <- bDF[bDF$STRAND == "-", "END"]   - 1

# Set SCORE = 1 per read for counting; drop ID so aggregation groups by position
bDF$SCORE <- 1
bDF <- bDF[, c("CHR", "START", "END", "STRAND", "SCORE")]

# b. BDF1: aggregate breaks at same position + strand → strand-specific BREAK BED
bDF1 <- aggregate(SCORE ~ ., bDF, sum)
bDF1 <- bDF1[order(bDF1$CHR, bDF1$START), ]
cat("Unique strand-specific DSB sites:", nrow(bDF1), "\n")

# c. BDF2: aggregate both strands together → BREAK BEDGRAPH
bDF2 <- aggregate(SCORE ~ ., bDF1[, c("CHR", "START", "END", "SCORE")], sum)
bDF2 <- bDF2[order(bDF2$CHR, bDF2$START), ]
cat("Unique DSB sites (both strands):", nrow(bDF2), "\n")

# Write outputs
write.table(bDF1, out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(bDF2, out_bdg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Written BREAK BED:      ", out_bed, "\n")
cat("Written BREAK BEDGRAPH: ", out_bdg, "\n")
