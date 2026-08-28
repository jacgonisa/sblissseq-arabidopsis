#!/usr/bin/env python3
"""
Extract 5' end coordinates of deduplicated sBLISS reads as DSB sites.

R1 reads from the genomic break-site side; after deduplication, the 5' end
of each alignment = the exact DSB position at single-nucleotide resolution.
"""
import argparse
import sys
import pysam

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("bam")
    p.add_argument("--stats", default=None)
    return p.parse_args()

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")

    counts = {}
    total = 0
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        total += 1
        chrom = read.reference_name
        if read.is_reverse:
            pos = read.reference_end  # 5' end on minus strand
            strand = "-"
        else:
            pos = read.reference_start  # 5' end on plus strand (0-based)
            strand = "+"
        key = (chrom, pos, strand)
        counts[key] = counts.get(key, 0) + 1

    bam.close()

    unique_sites = len(counts)
    for (chrom, pos, strand), n in sorted(counts.items()):
        # BED: chrom, start (0-based), end (1-based), name, score, strand
        print(f"{chrom}\t{pos}\t{pos+1}\tDSB\t{n}\t{strand}")

    if args.stats:
        with open(args.stats, "w") as f:
            f.write(f"total_reads\t{total}\n")
            f.write(f"unique_dsb_sites\t{unique_sites}\n")
            f.write(f"mean_reads_per_site\t{total/unique_sites:.2f}\n" if unique_sites else "mean_reads_per_site\tNA\n")

if __name__ == "__main__":
    main()
