#!/usr/bin/env python3
"""
sBLISS read pre-processing — exact algorithm from Hidmi et al. 2024:

For each FASTQ entry in R1:
    UMI            = line2[:8]
    Internal_Index = line2[8:16]
    line1          = line1 + '+' + Internal_Index + '+' + UMI
    line2          = line2[16:]
    line4          = line4[16:]

Reads are accepted only if Hamming distance of Internal_Index to the
expected barcode is <= max_mismatches (default 1).

Output header format: @original_name+BARCODE+UMI
This is compatible with umi_tools dedup --umi-separator=+
(umi_tools takes everything after the last '+' as the UMI).
"""

import argparse
import gzip
import sys


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input",  help="Input R1 FASTQ (gz or plain)")
    p.add_argument("output", help="Output FASTQ.gz")
    p.add_argument("--barcode", required=True,
                   help="Expected Internal_Index (8 bp), e.g. CATCACGC for BA1")
    p.add_argument("--max-mm", type=int, default=1,
                   help="Max Hamming mismatches allowed in barcode (default 1)")
    p.add_argument("--stats", default=None, help="Optional stats output file")
    return p.parse_args()


def main():
    args = parse_args()
    target = args.barcode
    if len(target) != 8:
        sys.exit(f"ERROR: barcode must be 8 bp, got {len(target)}: {target}")

    total = accepted = rejected_mm = rejected_short = 0

    with open_fastq(args.input) as fin, \
         gzip.open(args.output, "wt") as fout:

        while True:
            line1 = fin.readline()
            if not line1:
                break
            line2 = fin.readline().rstrip("\n")
            line3 = fin.readline()
            line4 = fin.readline().rstrip("\n")

            total += 1

            if len(line2) < 16:
                rejected_short += 1
                continue

            UMI            = line2[:8]
            Internal_Index = line2[8:16]

            discrepancies = hamming(Internal_Index, target)
            if discrepancies > args.max_mm:
                rejected_mm += 1
                continue

            # Modify header: append +Internal_Index+UMI to the read NAME
            # (before the space), so it survives SAM/BAM QNAME encoding.
            # umi_tools dedup --umi-separator=+ reads UMI from after last '+'.
            raw = line1[1:].rstrip("\n")   # strip '@' and newline
            if " " in raw:
                name, desc = raw.split(" ", 1)
                new_header = "@" + name + "+" + Internal_Index + "+" + UMI + " " + desc + "\n"
            else:
                new_header = "@" + raw + "+" + Internal_Index + "+" + UMI + "\n"
            new_seq    = line2[16:] + "\n"
            new_qual   = line4[16:] + "\n"

            fout.write(new_header)
            fout.write(new_seq)
            fout.write(line3)       # '+' line unchanged
            fout.write(new_qual)

            accepted += 1

    if args.stats:
        with open(args.stats, "w") as f:
            f.write(f"total_reads\t{total}\n")
            f.write(f"accepted\t{accepted}\n")
            f.write(f"rejected_barcode_mismatch\t{rejected_mm}\n")
            f.write(f"rejected_too_short\t{rejected_short}\n")
            f.write(f"barcode_pass_rate\t{accepted/total:.4f}\n")

    print(f"Total: {total}  Accepted: {accepted}  "
          f"Rejected (barcode): {rejected_mm}  Rejected (short): {rejected_short}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
