#!/usr/bin/env python3
"""
How are DSB HOTSPOTS distributed across features? (count/occupancy view)

Complements te_vs_gene_enrichment.py (which is per-feature *density*). Uses the
genuine (non-rDNA) robust 3-way hotspots and asks, for each feature class:
  - occupancy fold = (% of hotspots overlapping the feature) / (% of genome it covers)
  - tail metric    = % of that feature class that hosts a hotspot
Key point: median density (TE ~ gene) and hotspot occupancy (TE-enriched) can
disagree because hotspots are the upper TAIL of the density distribution.

Inputs : hotspots/robust_hotspots.annotated.tsv (col6 is_NOR), regions/*.bed6,
         reference/Col-CC.chrom.sizes.  Needs bedtools on PATH (conda activate bliss).
Outputs: results_TAIR12/analysis/hotspot_distribution.csv
         results_TAIR12/analysis/figures/hotspot_distribution.png
"""
import os, subprocess, tempfile
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.environ.get("BLISS_ROOT", "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026")
ANA  = f"{ROOT}/results_TAIR12/analysis"
REG  = os.environ.get("BLISS_REGDIR", f"{ANA}/regions")
FIG  = os.environ.get("BLISS_FIGDIR", f"{ANA}/figures")
HS   = f"{ANA}/hotspots/robust_hotspots.annotated.tsv"
SIZES= f"{ROOT}/reference/Col-CC.chrom.sizes"
FEATURES = [("genes","genes.bed6","#8c564b"), ("TEs","TEs.bed6","#9467bd"),
            ("lncRNA","lncRNA.bed","#2ca02c"), ("satellites","trash_satellites.bed","#7f7f7f")]

def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True, check=True).stdout
def nlines(p): return sum(1 for _ in open(p))

def main():
    os.makedirs(FIG, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        # genuine (non-NOR) robust hotspots + compartment tally
        hs = f"{tmp}/hs.bed"
        sh(f"awk -F'\\t' 'NR>1 && $6==\"FALSE\"{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}}' {HS} "
           f"| sort -k1,1 -k2,2n > {hs}")
        NH = nlines(hs)
        comp = pd.read_csv(hs, sep="\t", header=None, names=["c","s","e","region"])["region"].value_counts().to_dict()
        G = sum(int(l.split()[1]) for l in open(SIZES) if l.split()[0].startswith("Chr") and l.split()[0][3:].isdigit())

        rows = []
        for name, bed, col in FEATURES:
            src = f"{REG}/{bed}"
            nuc = f"{tmp}/{name}.bed"
            sh(f"awk '$1~/^Chr[1-5]$/' {src} | sort -k1,1 -k2,2n > {nuc}")
            nfeat = nlines(nuc)
            bp = int(sh(f"bedtools merge -i {nuc} | awk '{{s+=$3-$2}}END{{print s+0}}'").strip())
            frac = bp / G
            hs_in = int(sh(f"bedtools intersect -a {hs} -b {nuc} -u | wc -l").strip())
            feat_hit = int(sh(f"bedtools intersect -a {nuc} -b {hs} -u | wc -l").strip())
            rows.append(dict(feature=name, colour=col, genome_frac=frac, n_feature=nfeat,
                             pct_hotspots_in=hs_in/NH, occupancy_fold=(hs_in/NH)/frac,
                             n_feat_with_hotspot=feat_hit, pct_feat_hotspot=100*feat_hit/nfeat))
    out = pd.DataFrame(rows)
    out.to_csv(f"{ANA}/hotspot_distribution.csv", index=False, float_format="%.5g")
    print(f"genuine (non-rDNA) robust hotspots: {NH}  | compartments: {comp}")
    print(out[["feature","genome_frac","pct_hotspots_in","occupancy_fold","pct_feat_hotspot"]].to_string(index=False))

    # ---- figure ----
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(14, 4.6))
    # A: compartment distribution of genuine hotspots
    ckeys = ["arm","pericentromere","centromere"]
    cvals = [comp.get(k,0) for k in ckeys]
    a0.bar(range(3), cvals, color=["#1f6feb","#ff7f0e","#888"])
    for i,v in enumerate(cvals): a0.text(i, v, str(v), ha="center", va="bottom", fontweight="bold")
    a0.set_xticks(range(3)); a0.set_xticklabels(["arm","peri","centro"])
    a0.set_ylabel("genuine (non-rDNA) hotspots"); a0.set_title(f"A. Compartment ({NH} hotspots)\n(+ {int(nlines(HS))-1-NH} rDNA excluded)", fontsize=10, fontweight="bold")
    a0.grid(axis="y", alpha=.25)
    # B: occupancy fold
    x = np.arange(len(out))
    a1.bar(x, out.occupancy_fold, color=out.colour, alpha=.85)
    a1.axhline(1, color="k", ls="--", lw=.8)
    for xi,v in zip(x,out.occupancy_fold): a1.text(xi, v, f"{v:.2f}×", ha="center", va="bottom", fontsize=9, fontweight="bold")
    a1.set_yscale("log"); a1.set_xticks(x); a1.set_xticklabels(out.feature, rotation=20)
    a1.set_ylabel("occupancy fold  (%hotspots / %genome)")
    a1.set_title("B. Where hotspots fall\n>1 enriched, <1 depleted", fontsize=10, fontweight="bold")
    a1.grid(axis="y", alpha=.25)
    # C: tail metric — % of each feature class that is a hotspot
    a2.bar(x, out.pct_feat_hotspot, color=out.colour, alpha=.85)
    for xi,v in zip(x,out.pct_feat_hotspot): a2.text(xi, v, f"{v:.2f}%", ha="center", va="bottom", fontsize=9, fontweight="bold")
    a2.set_xticks(x); a2.set_xticklabels(out.feature, rotation=20)
    a2.set_ylabel("% of features hosting a hotspot")
    a2.set_title("C. Tail: how often a feature is a hotspot\n(TEs ~15× more often than genes)", fontsize=10, fontweight="bold")
    a2.grid(axis="y", alpha=.25)
    fig.suptitle("DSB hotspot distribution across features — occupancy/tail view "
                 "(complements the per-feature density test)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.94])
    png = f"{FIG}/hotspot_distribution.png"
    fig.savefig(png, dpi=140, bbox_inches="tight"); print("wrote", png)

if __name__ == "__main__":
    main()
