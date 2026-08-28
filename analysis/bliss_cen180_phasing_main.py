#!/usr/bin/env python3
"""
BLISS DSB phasing on the CEN180 dimer — style of spo11 cen180_phasing_main.png.
Maps already done (bliss_cen180_dimer.sh): BLISS 5' DSB profiles on 356bp dimer.
Overlays MNase nucleosome occupancy (bulk) so DSB position can be read relative
to the nucleosome dyad. GC track below.
"""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

P      = "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
SPO11  = "/mnt/ssd-8tb/spo11-mnaseq-project"
DIMER  = f"{P}/results_TAIR12/analysis/cen180_dimer"
MNASE  = f"{SPO11}/centromere/results/phasing"
OUT    = f"{P}/results_TAIR12/analysis/figures"
LEN, MONOMER = 356, 178
SAMPLES = ["BA1","BA2","old_BA1_BA2"]
SCOL = {"BA1":"#1f77b4","BA2":"#2ca02c","old_BA1_BA2":"#d62728"}

CEN180_SEQ = (
 "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCC"
 "AAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGG"
 "AGTCGAATATGACTTGATGTCATGTGTATGATTG"
 "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCC"
 "AAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGG"
 "AGTCGAATATGACTTGATGTCATGTGTATGATTG")

def load_norm(path):
    d = pd.read_csv(path, sep="\t", header=None, names=["pos","v"]).sort_values("pos")
    a = d["v"].values.astype(float)
    return a/a.mean() if a.mean()>0 else a

def avg_norm(pattern):
    fs=glob.glob(pattern); arr=[load_norm(f) for f in fs if len(open(f).readlines())==LEN]
    return np.mean(arr,axis=0) if arr else None

def smooth(a,w=7):
    k=np.ones(w)/w; return np.convolve(a,k,mode="same")

def gc(seq,win=10):
    return np.array([(seq[i:i+win].count("G")+seq[i:i+win].count("C"))/win
                     for i in range(len(seq)-win+1)])

pos = np.arange(1,LEN+1)
bliss = {s: load_norm(f"{DIMER}/{s}_cen180_se_r1_5prime.tsv") for s in SAMPLES}
mnase = avg_norm(f"{MNASE}/*_cen180_bulk_depth.tsv")   # nucleosome occupancy reference
gcv   = gc(CEN180_SEQ); gcx = np.arange(1,len(gcv)+1)

fig = plt.figure(figsize=(12,6.5))
gs = fig.add_gridspec(2,1,height_ratios=[5,1],hspace=0.05)
ax = fig.add_subplot(gs[0]); axg = fig.add_subplot(gs[1], sharex=ax)
fig.suptitle("sBLISS DSB phasing on the CEN180 dimer (2×178 bp consensus)\n"
             "DSB 5′-end density vs MNase nucleosome occupancy — minimap2 -ax sr dimer mapping",
             fontsize=11, fontweight="bold")

for s in SAMPLES:
    ax.plot(pos, smooth(bliss[s]), color=SCOL[s], lw=1.4, alpha=0.9, label=f"sBLISS DSB 5′ — {s}")
ax2 = ax.twinx()
if mnase is not None:
    ax2.plot(pos, mnase, color="#555555", lw=2.2, alpha=0.6, label="MNase occupancy (nucleosome)")
    ax2.set_ylabel("MNase occupancy (ratio to mean)", fontsize=8, color="#555")
    ax2.tick_params(axis="y", labelcolor="#555", labelsize=7)
ax.axvline(MONOMER, color="black", lw=0.8, ls="--", alpha=0.5)
ax.set_ylabel("sBLISS DSB 5′ density (ratio to mean)", fontsize=9)
ax.set_xlim(1,LEN); ax.spines[["top"]].set_visible(False)
ax.xaxis.set_major_locator(MultipleLocator(50)); ax.xaxis.set_minor_locator(MultipleLocator(10))
l1,lab1=ax.get_legend_handles_labels(); l2,lab2=ax2.get_legend_handles_labels()
ax.legend(l1+l2, lab1+lab2, fontsize=8, loc="upper right", framealpha=0.9)
plt.setp(ax.get_xticklabels(), visible=False)

axg.fill_between(gcx, gcv, color="#888", alpha=0.5, lw=0)
axg.axhline(0.5, color="black", lw=0.3, ls="--", alpha=0.4)
axg.axvline(MONOMER, color="black", lw=0.8, ls="--", alpha=0.5)
axg.set_ylabel("GC", fontsize=7); axg.set_ylim(0,1); axg.set_xlim(1,LEN)
axg.set_xlabel("Position on CEN180 dimer (bp)", fontsize=9)
axg.spines[["top","right"]].set_visible(False); axg.tick_params(labelsize=6)

plt.savefig(f"{OUT}/cen180_phasing_main_bliss.png", dpi=200, bbox_inches="tight")
plt.savefig(f"{OUT}/cen180_phasing_main_bliss.pdf", bbox_inches="tight")
print("Saved cen180_phasing_main_bliss.png")

# correlation of BLISS DSB vs MNase occupancy on the dimer (monomer 1: 1-178)
if mnase is not None:
    m1 = slice(0,MONOMER)
    for s in SAMPLES:
        r = np.corrcoef(smooth(bliss[s])[m1], mnase[m1])[0,1]
        print(f"  {s}: corr(DSB 5′, MNase occupancy) on monomer = {r:+.3f}")
