#!/usr/bin/env python3
"""
sBLISS DSB CEN178 metaplot — replicates the spo11-mnaseq-project approach
(scripts/centromere/19_cen178_metaplot.py): anchor at CEN180 repeat START
(strand-aware), +/-1068 bp window, per-window linear detrend, per-chromosome
average then average across chromosomes. Overlays sBLISS DSB density with the
saved CENH3 / MNase / H1 nucleosome metaplots from that project.

Signal = sBLISS break 5' density per bp (break.bedgraph, all reads — multimappers
placed at one position each, which recovers within-monomer phase across 66,131 copies).
"""
import os, numpy as np
from collections import defaultdict, Counter
from scipy.signal import savgol_filter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

BASE   = "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
SPO11  = "/mnt/ssd-8tb/spo11-mnaseq-project"
HITS   = f"{SPO11}/centromere/results/cen180_hits_colcen.bed"   # 66131, Chr1-5, +strand col4
NUCDIR = f"{SPO11}/centromere/results/metaplot"
FAI    = "/mnt/ssd-4tb/crisanto_project/genome/TAIR12/GCA_028009825.2_Col-CC_genomic_withorganelles.fna.fai"
OUTFIG = f"{BASE}/results_TAIR12/analysis/figures"
OUTDIR = f"{BASE}/results_TAIR12/analysis/metaplot"
os.makedirs(OUTDIR, exist_ok=True)

HALF_WIN, MONOMER = 1068, 178
WIN_LEN = 2*HALF_WIN + 1
CHROMS  = ["Chr1","Chr2","Chr3","Chr4","Chr5"]
SAMPLES = ["BA1","BA2","old_BA1_BA2"]
SCOL    = {"BA1":"#1f77b4","BA2":"#2ca02c","old_BA1_BA2":"#d62728"}

chrom_lens = {}
for ln in open(FAI):
    p = ln.split("\t");  chrom_lens[p[0]] = int(p[1])

# CEN180 hits grouped by chr
hits_by_chr = defaultdict(list)
for ln in open(HITS):
    c,s,e,strand = ln.split()[:4]
    hits_by_chr[c].append((int(s),int(e),strand))
print("CEN180 hits:", sum(len(v) for v in hits_by_chr.values()), dict(Counter({k:len(v) for k,v in hits_by_chr.items()})))

def metaplot_from_depth(depth_by_chr):
    """Ian-style: per-window detrend, per-chr mean, then mean across chr."""
    chr_profiles = []
    for chrom, hh in hits_by_chr.items():
        if chrom not in depth_by_chr: continue
        depth = depth_by_chr[chrom]; clen = chrom_lens[chrom]
        acc = np.zeros(WIN_LEN); n = 0
        for (hs, he, strand) in hh:
            anchor = hs if strand == "+" else he
            ws, we = anchor-HALF_WIN, anchor+HALF_WIN+1
            if ws < 0 or we > clen: continue
            w = depth[ws:we].copy()
            if w.sum() < 1: continue
            if strand == "-": w = w[::-1]
            x = np.arange(WIN_LEN)
            sl, ic = np.polyfit(x, w, 1); w = w - (sl*x + ic)   # detrend
            acc += w; n += 1
        if n: chr_profiles.append(acc/n)
    return np.mean(chr_profiles, axis=0) if chr_profiles else None

# build sBLISS break-depth arrays + metaplot per sample
bliss = {}
for s in SAMPLES:
    bg = f"{BASE}/results_TAIR12/breaks/{s}.break.bedgraph"
    total = 0; depth_by_chr = {c: np.zeros(chrom_lens[c]) for c in CHROMS}
    for ln in open(bg):
        c,st,en,v = ln.split("\t")
        if c in depth_by_chr:
            depth_by_chr[c][int(st)] += float(v); total += float(v)
    for c in depth_by_chr: depth_by_chr[c] *= 1e6/total   # CPM (per total breaks)
    prof = metaplot_from_depth(depth_by_chr)
    np.save(f"{OUTDIR}/bliss_cen178_{s}.npy", prof)
    bliss[s] = prof
    print(f"{s}: total breaks={total:.0f}, profile done")

# load nucleosome reference profiles (already detrended CPM, same anchor/window)
nuc = {}
for nm,f in [("CENH3 ChIP","metaplot_CENH3_ChIP.npy"),
             ("MNase-seq","metaplot_MNase-seq.npy"),
             ("H1 ChIP","metaplot_H1_ChIP.npy")]:
    p = f"{NUCDIR}/{f}"
    if os.path.isfile(p): nuc[nm] = np.load(p)
NCOL = {"CENH3 ChIP":"#C0392B","MNase-seq":"#2471A3","H1 ChIP":"#8E44AD"}

x = np.arange(-HALF_WIN, HALF_WIN+1)
sm = lambda a,w=21: savgol_filter(a, w, 3)
def zn(a, mask):           # normalise by std within mask (sigma units)
    sub=a[mask]; sd=sub.std(); return sub/sd if sd>0 else sub
mono = np.arange(-6*MONOMER, 6*MONOMER+1, MONOMER)

# ── Fig 1: sBLISS DSB metaplot (3 samples), full window ──────────────────────
fig, ax = plt.subplots(figsize=(14,5))
for s in SAMPLES:
    if bliss[s] is not None: ax.plot(x, sm(bliss[s]), color=SCOL[s], lw=1.6, label=s)
for v in mono: ax.axvline(v, color="#bbb", lw=0.5, ls="--", zorder=0)
ax.axvline(0, color="#333", lw=1.0)
ax.xaxis.set_major_locator(MultipleLocator(178)); ax.xaxis.set_minor_locator(MultipleLocator(89))
ax.set_xlabel("Distance from CEN180 repeat start (bp)"); ax.set_ylabel("sBLISS DSB density (detrended CPM)")
ax.set_title("sBLISS DSB metaplot over CEN180 (66,131 copies, anchor=start, ±1068 bp, detrended)\n"
             "method = spo11-mnaseq-project script 19", fontweight="bold", fontsize=10)
ax.legend(); ax.spines[["top","right"]].set_visible(False)
plt.tight_layout(); plt.savefig(f"{OUTFIG}/cen178_metaplot_bliss.png", dpi=200, bbox_inches="tight"); plt.close()
print("Saved cen178_metaplot_bliss.png")

# ── Fig 2: sBLISS (BA1) vs CENH3 / MNase / H1, zoomed ±534 (phase comparison) ──
zb = 534; zmask = (x>=-zb)&(x<=zb); xz = x[zmask]
monoz = np.arange(-3*MONOMER, 3*MONOMER+1, MONOMER)
ref = "BA1"
fig, axes = plt.subplots(len(nuc), 1, figsize=(12, 3.2*len(nuc)), gridspec_kw={"hspace":0.4})
if len(nuc)==1: axes=[axes]
fig.suptitle(f"sBLISS DSBs ({ref}) vs centromeric chromatin on CEN180 (±534 bp = 3 monomers, σ-normalised)",
             fontweight="bold", fontsize=11)
for ax,(nm,prof) in zip(axes, nuc.items()):
    ax2 = ax.twinx()
    ax.plot(xz, sm(zn(bliss[ref], zmask),11), color=SCOL[ref], lw=1.8, label=f"sBLISS {ref}")
    ax2.plot(xz, sm(zn(prof, zmask),11), color=NCOL[nm], lw=1.5, label=nm)
    for v in monoz: ax.axvline(v, color="#ccc", lw=0.5, ls="--", zorder=0)
    ax.axvline(0, color="#555", lw=0.8)
    ax.set_title(f"sBLISS DSB  vs  {nm}", fontsize=9, fontweight="bold", loc="left")
    ax.set_ylabel("sBLISS (σ)", color=SCOL[ref], fontsize=8)
    ax2.set_ylabel(nm+" (σ)", color=NCOL[nm], fontsize=8)
    ax.tick_params(axis="y", labelcolor=SCOL[ref]); ax2.tick_params(axis="y", labelcolor=NCOL[nm])
    ax.xaxis.set_major_locator(MultipleLocator(178)); ax.spines[["top"]].set_visible(False)
axes[-1].set_xlabel("Distance from CEN180 repeat start (bp)")
plt.savefig(f"{OUTFIG}/cen178_metaplot_bliss_vs_chromatin.png", dpi=200, bbox_inches="tight"); plt.close()
print("Saved cen178_metaplot_bliss_vs_chromatin.png")

# ── phase quantification: cross-correlation sBLISS(BA1) vs CENH3 within zoom ──
if "CENH3 ChIP" in nuc:
    a = sm(zn(bliss["BA1"], zmask),11); b = sm(zn(nuc["CENH3 ChIP"], zmask),11)
    a-=a.mean(); b-=b.mean()
    lags = np.arange(-MONOMER, MONOMER+1)
    cc = [np.corrcoef(a, np.roll(b, L))[0,1] for L in lags]
    best = lags[int(np.argmax(np.abs(cc)))]; r = cc[int(np.argmax(np.abs(cc)))]
    print(f"\nsBLISS vs CENH3: peak |corr| {r:+.3f} at lag {best} bp "
          f"({'IN-PHASE' if abs(best)<MONOMER/4 else 'ANTI-PHASE' if abs(abs(best)-MONOMER/2)<MONOMER/4 else 'offset'})")
print("Done.")
