#!/usr/bin/env python3
"""
gDNA-normalized BLISS DSB phasing on the CEN180 dimer — controls for
within-monomer mappability/sequence bias (the artefact in the raw dimer approach).
Normalised DSB = (BLISS 5' / mean) / (gDNA coverage / mean). If the nucleosome-
phased peak survives normalisation it is biological; if it flattens it was technical.
"""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

P     = "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
SPO11 = "/mnt/ssd-8tb/spo11-mnaseq-project"
DIMER = f"{P}/results_TAIR12/analysis/cen180_dimer"
MNASE = f"{SPO11}/centromere/results/phasing"
OUT   = f"{P}/results_TAIR12/analysis/figures"
LEN, MONOMER = 356, 178
SAMPLES=["BA1","BA2","old_BA1_BA2"]; SCOL={"BA1":"#1f77b4","BA2":"#2ca02c","old_BA1_BA2":"#d62728"}
CEN180_SEQ=("AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCC"
 "AAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGG"
 "AGTCGAATATGACTTGATGTCATGTGTATGATTG")*2

def norm(a): a=np.asarray(a,float); return a/a.mean() if a.mean()>0 else a
def load2(path):  # pos,count
    d=pd.read_csv(path,sep="\t",header=None).sort_values(0); return d.iloc[:,1].values.astype(float)
def load_depth3(path):  # chrom,pos,depth
    d=pd.read_csv(path,sep="\t",header=None); a=np.zeros(LEN)
    for _,r in d.iterrows():
        p=int(r[1])-1
        if 0<=p<LEN: a[p]=r[2]
    return a
def smooth(a,w=7): return np.convolve(a,np.ones(w)/w,mode="same")
def gc(seq,win=10): return np.array([(seq[i:i+win].count("G")+seq[i:i+win].count("C"))/win for i in range(len(seq)-win+1)])

GDNA_FILE = os.environ.get("BLISS_GDNA", f"{DIMER}/gdna_cen180_depth.tsv")
GDNA_LABEL= os.environ.get("BLISS_GDNA_LABEL", "MNase-input gDNA")
TAG       = os.environ.get("BLISS_TAG", "gdnanorm")
pos=np.arange(1,LEN+1)
gdna=norm(load_depth3(GDNA_FILE))
gdna_safe=np.where(gdna>0.05, gdna, np.nan)   # avoid divide-by-tiny
bliss={s:norm(load2(f"{DIMER}/{s}_cen180_se_r1_5prime.tsv")) for s in SAMPLES}
mn=glob.glob(f"{MNASE}/*_cen180_bulk_depth.tsv"); mnase=norm(np.mean([load2(f) for f in mn],axis=0)) if mn else None
gcv=gc(CEN180_SEQ)

fig=plt.figure(figsize=(12,9)); gs=fig.add_gridspec(3,1,height_ratios=[4,4,1],hspace=0.12)
axA=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1],sharex=axA); axg=fig.add_subplot(gs[2],sharex=axA)
fig.suptitle("sBLISS DSB phasing on CEN180 dimer — raw vs gDNA-normalised\n"
             f"(baseline = {GDNA_LABEL}; controls mappability/sequence bias)", fontsize=11, fontweight="bold")

# A: raw BLISS 5' + gDNA baseline
for s in SAMPLES: axA.plot(pos, smooth(bliss[s]), color=SCOL[s], lw=1.3, label=f"sBLISS DSB 5′ {s}")
axA.plot(pos, smooth(gdna), color="#999", lw=2.0, ls=":", label=f"{GDNA_LABEL} coverage (baseline)")
if mnase is not None: axA.plot(pos, mnase, color="#000", lw=1.5, alpha=0.4, label="MNase occupancy")
axA.axvline(MONOMER,color="k",lw=0.8,ls="--",alpha=0.5); axA.set_ylabel("ratio to mean",fontsize=9)
axA.set_title("A) Raw DSB 5′ density, gDNA baseline, MNase occupancy", fontsize=9, loc="left", fontweight="bold")
axA.legend(fontsize=7,loc="upper right",ncol=2); axA.spines[["top"]].set_visible(False)
plt.setp(axA.get_xticklabels(),visible=False)

# B: gDNA-normalised BLISS (mappability-corrected)
for s in SAMPLES:
    r=bliss[s]/gdna_safe; r=r/np.nanmean(r)
    axB.plot(pos, smooth(np.nan_to_num(r,nan=np.nanmean(r))), color=SCOL[s], lw=1.4, label=f"{s} DSB/gDNA")
if mnase is not None:
    ax2=axB.twinx(); ax2.plot(pos,mnase,color="#000",lw=1.5,alpha=0.35,label="MNase"); ax2.set_ylabel("MNase",fontsize=7,color="#555"); ax2.tick_params(labelcolor="#555",labelsize=6)
axB.axvline(MONOMER,color="k",lw=0.8,ls="--",alpha=0.5); axB.set_ylabel("DSB/gDNA (ratio to mean)",fontsize=9)
axB.set_title("B) gDNA-normalised DSB density (mappability-corrected)", fontsize=9, loc="left", fontweight="bold")
axB.legend(fontsize=7,loc="upper right"); axB.spines[["top"]].set_visible(False)
plt.setp(axB.get_xticklabels(),visible=False)

axg.fill_between(np.arange(1,len(gcv)+1),gcv,color="#888",alpha=0.5,lw=0); axg.axhline(0.5,color="k",lw=0.3,ls="--",alpha=0.4)
axg.axvline(MONOMER,color="k",lw=0.8,ls="--",alpha=0.5); axg.set_ylim(0,1); axg.set_xlim(1,LEN)
axg.set_ylabel("GC",fontsize=7); axg.set_xlabel("Position on CEN180 dimer (bp)",fontsize=9)
axg.spines[["top","right"]].set_visible(False); axg.tick_params(labelsize=6)
axg.xaxis.set_major_locator(MultipleLocator(50)); axg.xaxis.set_minor_locator(MultipleLocator(10))

plt.savefig(f"{OUT}/cen180_phasing_main_bliss_{TAG}.png", dpi=200, bbox_inches="tight")
print(f"Saved cen180_phasing_main_bliss_{TAG}.png")
# correlations on monomer1 after normalisation
if mnase is not None:
    m=slice(0,MONOMER)
    for s in SAMPLES:
        r=bliss[s]/gdna_safe
        rr=np.corrcoef(smooth(np.nan_to_num(r,nan=np.nanmean(r)))[m], mnase[m])[0,1]
        raw=np.corrcoef(smooth(bliss[s])[m], mnase[m])[0,1]
        print(f"  {s}: corr(DSB,MNase) raw={raw:+.3f}  gDNA-norm={rr:+.3f}")
