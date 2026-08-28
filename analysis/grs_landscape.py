#!/usr/bin/env python3
# ponytail: compact congress landscape. Reuses existing CPM bigwigs + region beds.
import pyBigWig, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE="/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12/analysis"
BW=f"{BASE}/bigwig_cpm"; REG=f"{BASE}/regions"
OUT="/mnt/ssd-8tb/HIFI_NAMIL/GRS_presentation/DSB_landscape.pdf"
CHR=[("Chr1",32637894),("Chr2",27782540),("Chr3",26149336),("Chr4",25768865),("Chr5",30142610)]
BIN=100_000
WT="BA1_WT.cpm.bw"; OX=["BA2_OX_40.cpm.bw","BA2_OX_80_1.cpm.bw","BA2_OX_80_2.cpm.bw"]

def load(bwfile, c, L):
    bw=pyBigWig.open(f"{BW}/{bwfile}"); n=L//BIN
    v=np.array(bw.stats(c,0,n*BIN,nBins=n,type="mean"),dtype=float); bw.close()
    return np.nan_to_num(v)
def track(files, c, L):
    return np.mean([load(f,c,L) for f in files],axis=0)

# --- concatenated genome coords ---
off={}; x=0
for c,L in CHR: off[c]=x; x+=L//BIN
TOT=x

def beds(name):  # -> list (chr,s,e)
    out=[]
    for ln in open(f"{REG}/{name}"):
        f=ln.split()
        if f[0] in off: out.append((f[0],int(f[1]),int(f[2])))
    return out
def span_bins(regs):
    xs=[]
    for c,s,e in regs: xs.append((off[c]+s//BIN, off[c]+e//BIN))
    return xs

wt=np.concatenate([track([WT],c,L) for c,L in CHR])
ox=np.concatenate([track(OX,c,L) for c,L in CHR])
# ponytail: median-scale OX to WT on arm background so baselines coincide (library offset, not biology)
bg=(wt<0.25)&(ox<0.25)&(wt>0.05)
ox*=np.median(wt[bg])/np.median(ox[bg])

fig,(ax,ax2)=plt.subplots(2,1,figsize=(13,6.2),gridspec_kw={"height_ratios":[3,1.4],"hspace":0.55})

xa=np.arange(TOT)
YTOP=np.percentile(np.concatenate([wt,ox]),99.7); YBOT=0.15
ax.set_xlim(0,TOT); ax.set_ylim(YBOT,YTOP)
ax.plot(xa,wt,color="#1f6feb",lw=.7,alpha=.9,label="WT")
ax.plot(xa,ox,color="#d62728",lw=.7,alpha=.85,label="CENH3-OX")
# centromere shading
for s,e in span_bins(beds("centromere.bed")):
    ax.axvspan(s,e,color="#888",alpha=.18,lw=0)
# NOR / rDNA shading + label (spikes clipped at top -> arrow)
for s,e in span_bins(beds("nor_45s.bed")):
    ax.axvspan(s,e,color="#2ca02c",alpha=.18,lw=0)
    ax.annotate("rDNA",xy=((s+e)/2,YTOP),xytext=((s+e)/2,YTOP*.86),ha="center",
                color="#2ca02c",fontweight="bold",fontsize=12,
                arrowprops=dict(arrowstyle="-|>",color="#2ca02c",lw=1.5))
# chr boundaries + labels
for c,L in CHR:
    ax.axvline(off[c],color="k",lw=.4,alpha=.3)
    ax.text(off[c]+(L//BIN)/2,YBOT+(YTOP-YBOT)*.02,c,ha="center",va="bottom",fontsize=10,color="#333")
ax.set_ylabel("DSB density\n(CPM / 100 kb)",fontsize=11)
ax.set_xticks([]); ax.spines[["top","right","bottom"]].set_visible(False)
ax.set_title("Genome-wide DSB landscape — flat across arms & centromeres, spikes at rDNA",fontsize=13,fontweight="bold")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],color="#1f6feb",lw=1.5,label="WT"),
                   Line2D([],[],color="#d62728",lw=1.5,label="CENH3-OX"),
                   Patch(fc="#2ca02c",alpha=.3,label="rDNA (45S NOR)"),
                   Patch(fc="#888",alpha=.3,label="centromere")],
          loc="upper right",frameon=False,fontsize=10,ncol=2)

# --- enrichment bar: per-Mb DSB fold, multimapper-inclusive (BA1.mm2raw) ---
# ponytail: precomputed via bedtools intersect. self = ÷genome-mean DSB;
# gDNA = additionally ÷WGA coverage fold (copy-number corrected, wga.cov.cpm.bw).
names =["rDNA","TEs","genes","arms","centro"]
self_ =[4.33, 1.26, 0.90, 0.74, 0.25]
gdna  =[6.48, 1.31, 0.87, 0.71, 0.32]
col   =["#2ca02c","#9467bd","#8c564b","#1f6feb","#888"]
x=np.arange(len(names)); w=0.38
ax2.bar(x-w/2,self_,w,color=col,alpha=.5,label="÷ genome mean")
ax2.bar(x+w/2,gdna ,w,color=col,label="÷ gDNA (copy-number corr.)")
ax2.axhline(1,color="k",lw=.7,ls="--")
ax2.set_xticks(x); ax2.set_xticklabels(names)
ax2.set_ylabel("DSB / Mb fold",fontsize=10)
ax2.set_yscale("log"); ax2.set_ylim(0.2,max(gdna)*1.7)
ax2.set_yticks([0.25,0.5,1,2,4]); ax2.set_yticklabels(["0.25","0.5","1","2","4"])
for xi,v in zip(x-w/2,self_): ax2.text(xi,v*1.07,f"{v:.1f}",ha="center",fontsize=8,color="#555")
for xi,v in zip(x+w/2,gdna ): ax2.text(xi,v*1.07,f"{v:.1f}×",ha="center",fontsize=9,fontweight="bold")
ax2.legend(loc="upper right",frameon=False,fontsize=8,ncol=2)
ax2.spines[["top","right"]].set_visible(False); ax2.tick_params(labelsize=11)
ax2.set_title("DSB enrichment per Mb (multimapper-inclusive): rDNA & transposons up, arms/centromere down",fontsize=11)

fig.savefig(OUT,dpi=200,bbox_inches="tight")
print("wrote",OUT,"| gDNA-corr fold:",dict(zip(names,gdna)))
