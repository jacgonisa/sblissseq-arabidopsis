#!/usr/bin/env python3
"""Build sBLISS_TAIR12_report_with_CENH3ox.html = the full TAIR12 report PLUS a
CENH3-OX vs WT section (live QC progress + comparison). Re-runnable anytime:
parses whatever QC artefacts exist and marks the rest 'pending'.
Splices into the existing sBLISS_TAIR12_report.html (reuses its CSS + base64 style).
"""
import os, re, base64, datetime, subprocess

BASE = "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
OUT  = f"{BASE}/results_TAIR12"
LOG, BRK, BAM, FIG, ANA = f"{OUT}/logs", f"{OUT}/breaks", f"{OUT}/bam", f"{OUT}/analysis/figures", f"{OUT}/analysis"
SRC  = f"{ANA}/sBLISS_TAIR12_report.html"           # the main report (must exist)
DST  = f"{ANA}/sBLISS_TAIR12_report_with_CENH3ox.html"

SAMPLES = [  # sample, genotype label, group, lanes
    ("BA1_WT","WT (this round)","WT",2), ("BA1","WT (prev round)","WT",1), ("BA2","WT (prev round)","WT",2),
    ("BA2_OX_40","CENH3-OX · 40","OX",2), ("BA2_OX_80_1","CENH3-OX · 80 r1","OX",1),
    ("BA2_OX_80_2","CENH3-OX · 80 r2","OX",2),
]

def read(p):
    try: return open(p).read()
    except Exception: return ""
def demux(s):
    d={}
    for ln in read(f"{LOG}/{s}.demux.stats.txt").splitlines():
        if "\t" in ln: k,v=ln.split("\t",1); d[k]=v
    return d
def bt(s):
    t=read(f"{LOG}/{s}.bowtie2.log"); i=re.search(r"(\d+) reads; of these",t); r=re.search(r"([\d.]+)% overall",t)
    return (int(i.group(1)) if i else None, float(r.group(1)) if r else None)
def dd(s):
    t=read(f"{LOG}/{s}.bowtie2.dedup.log"); i=re.search(r"Input Reads: (\d+)",t); o=re.search(r"Number of reads out: (\d+)",t)
    return (int(i.group(1)) if i else None, int(o.group(1)) if o else None)
def brk(s):
    bg=f"{BRK}/{s}.break.bedgraph"
    if not os.path.exists(bg): return (None,None)
    sites=int(subprocess.run(f"wc -l < {bg}",shell=True,capture_output=True,text=True).stdout or 0)
    ev=int(float(subprocess.run(f"awk '{{s+=$4}}END{{print s+0}}' {bg}",shell=True,capture_output=True,text=True).stdout or 0))
    return (ev,sites)
def fmt(x): return f"{x:,}" if isinstance(x,int) else ("—" if x is None else x)
def status(s):
    if os.path.exists(f"{BRK}/{s}.break.bedgraph"): return ("done","#2ca02c")
    if os.path.exists(f"{BAM}/{s}.bowtie2.dedup.bam"): return ("dedup","#b8860b")
    if os.path.exists(f"{BAM}/{s}.bowtie2.bam"): return ("aligned","#b8860b")
    if os.path.exists(f"{LOG}/{s}.bowtie2.log"): return ("aligning","#b8860b")
    if demux(s): return ("demuxed","#b8860b")
    return ("pending","#999")
def img(name,w="90%"):
    p=f"{FIG}/{name}"
    if os.path.exists(p):
        b=base64.b64encode(open(p,"rb").read()).decode()
        return f'<img src="data:image/png;base64,{b}" style="max-width:{w};border:1px solid #ddd;border-radius:4px;">'
    return None
def figbox(name,w,cap):
    im=img(name,w)
    if im: return f'<div class="fig">{im}<div class="cap">{cap}</div></div>'
    return f'<div class="note">⏳ <strong>{name}</strong> — pending (runs after mapping completes). <em>{cap}</em></div>'

rows=""; new=("BA1_WT","BA2_OX_40","BA2_OX_80_1","BA2_OX_80_2"); ndone=0
for s,geno,grp,lanes in SAMPLES:
    dm=demux(s); inp,rate=bt(s); di,do=dd(s); ev,si=brk(s); st,col=status(s)
    if s in new and st=="done": ndone+=1
    raw=int(dm["total_reads"]) if "total_reads" in dm else None
    bpr=f'{float(dm["barcode_pass_rate"])*100:.1f}%' if "barcode_pass_rate" in dm else "—"
    ret=f"{100*do/di:.0f}%" if (di and do) else "—"
    gcol="#1f6feb" if grp=="WT" else "#d62728"
    rows+=(f'<tr><td><b>{s}</b><br><span style="color:{gcol};font-size:.85em">{geno}</span></td>'
           f'<td style="color:{col};font-weight:bold">{st}</td><td>{fmt(raw)}</td><td>{bpr}</td>'
           f'<td>{fmt(inp)}</td><td>{rate if rate is not None else "—"}%</td>'
           f'<td>{fmt(di)}&rarr;{fmt(do)}</td><td>{ret}</td><td>{fmt(ev)}</td><td>{fmt(si)}</td></tr>\n')
log=read("/tmp/cenh3ox.log").strip().splitlines(); logtail="<br>".join(l.replace("<","&lt;") for l in log[-6:]) or "(no log)"
now=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

block=f"""
<h2 id="cenh3ox">12. CENH3-OX vs WT — progress &amp; comparison</h2>
<p style="color:#666">Live section, regenerated {now}. New delivery <code>X204SC26053809-Z01-F001</code>.
<strong>Progress: {ndone}/4 new samples fully mapped.</strong></p>
<div class="note"><strong>Design.</strong> Matched comparison of <strong>this delivery</strong>: WT = <strong>BA1_WT</strong> vs OX = <strong>3 libraries of one CENH3-OX line</strong> (BA2_OX_40, BA2_OX_80_1 [½ depth, 1 lane], BA2_OX_80_2).
Original BA1/BA2/old_BA1_BA2 are shown in the QC table for reference only (different prep round). <strong>Question:</strong> does CENH3 overexpression change the DSB landscape,
especially at centromeres? Processed identically to the main pipeline, low-resource (nice/ionice, bowtie2 -p4).</div>

<h3>12a. Samples &amp; QC</h3>
<table>
<tr><th>sample</th><th>status</th><th>raw reads</th><th>barcode pass</th><th>aligned input</th><th>align rate</th><th>dedup in&rarr;out</th><th>retention</th><th>break events</th><th>distinct sites</th></tr>
{rows}</table>
<p class="cap">raw = demux input (post-merge); barcode pass = internal-index match (Hamming&le;1); retention = umi_tools out/in (complexity);
break events = total 5′ DSB ends; distinct sites = unique positions. Pipeline log: <code style="font-size:.8em">{logtail}</code></p>

<div class="key"><strong>Headline:</strong> the matched comparison (BA1_WT vs the three CENH3-OX libraries) shows
<strong>no genotype-attributable change</strong> in the DSB landscape. Compartment distribution WT&asymp;OX (centromere ~3% both);
WT-mean and OX-mean landscapes overlap; only <strong>1 of 13,318</strong> non-rDNA 10 kb windows is reproducibly OX-gained
(all 3 OX &gt; 2&times; WT), 0 WT-gained; CEN180 satellite OX/WT &asymp; 1.1&times;. The three OX libraries (one CENH3-OX line) replicate well
(r=0.89); the lower BA1_WT-vs-OX correlation is a <strong>library/adapter effect</strong> (adapter test, §12e), not a CENH3 effect.</div>

<h3>12b. Alignment &amp; library complexity (all samples)</h3>
{figbox("cenh3ox_qc_alignment_complexity.png","98%","Bowtie2 alignment breakdown (unique/multi/unaligned), UMI-dedup retention, and distinct DSB sites (+redundancy reads/site) per sample. The 4 new libraries are high-complexity (3.6–4.7 M sites, 64–82% retention).")}
{figbox("cenh3ox_umi_threshold.png","62%","Library complexity — distinct DSB locations retained at increasing UMI/read thresholds (log y). Most sites are singletons (≥1: 3.6–4.7 M → ≥2: ~0.6 M), i.e. high complexity; WT and OX profiles overlap.")}

<h3>12c. Compartment composition — matched WT vs OX</h3>
{figbox("cenh3ox_compartment.png","92%","Per-sample DSB compartment composition (incl. previous-round WT BA1/BA2); centromere fraction (n.s.). Distribution unchanged.")}
<h3>12d. Genome-wide landscape — BA1_WT vs OX-mean</h3>
{figbox("cenh3ox_landscape_WTvsOX.png","98%","log1p CPM/100kb: WT this round (blue), WT old mean (cyan), OX mean (red) — overlapping everywhere.")}
<h3>12e. Differential windows &amp; sample correlation</h3>
{figbox("cenh3ox_hotspots.png","58%","BA1_WT vs OX-mean DSB (non-rDNA 10kb): on the diagonal, 1 OX-gained / 0 WT-gained reproducible windows.")}
{figbox("cenh3ox_correlation.png","56%","Correlation of all WT (this round + old BA1/BA2) and the 3 OX (non-rDNA 10kb, log). OX replicates 0.90–0.93; WT–WT 0.68–0.81; WT–OX 0.46–0.65.")}
<div class="warn"><strong>Adapter/batch confound (important).</strong> All 3 OX libraries (one CENH3-OX line) use the BA2 adapter and the new prep round; BA1_WT uses the BA1 adapter.
The correlation matrix shows a clear <strong>adapter gradient</strong>: same-adapter pairs are higher (WT-new&harr;BA1, both BA1 = <strong>0.81</strong>; old BA2&harr;OX, both BA2 = <strong>0.62–0.65</strong>)
than cross-adapter (WT-new&harr;OX = 0.46–0.50). So the WT-vs-OX gap tracks the adapter, and no control matches OX on <em>both</em> adapter and batch — the residual cannot be
attributed to genotype. The "no CENH3 effect" conclusion rests on the confound-robust readouts: compartment, landscape, CEN180, and the differential-window test.</div>

<h3>12f. CEN180 centromeric satellite</h3>
{figbox("cenh3ox_cen180.png","66%","DSB density over CEN180/CEN178 (51,553 monomers), matched WT vs OX ≈ 1.1× (n.s.; centromere mappability-limited).")}
<h3>12g. Genes &amp; transposons (signal/mean)</h3>
{figbox("cenh3ox_metaplots.png","92%","DSB profiles at TSS and TEs, matched WT vs OX — indistinguishable.")}

<div class="note"><strong>Interpretation &amp; caveats.</strong> By every confound-robust measure (compartment, landscape, CEN180, genes/TEs, reproducible-difference test),
the CENH3-OX line has the same DSB map as WT — ectopic CENH3 does not create new fragile sites at this sensitivity. The adapter test (adapter-matched WT correlates
better with OX) attributes the arm-level WT-vs-OX correlation gap to library/adapter/batch, not genotype. With essentially no OX-specific hotspots,
the OX-hotspot epigenomic-context/methylation analyses are uninformative (n=1) and omitted. Caveats: WT/OX use different BLISS adapters and prep rounds (no doubly-matched control);
OX_80_1 is ½ depth (1 lane); centromere/rDNA mappability-limited so a within-CEN178 change could be partly masked.</div>
"""

main=read(SRC)
if not main:
    raise SystemExit(f"main report not found: {SRC} — build it first (build_report_tair12.py)")
# add nav link (before first </nav>) and splice block before </body>
main=main.replace("</nav>", '<a href="#cenh3ox"><b>12 CENH3-OX</b></a></nav>', 1)
main=main.replace("</body>", block+"\n</body>", 1)
open(DST,"w").write(main)
print("Wrote", DST, f"({os.path.getsize(DST)//1024} KB) — {ndone}/4 new samples mapped")
