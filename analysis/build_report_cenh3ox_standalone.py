#!/usr/bin/env python3
"""STANDALONE CENH3-OX vs WT sBLISS report (all samples, both runs) — TAIR12.
Self-contained HTML; does NOT touch the main report. Re-runnable.
Output: results_TAIR12/analysis/sBLISS_CENH3ox_report.html
"""
import os, re, base64, subprocess, datetime
BASE="/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"; OUT=f"{BASE}/results_TAIR12"
LOG,BRK,HS,FIG=f"{OUT}/logs",f"{OUT}/breaks",f"{OUT}/analysis/hotspots",f"{OUT}/analysis/figures"
DST=f"{OUT}/analysis/sBLISS_CENH3ox_report.html"
# sample, run, adapter, genotype
SAMPLES=[("BA1","run1","BA1","WT"),("BA2","run1","BA2","WT"),("old_BA1_BA2","run1","BA2","WT (mixed)"),
         ("BA1_WT","run2","BA1","WT"),("BA2_OX_40","run2","BA2","CENH3-OX"),
         ("BA2_OX_80_1","run2","BA2","CENH3-OX"),("BA2_OX_80_2","run2","BA2","CENH3-OX")]
def g(p): return open(p).read() if os.path.exists(p) else ""
def sh(c): return subprocess.run(c,shell=True,capture_output=True,text=True).stdout.strip()
# UMI-dedup counts (mapped before dedup, unique after) cached from BAMs
DEDUP={}
for ln in g(f"{OUT}/analysis/dedup_counts.tsv").splitlines()[1:]:
    p=ln.split("\t")
    if len(p)==3: DEDUP[p[0]]=(int(p[1]),int(p[2]))
HOTC={}
for ln in g(f"{HS}/hotspot_counts.tsv").splitlines()[1:]:
    p=ln.split("\t")
    if len(p)==7: HOTC[(p[0],p[1])]=p[2:]
def row(s):
    t=g(f"{LOG}/{s}.bowtie2.log")
    inp=re.search(r"(\d+) reads",t);un=re.search(r"([\d.]+)%\) aligned 0",t);u1=re.search(r"([\d.]+)%\) aligned exactly",t);um=re.search(r"([\d.]+)%\) aligned >1",t)
    dm=g(f"{LOG}/{s}.demux.stats.txt");bc=re.search(r"barcode_pass_rate\t([\d.]+)",dm)
    mapped,dedup=DEDUP.get(s,(None,None))
    ret=100*dedup/mapped if (mapped and dedup) else None
    bg=f"{BRK}/{s}.break.bedgraph"
    ev=float(sh(f"awk '$1~/^Chr[1-5]$/{{x+=$4}}END{{print x+0}}' {bg}") or 0)
    sites=int(sh(f"awk '$1~/^Chr[1-5]$/' {bg}|wc -l") or 0)
    lam=ev/142500.0
    hf=f"{HS}/{s}_1kb.hotspots.bed"
    hc=sh("awk '{t++; if($10==\"TRUE\")n++; r[$9]++} END{printf \"%d %d %d %d %d\",t,n,r[\"pericentromere\"],r[\"arm\"]-n,r[\"centromere\"]}' "+hf).split() if os.path.exists(hf) else None
    return dict(inp=inp.group(1) if inp else "—",un=un.group(1) if un else "—",u1=u1.group(1) if u1 else "—",um=um.group(1) if um else "—",
                bc=f"{float(bc.group(1))*100:.0f}" if bc else "—",ret=f"{ret:.0f}" if ret else "—",sites=sites/1e6,lam=lam,hc=hc,
                mapped=mapped,dedup=dedup)
def img(name,w="90%"):
    p=f"{FIG}/{name}"
    if not os.path.exists(p): return f'<p><em>[missing: {name}]</em></p>'
    b=base64.b64encode(open(p,"rb").read()).decode()
    return f'<div class="fig"><img src="data:image/png;base64,{b}" style="max-width:{w};border:1px solid #ddd;border-radius:4px"></div>'
def hottable(method):
    r=""
    for s,run,adp,geno in SAMPLES:
        v=HOTC.get((s,method))
        if not v: continue
        gc="#1f6feb" if "WT" in geno else "#d62728"
        r+=(f'<tr><td><b>{s}</b> <span style="color:{gc};font-size:.8em">{geno}</span></td><td>{run}</td>'
            f'<td><b>{int(v[0]):,}</b></td><td>{int(v[1]):,}</td><td>{v[2]}</td><td>{v[3]}</td><td>{v[4]}</td></tr>')
    return f'<table><tr><th>sample</th><th>run</th><th>total</th><th>45S rDNA</th><th>pericentromere</th><th>arm</th><th>centromere</th></tr>{r}</table>'

# build master table
align="";hot=""
for s,run,adp,geno in SAMPLES:
    r=row(s);gc="#1f6feb" if "WT" in geno else "#d62728"
    md=f'{r["mapped"]/1e6:.1f} M' if r["mapped"] else "—"; dp=f'{r["dedup"]/1e6:.2f} M' if r["dedup"] else "—"
    align+=(f'<tr><td><b>{s}</b><br><span style="color:{gc};font-size:.85em">{geno}</span></td><td>{run}</td><td>{adp}</td>'
            f'<td>{int(r["inp"]):,}</td><td>{r["bc"]}%</td><td>{r["un"]}%</td><td>{r["u1"]}%</td><td><b>{r["um"]}%</b></td>'
            f'<td>{md}</td><td><b>{dp}</b></td><td>{r["ret"]}%</td><td>{r["sites"]:.2f} M</td></tr>\n')
    if r["hc"] and len(r["hc"])==5:
        tot,nor,peri,arm,cen=r["hc"]
        hot+=(f'<tr><td><b>{s}</b><br><span style="color:{gc};font-size:.85em">{geno}</span></td><td>{run}</td><td>{r["lam"]:.0f}</td>'
              f'<td><b>{int(tot):,}</b></td><td>{int(nor):,}</td><td>{peri}</td><td>{arm}</td><td>{cen}</td></tr>\n')
now=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
css="""body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1080px;margin:1.5em auto;padding:0 2em;color:#222;line-height:1.5}
h1{border-bottom:3px solid #6a1b9a;padding-bottom:.3em}h2{color:#6a1b9a;margin-top:1.6em;border-bottom:1px solid #ddd;padding-bottom:.2em}h3{color:#555}
table{border-collapse:collapse;margin:1em 0;font-size:.88em}th,td{border:1px solid #ccc;padding:5px 9px;text-align:right}th{background:#f3e5f5}td:first-child,th:first-child{text-align:left}
.fig{margin:1.2em 0;text-align:center}.cap{font-size:.86em;color:#666}.note{background:#fff8e1;border-left:4px solid #ffb300;padding:.7em 1.1em;margin:1em 0}
.key{background:#e8f5e9;border-left:4px solid #2ca02c;padding:.7em 1.1em;margin:1em 0}.warn{background:#ffebee;border-left:4px solid #e53935;padding:.7em 1.1em;margin:1em 0}
code{background:#f4f4f4;padding:1px 5px;border-radius:3px}nav{background:#fafafa;border:1px solid #eee;padding:.6em 1.2em}nav a{margin-right:1em;color:#6a1b9a;text-decoration:none}"""
html=f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>sBLISS CENH3-OX vs WT (TAIR12)</title><style>{css}</style></head><body>
<h1>sBLISS-seq — CENH3-OX vs WT (Col-CC / TAIR12)</h1>
<p>Standalone analysis of all 7 sBLISS libraries across both Novogene deliveries. Generated {now}.
<strong>run1</strong> = original (BA1, BA2, old_BA1_BA2); <strong>run2</strong> = this delivery <code>X204SC26053809</code> (BA1_WT + 3 CENH3-OX libraries, one OX line).</p>
<nav><a href="#qc">1 Sequencing &amp; alignment</a><a href="#cx">2 Complexity</a><a href="#un">3 Unaligned reads</a><a href="#hot">4 Hotspots</a><a href="#cmp">5 WT vs OX</a><a href="#dif">5b Differential</a><a href="#concl">6 Conclusion</a></nav>

<h2 id="qc">1. Sequencing &amp; alignment — all samples</h2>
<table>
<tr><th>sample</th><th>run</th><th>adapter</th><th>input reads</th><th>barcode pass</th><th>unaligned</th><th>unique</th><th>multi (&gt;1×)</th><th>mapped</th><th>after UMI-dedup</th><th>retention</th><th>distinct sites</th></tr>
{align}</table>
<div class="note">Unaligned tracks <strong>adapter + prep quality</strong>, not genotype: BA1-adapter libraries (BA1 13.7%, BA1_WT 14.4%) and old_BA1_BA2 (19.5%) are high; run-2 OX (BA2 adapter, cleaner prep) lowest (2.2–2.5%). OX's higher <strong>multimapper</strong> fraction (42–44% vs 30% WT) = more reads from repetitive (rDNA/centromere) DNA.</div>

<h2 id="cx">2. Library complexity</h2>
{img("cenh3ox_umi_threshold.png","62%")}
<p class="cap">Distinct DSB locations retained at increasing UMI/read thresholds (log y). Most sites are singletons (high complexity); WT and OX overlap.</p>
{img("cenh3ox_qc_alignment_complexity.png","98%")}
<p class="cap">Alignment breakdown, UMI-dedup retention, and distinct sites (+ redundancy reads/site) per sample.</p>

<h2 id="un">3. What are the (excess) unaligned reads?</h2>
<div class="note">In BA1_WT the 14.4% unaligned are <strong>BLISS adapter / empty-insert read-through</strong>, not genomic: ~9% carry the adapter backbone
(<code>CGATTGAGGCCGGTAATACGACTCACTATAGGGG…</code> = the <strong>T7 promoter</strong> of the BLISS adapter), and the rest are GC-rich (mean 50% GC vs 38% for
aligned reads), diverse, non-mapping sequences — the signature of empty/short-insert clusters + adapter dimer. This matches the original Novogene QC flag
(BA1 elevated GC +12.7% = adapter carry-through). It does <em>not</em> bias DSB results (CPM is on mapped reads) but marks the WT library as the dirtier prep.</div>

<h2 id="hot">4. Hotspot calling — method, normalisation &amp; results</h2>
<div class="key"><strong>Method (<code>hotspots_brgenomics.R</code>).</strong> Genome in <strong>1 kb windows</strong> → count <strong>raw break events</strong>
(5′ DSB ends, UMI-deduplicated, multimapper-inclusive MAPQ≥0) per window → <strong>Poisson</strong> test vs a flat null <strong>λ = that sample's own genome-mean
events/window</strong> → BH-adjust → call if <strong>padj&lt;0.001 AND fold (count/λ) ≥ 3</strong>. No blacklist; rDNA/centromere flagged.</div>
<p>Because the answer depends entirely on the normalisation, here are <strong>three versions</strong> for all samples:</p>
<h3>4a. MAPQ≥0, flat-Poisson (multimapper-inclusive) — the default</h3>
{hottable("MAPQ0")}
<p class="cap">λ = per-sample genome mean over <em>all</em> reads. rDNA dominates and tracks multimapper content (OX 42–44% multi → ~6,100 rDNA vs WT 30% → 4,750) — a library/alignment effect, not CENH3.</p>
<h3>4b. MAPQ≥10, flat-Poisson (unique reads only)</h3>
{hottable("MAPQ10")}
<p class="cap"><strong>rDNA collapses to ~20–40</strong> (it was multimapper-driven). With the rDNA pileup gone, λ drops, so arm/pericentromere windows now pass — the counts there are <em>not</em> comparable to 4a. The OX-vs-WT gap is small and inconsistent (OX_40 high, OX_80_2 below WT).</p>
<h3>4c. gDNA-normalised (÷ WGA gDNA per window) — the DNA-content-corrected version</h3>
{hottable("gDNA")}
<p class="cap">λ from the local WGA gDNA coverage (corrects mappability + copy-number), call padj&lt;0.001 &amp; fold≥3 &amp; ≥10 breaks. <strong>rDNA returns as genuine</strong> (it is ~12× over gDNA, i.e. a real DSB hotspot, not just copy number). Non-rDNA (arm/peri/cen) are modest and similar across genotypes.</p>
<h3>4d. Matched CENH3-OX <em>input</em>-normalised (÷ the OX line's own gDNA per window)</h3>
{hottable("input")}
<p class="cap">Same Poisson-vs-baseline call but λ from the <strong>matched OX input/gDNA</strong> (this study) instead of the external WGA. <strong>Difference vs flat-Poisson (4a) is the expected one:</strong>
rDNA <em>up</em> (genuine 12× enrichment survives), <strong>pericentromere down ~half</strong> (794→338 — those were coverage/copy-number-driven, not DSB-specific). And it agrees closely with the WGA version (4c: BA1_WT 7,419 vs 7,443), confirming the gDNA correction is robust. WT vs OX still shows no consistent difference.</p>
<div class="warn"><strong>Take-home on normalisation.</strong> The three callers give very different counts from the <em>same</em> data — MAPQ≥0 is multimapper-/depth-confounded (why I'm now showing the others), MAPQ≥10 removes the repeats entirely, and ÷gDNA is the confound-robust one (rDNA genuine, arms modest). Under <em>none</em> of them is there a consistent CENH3-OX vs WT hotspot difference; the apparent gaps track adapter/prep/multimapper content.</div>

<h3>4e. Hotspot overlap (Venn) — gDNA-normalised</h3>
{img("cenh3ox_hotspot_venn.png","94%")}
<p class="cap">BA1_WT vs OX-40 vs OX-80 (rep union). <strong>Left:</strong> all hotspots (rDNA dominates and is shared). <strong>Right:</strong> non-rDNA only —
359 shared by all three, but <strong>WT-only (268) ≈ OX-only (269)</strong>: symmetric per-library variability, <em>not</em> an OX-specific gain.</p>

<h3>4f. Genome-wide hotspot map</h3>
{img("cenh3ox_hotspot_map.png","98%")}
<p class="cap">Each circle = one 1 kb gDNA-normalised hotspot; 7 sample tracks per chromosome, coloured by region (rDNA purple / centromere red / pericentromere orange / arm blue).
The pattern is the same across all samples — dense rDNA in the NORs, pericentromeric flanking the centromeres, scattered on the arms; WT and OX hotspots fall in the same places.</p>

<h2 id="cmp">5. WT vs CENH3-OX comparison</h2>
{img("cenh3ox_genomewide_persample.png","98%")}
<p class="cap"><strong>Genome-wide DSB density, every sample</strong> (log1p CPM/100 kb; Chr1–5 concatenated; grey=centromere, purple=NOR, lines=chrom boundaries).
NOR rDNA spikes in all; arms flat; WT, OX and the run-1 WTs are visually indistinguishable.</p>
{img("cenh3ox_compartment.png","92%")}
<p class="cap">DSB compartment composition per sample — distribution unchanged by CENH3-OX (centromere ~3% in all).</p>
{img("cenh3ox_landscape_allsamples.png","98%")}
<p class="cap">Genome-wide DSB density per chromosome — <strong>all 7 samples overlaid</strong> (log1p CPM/100 kb; grey=cen, purple=NOR). WT and OX overlap everywhere at this scale.</p>
{img("cenh3ox_hotspots.png","56%")}
{img("cenh3ox_correlation.png","56%")}
<p class="cap">Differential windows (1 OX-gained / 0 WT-gained of 13,318) and the sample correlation: OX replicates 0.90–0.93; WT–WT 0.68–0.81; WT–OX 0.46–0.65.</p>
<div class="warn"><strong>Adapter/batch confound.</strong> The 3 OX libraries (one line) use the BA2 adapter + new prep; BA1_WT uses BA1. The correlation shows an adapter gradient
(same-adapter WT-new↔BA1 = 0.81; BA2↔OX = 0.62–0.65; cross-adapter WT-new↔OX = 0.46–0.50). No control matches OX on both adapter and batch, so the WT-vs-OX gap is not attributable to genotype.</div>
{img("cenh3ox_cen180.png","66%")}
{img("cenh3ox_metaplots.png","92%")}
<p class="cap">CEN180 satellite DSB density (≈1.1×, n.s.) and TSS/TE profiles — WT and OX indistinguishable.</p>

<h2 id="dif">5b. Differential enrichment (1 kb windows, edgeR) &amp; MDS</h2>
<p>Formal per-window test (edgeR quasi-likelihood, 1 kb windows, FDR&lt;0.05 &amp; |log2FC|&gt;1) — the rigorous version of the BiCroLab descriptive comparison.</p>
{img("cenh3ox_MDS.png","52%")}
<p class="cap"><strong>MDS of 1 kb window counts.</strong> Dim-1 (59%) separates the 3 OX (tight cluster) from <em>all</em> WT. Notably this is <strong>not</strong> explained by adapter (BA2-adapter WT <code>BA2</code> sits with WT, not OX) <em>or</em> batch (run-2 WT <code>BA1_WT</code> sits with WT, not the run-2 OX) → a real, reproducible OX-vs-WT difference exists at fine scale.</p>
{img("cenh3ox_edgeR_MA_powered_3v3.png","98%")}
<p class="cap">MA + volcano (WT n=3 vs OX n=3). edgeR calls ~5,500 windows at MAPQ≥0 but <strong>rDNA-dominated (4,097)</strong> = multimapper confound; at <strong>MAPQ≥10 it drops to ~845</strong> (arm/peri, rDNA gone). A WT-vs-WT null split gives only 16 — so the OX-vs-WT signal is real but small.</p>
{img("cenh3ox_differential_map.png","98%")}
<p class="cap">The 845 MAPQ≥10 differential windows (542 OX-up red, 303 OX-down blue) across the genome. <strong>They are scattered and show NO compartment preference</strong> (region enrichment ≈1.0 everywhere — <em>not</em> centromere/pericentromere-biased, which is what a CENH3-driven effect should look like). The strongest signals are technical: the Chr5:8.6 Mb OX-down cluster has <strong>normal WGA gDNA (no CNV)</strong> and similar density at MAPQ≥0 — it's a MAPQ-filtering (unique-mapping) artifact, not lost breaks.</p>
<div class="key"><strong>Even taking the difference as genotype (adapter assumed not a confounder):</strong> it does <em>not</em> have the signature of a CENH3 mechanism — there is <strong>no centromere/pericentromere enrichment</strong> (CENH3 is a centromeric protein). So there is no compelling evidence that CENH3 overexpression remodels the physiological DSB landscape.</div>

<h3>5c. Matched gDNA input control — is the difference copy-number/structural?</h3>
<p>A matched <strong>input/gDNA library for the CENH3-OX line</strong> (<code>cenh3ox_input</code>, 22.8 M reads, 99.8% aligned) was mapped to TAIR12 to test whether the differential loci are DNA-content/structural differences.</p>
{img("cenh3ox_input_CNV.png","98%")}
<p class="cap">log2(OX input / Col-0 WGA), 10 kb. Genome-wide flat (the uniform offset is the input's repeat-heavy read distribution; green at centromeres = repeat enrichment) → <strong>no widespread copy-number changes</strong> in the OX line.</p>
{img("cenh3ox_diff_vs_input.png","56%")}
<p class="cap"><strong>The differential does NOT track DNA content.</strong> Across the 845 differential windows, matched-input (MAPQ≥10) coverage is identical for OX-up, OX-down and arm baseline (1.46 / 1.48 / 1.38), and Spearman(logFC, input) = <strong>0.00</strong>. So the OX-vs-WT difference is <em>not</em> a copy-number/structural artifact.</p>
<div class="key"><strong>Net, with the matched input:</strong> the reproducible OX-vs-WT difference is <strong>not explained by copy number/structure</strong> (input is uniform) — but it also has <strong>no CENH3 signature</strong> (not centromeric). Its origin (subtle systematic difference between the OX and WT BLISS libraries vs a genuine non-centromeric genotype effect) cannot be resolved with this design; the gDNA input is a different library type so it doesn't control BLISS-specific biases. A same-batch, same-adapter, replicated WT-vs-OX (ideally same line ±CENH3 induction) is needed to settle it.</div>
<div class="warn"><strong>The honest interpretation.</strong> There <em>is</em> a detectable, reproducible fine-scale difference between OX and WT libraries (MDS Dim-1; ~845 MAPQ≥10 windows). But it is <strong>confounded</strong>: the OX trio uniquely shares <em>run-2 batch + BA2 adapter + OX genotype</em>, and no WT control matches them on batch <em>and</em> adapter. The matched same-batch comparison has only <strong>1 WT replicate</strong> (BA1_WT) so its FDR is unreliable; the replicated 3v3 mixes batches. <strong>So the difference cannot be cleanly attributed to CENH3 overexpression.</strong> The coarse measures (compartment, chromosome landscape, CEN180, genes/TEs) show no large effect; the multimapper-driven part is technical.</div>

<h2 id="concl">6. Conclusion</h2>
<div class="key"><strong>No CENH3-OX effect can be cleanly attributed at this design's resolution.</strong> Coarse readouts (compartment, per-chromosome landscape, CEN180, genes/TEs)
are <strong>indistinguishable</strong> WT vs OX. A fine-scale (1 kb) analysis <em>does</em> separate the OX trio from all WT (MDS Dim-1; ~845 MAPQ≥10 differential windows) — but the design
<strong>confounds genotype with batch + adapter</strong> (OX = run-2 + BA2 adapter; the only same-batch WT, BA1_WT, has the other adapter and no replicates), so this difference cannot be assigned to CENH3.
Much of the raw signal (alignment rate, multimapper %, raw rDNA hotspot counts) is clearly technical. <strong>A definitive test needs WT and OX prepared in the same batch with the same adapter, replicated.</strong></div>
<p class="cap">Pipeline: demux (barcode-verify, UMI→header, 16 bp strip) → trim_galore q20 → bowtie2 --very-sensitive -N1 → umi_tools dedup → 5′-break extraction → 1 kb Poisson hotspots. All on TAIR12 (Col-CC, GCA_028009825.2).</p>
</body></html>"""
open(DST,"w").write(html)
print("Wrote",DST,f"({os.path.getsize(DST)//1024} KB)")
