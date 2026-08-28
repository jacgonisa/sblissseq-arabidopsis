#!/usr/bin/env python3
"""Assemble a self-contained HTML report for the sBLISS-seq analysis."""
import base64, csv, os

BASE = "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
FIG  = f"{BASE}/results/analysis/figures"
OUT  = f"{BASE}/results/analysis/sBLISS_report.html"

def img(path, width="100%"):
    if not os.path.exists(path):
        return f"<p><em>[missing: {os.path.basename(path)}]</em></p>"
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode()
    return f'<img src="data:image/png;base64,{b64}" style="max-width:{width};border:1px solid #ddd;border-radius:4px;">'

# ── pipeline summary ──────────────────────────────────────────────────────────
pipe = [
    # sample, raw R1, accepted, aln%, dedup_in, final, dup%, hotspots (1kb, no blacklist)
    ("BA1", "42,568,398", "36,696,057", "85.0%", "30.0M", "11,239,192", "70%", "2,284"),
    ("BA2", "35,548,814", "31,354,922", "90.2%", "28.0M", "6,509,495",  "80%", "2,276"),
    ("old_BA1_BA2", "52,727,272", "43,370,598", "79.6%", "34.5M", "6,263,466", "90%", "1,952"),
]

# ── region enrichment ─────────────────────────────────────────────────────────
regrows = {}
with open(f"{BASE}/results/dsb_enrichment.csv") as f:
    for r in csv.DictReader(f):
        regrows.setdefault(r["sample"], {})[r["region"]] = r

def reg_cell(s, region):
    d = regrows[s][region]
    return f'{float(d["dsb_per_1k_mapped"]):.0f}'

def reg_mb(s, region):
    d = regrows[s][region]
    return f'{float(d["dsb_sites"]) / float(d["genomic_bp"]) * 1e6:,.0f}'

# ── DSBs per Mb of uniquely-mapped reads (distinct sites / aligned bases) ──────
pmb = {}
with open(f"{BASE}/results/analysis/dsb_per_mapped_mb.csv") as f:
    for r in csv.DictReader(f):
        pmb.setdefault(r["sample"], {})[r["region"]] = r

def pmb_cell(s, region):
    return f'{float(pmb[s][region]["sites_per_mapped_mb"]):,.0f}'

# ── reconciliation table (consistent regions, NOR separated) ──────────────────
rec = {}
with open(f"{BASE}/results/analysis/density_reconciliation.csv") as f:
    for r in csv.DictReader(f):
        rec.setdefault(r["sample"], {})[r["region"]] = r

def rc(s, region, field):
    return f'{float(rec[s][region][field]):,.0f}'

# ── robust hotspots from the region/NOR-annotated table ───────────────────────
# columns: chr start end score region is_NOR top_gene
ann = []
with open(f"{BASE}/results/analysis/hotspots/robust_hotspots.annotated.tsv") as f:
    next(f)  # header
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 7: continue
        ann.append({"chr": c[0], "start": int(c[1]), "end": int(c[2]),
                    "score": int(float(c[3])), "region": c[4],
                    "is_NOR": c[5] == "TRUE", "gene": c[6]})

n_robust = len(ann)
n_nor    = sum(a["is_NOR"] for a in ann)
n_cen    = sum(a["region"] == "centromere" for a in ann)
n_peri   = sum(a["region"] == "pericentromere" for a in ann)
n_arm_real = sum(a["region"] == "arm" and not a["is_NOR"] for a in ann)

# Top real (non-NOR) hotspots — the actual biology, NOR artefacts excluded
real = sorted([a for a in ann if not a["is_NOR"]], key=lambda a: -a["score"])
hot_rows = "\n".join(
    f"<tr><td>{a['chr']}:{a['start']:,}</td><td>{a['region']}</td>"
    f"<td>{a['score']:,}</td><td>{a['gene'] if a['gene'] != '.' else '&ndash;'}</td></tr>"
    for a in real[:15])

# Top NOR artefact hotspots (shown separately so they're transparent, not hidden)
nor_top = sorted([a for a in ann if a["is_NOR"]], key=lambda a: -a["score"])
nor_rows = "\n".join(
    f"<tr><td>{a['chr']}:{a['start']:,}</td><td>{a['score']:,}</td>"
    f"<td>{a['gene'] if a['gene'] != '.' else '&ndash;'}</td></tr>"
    for a in nor_top[:5])

css = """
body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1100px;margin:2em auto;
padding:0 2em;color:#222;line-height:1.55}
h1{border-bottom:3px solid #1f77b4;padding-bottom:.3em}
h2{color:#1f77b4;margin-top:1.8em;border-bottom:1px solid #ddd;padding-bottom:.2em}
h3{color:#555}
table{border-collapse:collapse;margin:1em 0;font-size:.92em}
th,td{border:1px solid #ccc;padding:6px 10px;text-align:right}
th{background:#f0f4f8}td:first-child,th:first-child{text-align:left}
.fig{margin:1.5em 0;text-align:center}
.cap{font-size:.88em;color:#666;margin-top:.4em}
.note{background:#fff8e1;border-left:4px solid #ffb300;padding:.8em 1.2em;margin:1em 0}
.key{background:#e8f5e9;border-left:4px solid #2ca02c;padding:.8em 1.2em;margin:1em 0}
code{background:#f4f4f4;padding:1px 5px;border-radius:3px}
"""

html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>sBLISS-seq Report — Arabidopsis Col-0</title><style>{css}</style></head><body>
<h1>sBLISS-seq DSB Mapping Report</h1>
<p><strong>Organism:</strong> <em>Arabidopsis thaliana</em> Col-0 (ragtag assembly) &nbsp;|&nbsp;
<strong>Samples:</strong> BA1, BA2, old_BA1_BA2 &nbsp;|&nbsp;
<strong>Protocol:</strong> Hidmi et al. 2024 (STAR Protocols), adapted for Arabidopsis &nbsp;|&nbsp;
<strong>Aligner:</strong> HISAT2 (primary alignments only)</p>

<h2>1. Pipeline summary</h2>
<table>
<tr><th>Sample</th><th>Raw R1</th><th>Barcode-passed</th><th>Aln rate</th>
<th>Dedup input</th><th>Unique breaks</th><th>Dup rate</th><th>Hotspots (1kb)</th></tr>
{''.join(f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td><td>{r[4]}</td><td>{r[5]}</td><td>{r[6]}</td><td>{r[7]}</td></tr>" for r in pipe)}
</table>
<p>Pipeline: <code>bliss_demux</code> (UMI + barcode, Hamming&le;1) &rarr; Trim Galore &rarr;
HISAT2 &rarr; <code>umi_tools dedup</code> &rarr; 5&prime; break extraction &rarr; BREAK BED/BEDGRAPH/bigWig.</p>

<h2>2. Genome-wide DSB density landscape</h2>
<p>CPM-normalised DSB density (y = CPM per kb, comparable across bin sizes), <strong>no blacklist</strong>.
Shaded: centromere (red), pericentromere (orange), NOR 45S rDNA (purple); per-chromosome y-scale.
The tall purple spikes at Chr2/Chr4 0&ndash;1 Mb are the collapsed NOR 45S rDNA arrays
(repeat-mapping artefacts); pericentromeres (orange) show genuine DSB enrichment. Shown at three
bin sizes &mdash; coarser (100 kb) reveals the broad compartment trend, finer (1 kb) resolves
individual peaks.</p>
<div class="fig">{img(f"{FIG}/genome_density_100kb.png")}
<div class="cap">100 kb bins &mdash; chromosome-scale compartment view.</div></div>
<div class="fig">{img(f"{FIG}/genome_density_10kb.png")}
<div class="cap">10 kb bins &mdash; intermediate resolution; pericentromeric structure sharper.</div></div>
<div class="fig">{img(f"{FIG}/genome_density_1kb.png")}
<div class="cap">1 kb bins &mdash; peak-level resolution (matches the hotspot tile size). The NOR spikes
become very tall here; per-chromosome y-scaling keeps the arms readable.</div></div>
<div class="note">The Chr2/Chr4 NOR 45S rDNA spikes are so tall they compress the rest of the genome.
The plots below <strong>mask the NOR bins</strong> (Chr2/Chr4 0&ndash;1 Mb, purple) so arm and
pericentromeric structure is readable &mdash; centromeres are still shown.</div>
<div class="fig">{img(f"{FIG}/genome_density_noNOR_100kb.png")}
<div class="cap">NOR-masked, 100 kb bins &mdash; clean chromosome-scale view; pericentromeric
enrichment (orange) now obvious on every chromosome.</div></div>
<div class="fig">{img(f"{FIG}/genome_density_noNOR_10kb.png")}
<div class="cap">NOR-masked, 10 kb bins &mdash; the Chr2 ~4.5 Mb pericentromeric peak (top genuine
hotspot cluster) stands out clearly.</div></div>

<h2>3. DSB density by chromosomal compartment</h2>
<p>DSB density can be normalised three ways, all with the same numerator (distinct break sites):
<strong>(a)</strong> per Mb of DNA (raw spatial), <strong>(b)</strong> per 1,000 mapped reads, and
<strong>(c)</strong> per Mb of uniquely-mapped read sequence (aligned bases). Crucially, all regions
are defined with the <strong>NOR 45S rDNA broken out as its own category</strong> &mdash; never lumped
into "arms".</p>

<div class="note"><strong>Reconciling (b) and (c).</strong> They share the numerator and differ only
in the denominator (read count vs aligned bases = reads&times;read-length). Because mean read length
is nearly constant genome-wide (~125&ndash;132 bp), (c) = (b) &times; (1000/read-length) &mdash; a
near-constant factor of <strong>7.8</strong>. Across all compartments and samples the two are
correlated at <strong>r = 0.998</strong>: they are the <em>same</em> measurement and give the
<em>same</em> ranking.</p></div>

<div class="note" style="background:#ffebee;border-color:#e53935">
<strong>Correction to an earlier conclusion.</strong> An initial version reported
"centromere &gt; pericentromere &gt; arms (939 vs 549)" per 1,000 reads. That was an artefact: the
"arms" region had silently <strong>included the NOR rDNA</strong>, a huge low-complexity pile-up
(4.3M reads but only ~1.0M unique sites = ~226 sites/1k reads) that dragged the arms value down to
549. With the NOR separated, <strong>arms = 949</strong>. The metrics then agree.</div>

<h3>Reconciled compartment density (BA1; BA2 / old_BA1_BA2 concordant &mdash; see figure)</h3>
<table>
<tr><th>Compartment</th><th>(a) sites / Mb DNA</th><th>(b) sites / 1k reads</th><th>(c) sites / mapped-Mb</th></tr>
<tr><td>Arm (NOR-excl.)</td><td>{rc('BA1','arms_noNOR','per_region_Mb')}</td><td>{rc('BA1','arms_noNOR','per_1k_reads')}</td><td>{rc('BA1','arms_noNOR','per_mapped_Mb')}</td></tr>
<tr><td>Pericentromere</td><td>{rc('BA1','pericentromere','per_region_Mb')}</td><td>{rc('BA1','pericentromere','per_1k_reads')}</td><td>{rc('BA1','pericentromere','per_mapped_Mb')}</td></tr>
<tr><td>Centromere core</td><td>{rc('BA1','centromere','per_region_Mb')}</td><td>{rc('BA1','centromere','per_1k_reads')}</td><td>{rc('BA1','centromere','per_mapped_Mb')}</td></tr>
<tr><td>rDNA (NOR)</td><td>{rc('BA1','nor','per_region_Mb')}</td><td>{rc('BA1','nor','per_1k_reads')}</td><td>{rc('BA1','nor','per_mapped_Mb')}</td></tr>
</table>
<div class="fig">{img(f"{FIG}/dsb_per_mapped_mb.png")}
<div class="cap">Metric (c) on its own: distinct DSB sites per Mb of uniquely-mapped reads. Genuine
compartments cluster ~6,900&ndash;7,700; the NOR rDNA collapses to ~1,700&ndash;2,250 (redundant pile-up).</div></div>
<div class="fig">{img(f"{FIG}/density_reconciliation.png")}
<div class="cap">All three normalisations, all three samples, consistent NOR-separated regions.
Rows (b) and (c) are proportional (identical shape). The raw spatial view (a) looks different
<em>only</em> because heterochromatin recruits more sequencing depth per Mb (mappability), and
because rDNA has many distinct positions packed into 2 Mb.</div></div>

<div class="key"><strong>Reconciled result:</strong> per unit of successfully-mapped sequence, DSB
density is <strong>Arm &asymp; Centromere &gt; Pericentromere</strong> &mdash; fairly uniform across
genuine chromatin (within ~10%), with the centromere as break-prone as the arms (it is <em>not</em>
depleted, nor dramatically enriched). The <strong>NOR 45S rDNA is the only strong outlier</strong>
(~4&times; lower unique-site density = redundant repeat pile-up). The raw per-Mb view exaggerates
heterochromatin purely through differential mappability.</div>

<h3>3d. Centromere fine-structure: sBLISS vs CENH3 ChIP-seq</h3>
<p>Zooming into each centromere (core &plusmn;1.5 Mb) and overlaying sBLISS DSB density with CENH3
ChIP-seq (Col-0 ragtag, 10 kb).</p>
<div class="fig">{img(f"{FIG}/cenh3_zoom.png")}
<div class="cap">Top: sBLISS DSB density (3 samples). Bottom: CENH3 ChIP. Purple = centromere core.
CENH3 peaks within the core; sBLISS DSBs peak in the <em>flanking pericentromere</em> and dip over
the CENH3 core (Pearson r = &minus;0.21, Spearman &minus;0.34 within the zoom).</div></div>
<div class="note"><strong>Interpretation &amp; caveat.</strong> DSBs are anti-correlated with the
CENH3-defined functional core and concentrate in the flanking pericentromeric heterochromatin.
However, the CENH3 satellite core (CEN178 arrays) is the least mappable part of the genome for short
sBLISS reads, so part of the core dip is a mappability floor rather than true protection &mdash;
interpret the core values cautiously.</div>

<h2>4. DSB density at genes (TSS / gene body / TES)</h2>
<div class="fig">{img(f"{FIG}/profile_TSS.png","48%")} {img(f"{FIG}/profile_TES.png","48%")}</div>
<div class="fig">{img(f"{FIG}/profile_genebody.png","70%")}
<div class="cap">DSBs are mildly depleted within gene bodies relative to flanks, with a
characteristic peak just downstream of the TSS. BA1/BA2 are near-identical; old_BA1_BA2
runs higher (lower-complexity library, more residual signal).</div></div>
<div class="fig">{img(f"{FIG}/heatmap_genebody.png","60%")}
<div class="cap">Per-gene heatmap of DSB signal, TSS&rarr;TES scaled.</div></div>

<h2>5. DSB density at transposons</h2>
<div class="fig">{img(f"{FIG}/profile_TE.png","60%")}
<div class="cap">DSB signal over all annotated TEs (EDTA repeat_regions).</div></div>
<div class="fig">{img(f"{FIG}/profile_TE_superfamily.png")}
<div class="cap">BA1 DSB density resolved by TE superfamily. Copia LTRs show a sharp
3&prime;-end break peak; Mutator TIRs peak internally.</div></div>

<h2>6. DSB hotspots (BRGenomics)</h2>
<p>Genome tiled into 1 kb windows; <strong>raw</strong> break counts per tile; hotspots called by
Poisson test (BH-adjusted <em>p</em>&lt;0.001, &ge;3&times; genome-mean).
<strong>No blacklist is applied</strong> &mdash; centromeres and the NOR 45S rDNA arrays are
analysed like any other region. Instead, every hotspot is <strong>annotated</strong> with its
chromosomal compartment (arm / pericentromere / centromere) and a NOR flag, so repetitive-pileup
artefacts stay identifiable rather than being silently removed.</p>
<div class="note"><strong>The strongest raw hotspots are NOR 45S rDNA artefacts.</strong>
{n_nor:,} of the {n_robust:,} robust hotspots fall in the collapsed rDNA arrays at Chr2/Chr4
0&ndash;1 Mb (<code>is_NOR=TRUE</code>). These are repeat-mapping pile-ups, not physiological DSB
hotspots, and should be excluded when interpreting biology &mdash; they are listed separately below
for transparency.</div>
<ul>
<li><strong>{n_robust:,} robust hotspots</strong> reproducible across all three samples.</li>
<li><strong>{n_nor:,}</strong> are NOR 45S rDNA artefacts; the remaining <strong>{n_robust - n_nor:,}</strong>
are genuine.</li>
<li>By compartment: <strong>{n_peri:,}</strong> pericentromeric, <strong>{n_cen:,}</strong>
<em>centromeric</em> (now recovered &mdash; previously blacklisted), <strong>{n_arm_real:,}</strong>
on chromosome arms.</li>
</ul>
<h3>Top genuine (non-NOR) DSB hotspots</h3>
<table>
<tr><th>Location</th><th>Compartment</th><th>Breaks (BA1)</th><th>Gene</th></tr>
{hot_rows}
</table>
<p class="cap">The strongest genuine hotspots cluster in Chr2 pericentromeric heterochromatin
(~4.5 Mb). The top <em>gene</em>-level arm hotspot, AT5G10250 (Chr5:3.26 Mb, 2,137 breaks), is
analogous to the physiological break hotspots (e.g. MIR21) reported in the original sBLISS study.
Centromere cores now contribute {n_cen:,} robust hotspots (e.g. Chr4:7.03 Mb, 593 breaks) &mdash;
consistent with the reconciled compartment analysis (Section 3), where the centromere is as
break-prone as the arms per mappable read (not depleted).</p>
<h3>Top NOR 45S rDNA artefact hotspots (excluded from biology)</h3>
<table>
<tr><th>Location</th><th>Breaks (BA1)</th><th>Overlapping rDNA gene model</th></tr>
{nor_rows}
</table>
<div class="note"><strong>The rDNA paradox.</strong> The NOR is <em>enriched</em> in raw counts
(it dominates the table above) yet <em>depleted</em> per mapped read (Section 3: 226 vs 949
sites/1k reads). Both are artefacts, not contradictions. The 45S rDNA exists in hundreds&ndash;
thousands of copies but the assembly collapses them into a ~1 Mb stub: reads from every copy pile
onto it (&rarr; huge raw counts), while the tiny stub <strong>saturates</strong> &mdash; nearly every
position is already a break site and the 8 bp UMI saturates, so <code>umi_tools</code> merges genuinely
distinct breaks (4.4 reads/unique-site here vs 1.05 on the arms). New unique sites stop accumulating,
so per-read density looks low. <strong>rDNA fragility is real and consistent with the raw signal, but
it cannot be quantified with this assay+assembly</strong> (would need an rDNA-resolved reference and a
longer UMI) &mdash; which is exactly why the NOR is flagged as its own category throughout.</div>

<h3>Cross-sample reproducibility (Venn)</h3>
<div class="fig">{img(f"{FIG}/hotspots_venn.png")}
<div class="cap">3-way overlap of 1 kb hotspot tiles. <strong>Left:</strong> all hotspots &mdash; 1,886 shared
by all three samples. <strong>Right:</strong> rDNA (NOR) excluded &mdash; 581 genuine hotspots are
reproducible across all three. The bulk of the shared signal is the rDNA pile-up; the genuine
3-way core (581) is the high-confidence physiological set.</div></div>

<h3>Hotspot classification by compartment</h3>
<div class="fig">{img(f"{FIG}/hotspots_classification.png")}
<div class="cap">Every hotspot classified as Arm / Pericentromere / Centromere, with NOR 45S rDNA as a
separate artefact category. rDNA dominates the raw counts in every sample (~1,300); among genuine
hotspots, pericentromere is the largest class, centromeres are now recovered, and arm hotspots are
the rarest. old_BA1_BA2 has fewer genuine hotspots (low complexity).</div></div>

<h2>7. Sample reproducibility (correlation)</h2>
<p>Genome binned into 10 kb windows; per-bin DSB density (CPM) compared between every pair of
samples. Pearson <em>r</em> is computed on log<sub>10</sub>(CPM+1) with the NOR 45S rDNA bins
<strong>excluded</strong> (shown red) &mdash; the rDNA pile-up is a reproducible artefact that would
otherwise inflate every correlation toward ~1 and hide the real differences between libraries.</p>
<table>
<tr><th>Pair</th><th>Pearson (log<sub>10</sub>)</th><th>Spearman</th></tr>
<tr><td>BA1 vs BA2</td><td>0.896</td><td>0.859</td></tr>
<tr><td>BA1 vs old_BA1_BA2</td><td>0.797</td><td>0.693</td></tr>
<tr><td>BA2 vs old_BA1_BA2</td><td>0.819</td><td>0.713</td></tr>
</table>
<div class="fig">{img(f"{FIG}/sample_correlation.png")}
<div class="cap">Pairwise DSB-density scatter (10 kb bins, log&ndash;log). Dashed line = identity.
Red points are NOR 45S rDNA bins (excluded from <em>r</em>).</div></div>
<div class="key"><strong>Result:</strong> the two genuine libraries <strong>BA1 and BA2 are strongly
concordant (r = 0.90)</strong>. <strong>old_BA1_BA2 correlates noticeably worse</strong> (0.69&ndash;0.82,
most clearly in the rank-based Spearman) &mdash; consistent with its ~90% duplication / low complexity,
so it should be treated as a qualitative replicate only. The NOR bins (red) lie tightly on the
identity line in every panel, confirming the rDNA pile-up is a reproducible mapping artefact rather
than biology.</div>

<h2>8. Methods &amp; reproducibility notes</h2>
<ul>
<li><strong>Reference:</strong> Col-0 ragtag assembly; annotation = Liftoff (TAIR10 genes) + EDTA (TEs).</li>
<li><strong>Region definitions:</strong> centromeres from PoreC project BED; pericentromere = &plusmn;2 Mb flanks.
Centromere analysis is region-based (not gene/TE GFF) due to repetitive-mapping caveats.</li>
<li><strong>No blacklist:</strong> hotspot calling and the density landscape use <strong>no masking</strong>
&mdash; centromeres and the NOR2/NOR4 45S rDNA arrays (Chr2/Chr4 0&ndash;1 Mb) are analysed and instead
<strong>flagged</strong> (region label + <code>is_NOR</code>) so artefacts remain identifiable.
See <code>hotspots/robust_hotspots.annotated.tsv</code>.</li>
<li><strong>Normalisation:</strong> CPM (per total breaks) for cross-sample tracks &amp; metaplots;
per-mapped-read for compartment density.</li>
<li><strong>BA1 vs BA2</strong> are highly concordant; old_BA1_BA2 has 90% duplication
(low complexity) and should be treated as a qualitative replicate only.</li>
</ul>
</body></html>
"""

with open(OUT, "w") as f:
    f.write(html)
print("Wrote", OUT, f"({os.path.getsize(OUT)//1024} KB)")
