#!/usr/bin/env python3
"""Clean rebuild of the sBLISS TAIR12 (Col-CC) report — 11 sections, narrative-ordered.
Centerpiece: how DSB density normalisation + MAPQ handling change every conclusion."""
import base64, os
BASE = "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026"
FIG  = f"{BASE}/results_TAIR12/analysis/figures"
OUT  = f"{BASE}/results_TAIR12/analysis/sBLISS_TAIR12_report.html"

def img(path, width="100%"):
    if not os.path.exists(path): return f"<p><em>[missing: {os.path.basename(path)}]</em></p>"
    with open(path,"rb") as f: b=base64.b64encode(f.read()).decode()
    return f'<img src="data:image/png;base64,{b}" style="max-width:{width};border:1px solid #ddd;border-radius:4px;">'
def f(p,w="100%"): return img(f"{FIG}/{p}",w)

css="""body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1080px;margin:1.5em auto;padding:0 2em;color:#222;line-height:1.5}
h1{border-bottom:3px solid #d62728;padding-bottom:.3em}h2{color:#d62728;margin-top:1.6em;border-bottom:1px solid #ddd;padding-bottom:.2em}
h3{color:#555;margin-top:1.2em}table{border-collapse:collapse;margin:1em 0;font-size:.9em}th,td{border:1px solid #ccc;padding:5px 9px;text-align:right}
th{background:#fdecea}td:first-child,th:first-child{text-align:left}.fig{margin:1.2em 0;text-align:center}.cap{font-size:.86em;color:#666;margin-top:.3em}
.note{background:#fff8e1;border-left:4px solid #ffb300;padding:.7em 1.1em;margin:1em 0}.key{background:#e8f5e9;border-left:4px solid #2ca02c;padding:.7em 1.1em;margin:1em 0}
.warn{background:#ffebee;border-left:4px solid #e53935;padding:.7em 1.1em;margin:1em 0}code{background:#f4f4f4;padding:1px 5px;border-radius:3px}
nav{background:#fafafa;border:1px solid #eee;padding:.6em 1.2em;font-size:.9em}nav a{margin-right:1em;color:#1f6feb;text-decoration:none}</style>"""

html=f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>sBLISS TAIR12 (Col-CC) Report</title><style>{css}</head><body>
<h1>sBLISS-seq DSB mapping — <em>Arabidopsis</em> Col-CC (TAIR12, T2T)</h1>
<p><strong>Samples:</strong> BA1, BA2, old_BA1_BA2 &nbsp;|&nbsp; <strong>Assembly:</strong> GCA_028009825.2 Col-CC (rDNA-resolved)
&nbsp;|&nbsp; <strong>Aligners:</strong> bowtie2 <code>--very-sensitive</code> + minimap2 <code>-ax sr</code> &nbsp;|&nbsp;
<strong>gDNA control:</strong> WGA (E-MTAB-6257). Reads PE150, analysed single-end on R1 (sBLISS).</p>
<nav><b>Sections:</b>
<a href="#s1">1 Pipeline</a><a href="#s2">2 Normalisation</a><a href="#s3">3 Landscape</a><a href="#s4">4 Enrichment</a>
<a href="#s5">5 rDNA</a><a href="#s6">6 Mappability</a><a href="#s7">7 Genes/TEs</a><a href="#s8">8 Reproducibility</a>
<a href="#s9">9 Satellite phasing</a><a href="#s10">10 BiCroLab</a><a href="#s11">11 Hotspot context</a><a href="#s12">12 Methods</a></nav>

<div class="key"><strong>One-line summary.</strong> sBLISS DSBs are distributed fairly uniformly per mapped read across
arms, pericentromere and centromere; the apparent compartment "enrichment/depletion" is mostly a
mappability/coverage effect. The one genuine, gDNA-confirmed hotspot is the <strong>45S rDNA (~12–18×)</strong>.
Every such statement depends on two analysis choices — MAPQ handling and normalisation — which §2 and §4 make explicit.</div>

<h2 id="s1">1. Pipeline, alignment & complexity</h2>
<p>bliss_demux (UMI + barcode) &rarr; Trim Galore &rarr; align (bowtie2 &amp; minimap2) &rarr; umi_tools dedup &rarr;
5&prime; break extraction. Mean trimmed R1 ≈ 127 bp.</p>
<table>
<tr><th>Sample</th><th>bowtie2 mapped</th><th>minimap2 mapped</th><th>dedup (bt2/mm2)</th><th>dup rate</th></tr>
<tr><td>BA1</td><td>30,444,499</td><td>30,708,472</td><td>11.4M / 13.1M</td><td>~60%</td></tr>
<tr><td>BA2</td><td>28,445,689</td><td>28,915,279</td><td>6.6M / 7.7M</td><td>~75%</td></tr>
<tr><td>old_BA1_BA2</td><td>34,852,498</td><td>35,026,503</td><td>6.4M / 8.5M</td><td>~82%</td></tr>
</table>
<p class="cap">The two aligners are highly concordant (mapped counts within ~1%); minimap2 is much faster.
TAIR12 alignment rates are slightly above ragtag (recovered rDNA reads).</p>
<h3>1a. Library complexity (UMIs per break location)</h3>
<div class="fig">{f("bicrolab_umi_threshold.png","62%")}
<div class="cap">Distinct break sites retained at increasing UMI thresholds. BA1 is the most complex; old_BA1_BA2 the least.</div></div>

<h2 id="s2">2. How DSB density is normalised <em>(read this first)</em></h2>
<div class="note">Every "enrichment" number in this report is set by <strong>two independent choices</strong>:
<ul style="margin:.3em 0">
<li><strong>Numerator:</strong> break <em>events</em> (&approx; mapped reads; inflated by repeat pile-ups) vs
<em>distinct sites</em> (deduplicated positions; robust in repeats).</li>
<li><strong>Denominator:</strong>
 <b>(1) per library</b> (CPM) &mdash; corrects depth <em>between samples</em>, no effect on within-genome ranking;
 <b>(2) per region length</b> (per Mb of DNA) &mdash; spatial <em>amount</em>, confounded by mappability;
 <b>(3) per mapped read / aligned base</b> &mdash; intrinsic <em>propensity</em>, corrects mappability;
 <b>(4) &divide; gDNA</b> &mdash; corrects mappability <em>and</em> DNA-content/copy-number (the unbiased one).</li>
</ul>
And a third axis &mdash; <strong>MAPQ handling</strong>: keep multimappers (MAPQ&ge;0, each at its best position)
vs drop them (MAPQ&ge;1). Repeats (centromere, rDNA) are ~2&ndash;4% uniquely mappable, so this flips them between
"enriched" and "invisible". <strong>&sect;4 shows all of this on one figure.</strong></div>

<h3>2a. The exact formulas used</h3>
<p>Notation: for genomic bin/region <em>i</em>, let <code>b<sub>i</sub></code> = sBLISS break events (5&prime; ends of
deduplicated reads), <code>g<sub>i</sub></code> = WGA gDNA coverage, <code>L<sub>i</sub></code> = length (bp).
Totals over the library: <code>B = &Sigma;b<sub>i</sub></code>, <code>G = &Sigma;g<sub>i</sub></code>.</p>
<pre style="background:#f6f8fa;border:1px solid #ddd;border-radius:6px;padding:12px 14px;font-size:13px;line-height:1.7;overflow-x:auto">
(1) Per library — CPM (depth, between-sample):
      CPM_i  =  b_i / B  &times; 10^6
      density shown as  CPM per kb  =  CPM_i / (L_i / 1000)

(2) Per region length — amount of DNA (per Mb):
      density_i  =  b_i / (L_i / 10^6)          [breaks per Mb]

(3) Per mapped read / aligned base — intrinsic propensity (mappability-corrected):
      sites_per_1k_i  =  (distinct break sites in i) / (mapped reads in i) &times; 1000
      (equivalently breaks per Mb of aligned bases, via samtools bedcov)

(4) &divide; gDNA — DNA-content / copy-number corrected (the unbiased enrichment):
      enrichment_i  =  (b_i / B)  /  (g_i / G)
                    =  fraction of DSBs in i  /  fraction of gDNA in i
      &gt; 1  =  enriched beyond DNA content & mappability ;  &lt; 1  =  depleted

(5) Hotspot call vs gDNA baseline (per-window Poisson, &sect;8):
      expected  &lambda;_i  =  g_i &times; (B / G)        [local gDNA scaled to BLISS depth]
      call if   Poisson padj(b_i ; &lambda;_i) &lt; 0.001   AND   fold = b_i / &lambda;_i &ge; 3   AND   b_i &ge; 10
      (flat-Poisson variant instead uses a constant &lambda; = B / n_windows)

(6) gDNA-controlled IGV track — matched-unit log2 ratio:
      log2ratio_i  =  log2( (CPM_BLISS,i + c) / (CPM_gDNA,i + c) ) ,   c = 1 (pseudocount)
      both tracks as READ COVERAGE in 500 bp bins, each CPM-normalised
      &gt; 0  =  DSB-enriched over gDNA ;  ~&minus;0.6  =  typical arm ;  strongly &minus;  =  coverage artifact

(7) Metaplot shape normalisation (signal/mean, &sect;7b):
      shape(x)  =  profile(x) / mean_x[ profile ]      per track, then overlay BLISS vs gDNA
      (removes the magnitude/unit gap so DSB-specific structure is visible against flat gDNA)
</pre>
<p class="cap"><strong>Why (4) and (6) are the honest ones.</strong> (1)&ndash;(3) describe the BLISS signal alone and are all
confounded to some degree by mappability or copy number; dividing by the WGA gDNA control (4 at the compartment level,
6 per-bin for IGV) removes both, so what remains is genuine DSB propensity. Form (6) was rebuilt in matched units
(coverage-vs-coverage, CPM, pseudocount) because the naive log2 of sparse break-points over dense coverage is negative
everywhere and only its <em>shape</em> is meaningful.</p>

<h2 id="s3">3. Genome-wide DSB landscape</h2>
<div class="fig">{f("genome_density_100kb.png")}
<div class="cap">CPM/kb, 100 kb windows, no blacklist. The 45S rDNA (purple, Chr2/Chr4 0&ndash;5 Mb) is a broad plateau across the
resolved array; pericentromeres (orange) are elevated. (y = CPM per kb; events-based, so repeats are inflated.)</div></div>
<div class="fig">{f("genome_density_noNOR_100kb.png")}
<div class="cap">NOR masked, <strong>fixed y-axis 0&ndash;20 CPM/kb</strong> (comparable across chromosomes). Pericentromeric
peaks clip at 20; arms sit ~2&ndash;5. Pericentromeres are clearly the elevated non-rDNA compartment.</div></div>
<h3>3-i. Finer window resolution (50 kb and 10 kb)</h3>
<p>Same landscape at higher resolution. Coarser 100 kb windows give the compartment overview; 50 kb and 10 kb
windows resolve individual pericentromeric clusters and the structure of the rDNA arrays.</p>
<div class="fig">{f("genome_density_noNOR_50kb.png")}
<div class="cap"><strong>50 kb windows</strong>, NOR masked, fixed y 0&ndash;20 CPM/kb.</div></div>
<div class="fig">{f("genome_density_noNOR_10kb.png")}
<div class="cap"><strong>10 kb windows</strong>, NOR masked, fixed y 0&ndash;20 CPM/kb &mdash; sharp pericentromeric peaks
(e.g. the Chr2 ~8.8&ndash;9.3 Mb cluster) become individually resolved.</div></div>
<h3>3a. Effect of multimapper handling (MAPQ≥0 vs MAPQ≥1)</h3>
<div class="fig">{f("genome_density_mm2_processing.png","85%")}
<div class="cap">BA1, minimap2. raw &amp; dedup (keep MAPQ0) carry signal at the centromere cores and NOR; dedup+MAPQ&ge;1
(red) collapses to ~0 there &mdash; multimappers dropped &rArr; repeats become deserts. Arms unaffected.</div></div>
<h3>3b. sBLISS vs WGA gDNA, genome-wide</h3>
<div class="fig">{f("genome_bliss_vs_wga.png","85%")}
<div class="cap">WGA gDNA (black) is ~uniform &mdash; including across the NOR; sBLISS (red) spikes ~100 CPM/kb at the
Chr2/Chr4 NOR. So the rDNA signal is <strong>DSB-specific, not DNA-content</strong>.</div></div>

<h2 id="s4">4. Compartment enrichment — the master comparison</h2>
<p>The central result. Same BA1 data, four ways (MAPQ handling × normalisation), fold vs genome expectation:</p>
<div class="fig">{f("master_comparison.png")}
<div class="cap">Follow the <strong>45S rDNA row</strong>: 0.14 (drop multimappers) &rarr; 8.24 (keep) &rarr; <strong>12.3&times;</strong>
(&divide;gDNA). The "answer" depends entirely on the two choices; <strong>&divide;gDNA (rightmost) is the unbiased one</strong>.</div></div>
<h3>4a. Per region length vs per mapped read (MAPQ≥1, uniquely-mapped)</h3>
<div class="fig">{f("compartment_enrichment_mm2q1_norm.png")}
<div class="cap"><strong>A</strong> per Mb of DNA: pericentromere&gt;arms, cen/NOR low (mappability). <strong>B</strong> per mapped read:
~uniform (centromere&asymp;arms&asymp;peri, ~940&ndash;985 sites/1k reads) &mdash; intrinsic break propensity is flat across genuine chromatin.</div></div>
<h3>4b. ÷ WGA gDNA (DNA-content corrected — the unbiased enrichment)</h3>
<div class="fig">{f("bliss_vs_gdna_minimap2.png")}
<div class="cap">A (no gDNA) vs B (÷gDNA). Only the rDNA survives as genuinely enriched: <strong>45S ~12–18×, 5S ~1.6–1.9×</strong>;
arms/peri/cen sit &lt;1 because the rDNA captures a large share of all DSB reads.
<strong>bowtie2 cross-check (matched aligner):</strong> NOR 12.25&times;, 5S 1.53&times;, peri 0.70, arms 0.46, cen 0.26 &mdash; identical to minimap2.</div></div>
<div class="key"><strong>Conclusion.</strong> (i) Per mapped read, DSB propensity is <strong>uniform across arms, pericentromere
and centromere</strong> — heterochromatin is not a cold spot; its low per-Mb values are coverage/mappability.
(ii) The one genuine, gDNA-confirmed DSB hotspot is the <strong>45S rDNA (~12–18×)</strong>, with 5S rDNA mildly enriched.
(iii) The "pericentromere &gt; arms" seen per-Mb is a real but modest amount difference driven by coverage, not propensity.</div>

<h2 id="s5">5. The 45S rDNA — resolved assembly + gDNA confirmation</h2>
<p>The ragtag assembly collapsed the 45S rDNA into a ~2 Mb stub; Col-CC resolves it as 9.26 Mb of array.</p>
<div class="fig">{f("rdna_resolution_compare.png")}
<div class="cap">Resolving the assembly recovers ~half the per-read signal (226&rarr;427 sites/1k reads) and halves read
redundancy (4.4&rarr;2.3/site) vs the collapsed ragtag stub.</div></div>
<div class="key"><strong>Is the rDNA truly break-enriched, or just copy number?</strong> The WGA gDNA control answers it:
gDNA sits at ~0.67× across the (resolved) NOR — random DNA does <em>not</em> pile there — while sBLISS is 8× by events,
giving <strong>BLISS/gDNA ≈ 12–18×</strong>. So the rDNA is a <strong>genuine DSB hotspot</strong> (consistent with rDNA being a
recombination/replication-stress fragile locus), not a copy-number/mappability artefact. (On the <em>collapsed</em> ragtag
stub gDNA would also pile up, hiding this — which is why the resolved assembly + gDNA control were both needed.)</div>

<h2 id="s6">6. Mappability (why the repeats are a blind spot)</h2>
<div class="fig">{f("mappability_compartments.png","78%")}
<div class="cap">GenMap k=127 (=read length): uniquely mappable = arms 97.5%, pericentromere 80.5%, <strong>centromere 3.9%,
45S NOR 2.0%</strong>. Panel A: rDNA is the most break-dense per Mb of DNA (~8× arms, all reads). Panel C: within the
~2&ndash;4% observable satellite, density is arm-like.</div></div>
<div class="fig">{f("genome_coverage_q10_100kb.png","85%")}
<div class="cap">MAPQ&ge;10 read coverage: the centromere cores and NOR are coverage deserts &mdash; the unique-mapping blind spot.</div></div>

<h2 id="s7">7. DSB density at genes &amp; transposons</h2>
<div class="fig">{f("profile_TSS.png","48%")} {f("profile_TES.png","48%")}</div>
<div class="fig">{f("profile_genebody.png","68%")}<div class="cap">Scaled gene body (TAIR12 annotation, 26,867 protein-coding genes), ±2 kb flanks.</div></div>
<div class="fig">{f("heatmap_genebody.png","55%")}</div>
<div class="fig">{f("profile_TE.png","55%")} {f("profile_TE_superfamily.png","55%")}<div class="cap">TEs (EDTA), all + by superfamily.</div></div>
<h3>7a. The gDNA baseline: naked DNA is not flat at the TSS</h3>
<p>The raw profiles above could reflect DSB propensity <em>or</em> merely how accessible/mappable naked DNA is at
these features. The WGA gDNA control separates the two. Crucially, gDNA is <em>not</em> perfectly flat &mdash; random
naked DNA already peaks ~5% at the TSS from accessibility/GC/mappability &mdash; so the DSB profile must be judged
against <em>this</em> baseline, not against zero.</p>
<div class="fig">{f("profile_TSS_gdna_only.png","48%")}
<div class="cap"><strong>WGA gDNA coverage alone</strong> peaks ~5% at the TSS. This mild baseline is what the sBLISS
profile is compared against in &sect;7b.</div></div>

<h3>7b. Reconciled: signal/mean overlay <em>and</em> gDNA-normalised log2 ratio are the same signal</h3>
<p>Earlier drafts split this into a &ldquo;log2(BLISS/gDNA) looks modest/noisy&rdquo; panel and a &ldquo;signal/mean
overlay looks clearly enriched&rdquo; panel &mdash; which read as a contradiction. They are <strong>not</strong>: both use
the identical two profiles and differ only in how they are combined, and algebraically
<code>log2(BLISS/gDNA) = log2(shape<sub>BLISS</sub>/shape<sub>gDNA</sub>) + const</code>. The old log2 <em>track</em> looked
noisy only because it divided a <em>sparse 5&prime;-break</em> CPM by <em>dense fragment-coverage</em> CPM with a &plus;1
pseudocount &mdash; a unit mismatch, not weak biology. Computing the ratio on the <strong>mean-scaled shapes</strong>
removes that artifact and yields one honest curve.</p>
<div class="fig">{f("metaplot_gene_te_reconciled.png","99%")}
<div class="cap"><strong>Top row (signal/mean):</strong> each track scaled to its own mean; sBLISS (red) departs from the
near-flat gDNA baseline (grey). <strong>Bottom row (log2 sBLISS/gDNA on those shapes):</strong> the same departure,
quantified = DSB enrichment <em>above the naked-DNA expectation</em>. <b>TSS</b> +0.43 (<strong>1.35&times;</strong>) just
downstream, with an upstream dip; <b>gene body</b> +0.28 (<strong>1.22&times;</strong>) at the 5&prime; end declining to
3&prime; depletion; <b>transposons</b> a sharp 5&prime;-end spike +0.40 (<strong>1.33&times;</strong>) over a noisier
(short-feature) body. Row&nbsp;2 is literally the vertical gap between the row-1 curves in log space &mdash; same numbers,
two views.</div></div>
<div class="key"><strong>Conclusion.</strong> Beyond the coverage/copy-number signal, physiological DSBs are
<strong>genuinely but modestly enriched (~1.2&ndash;1.4&times;)</strong> at TSS-proximal and 5&prime; feature ends, above the
naked-DNA (gDNA) baseline &mdash; a real, position-specific pattern. The previous &ldquo;modest/noisy&rdquo; log2 track was a
plotting artifact of sparse-break-over-dense-coverage units, now resolved by the shape-normalised ratio.</div>
<h3>7c. GO enrichment of high-DSB protein-coding genes &mdash; nothing</h3>
<div class="note">We ranked the 26,715 protein-coding genes (&ge;5 breaks, &ge;20&times; gDNA) by DSB content and tested the
top 5% (1,336 genes) for GO enrichment with <code>clusterProfiler</code> + <code>org.At.tair.db</code>, using <em>two</em>
ranking metrics (DSB-specific BLISS/gDNA enrichment, and raw DSB density per kb) &times; all three ontologies (BP, MF, CC),
against the protein-coding background. <strong>Zero significantly enriched terms in all six tests.</strong> Physiological
DSBs at protein-coding genes show <em>no</em> functional-category preference &mdash; consistent with breaks being governed by
chromatin context and genomic location (rDNA, pericentromeric heterochromatin) rather than by gene function.</div>

<h3>7d. Are transposons enriched for DSBs relative to genes? &mdash; a formal test</h3>
<p>The occupancy/per-Mb view above suggests DSBs favour TEs, but that mixes in the compartment the TEs live in. Here we
test it <em>per feature</em>: DSB density (break events / kb) &divide; WGA gDNA coverage over the same feature, <strong>gene
vs TE</strong>, with a Mann&ndash;Whitney U test (effect size = Cliff&rsquo;s &delta;; bootstrap CI on the median ratio),
stratified by compartment (rDNA excluded).</p>
<div class="fig">{f("te_vs_gene_enrichment.png","99%")}
<div class="cap"><strong>Per-feature gDNA-normalised DSB density, genes vs TEs, by compartment</strong> (rows: genome-wide,
arms only, pericentromere only; columns: the three libraries). Titles give TE/gene median ratio, Cliff&rsquo;s &delta; and
significance.</div></div>
<table>
<tr><th>Within-compartment (TE/gene, &divide;gDNA)</th><th>BA1</th><th>BA2</th><th>reading</th></tr>
<tr><td>Genome-wide (rDNA excl.)</td><td>0.99 (&delta;&minus;0.03)</td><td>0.98 (&delta;&minus;0.04)</td><td>TE &asymp; gene &mdash; no enrichment</td></tr>
<tr><td>Chromosome arms only</td><td><strong>0.69</strong> (&delta;&minus;0.29)</td><td><strong>0.76</strong> (&delta;&minus;0.22)</td><td>TEs <strong>depleted</strong> vs genes (least-fragile features)</td></tr>
<tr><td>Pericentromere only</td><td>0.99 (&delta;&minus;0.03)</td><td>0.94 (&delta;&minus;0.08)</td><td>TE &asymp; gene &mdash; both elevated</td></tr>
</table>
<table>
<tr><th>Across compartments (pericentromere &divide; arm, same feature)</th><th>BA1</th><th>BA2</th><th>old</th></tr>
<tr><td><strong>TEs</strong>, peri/arm density</td><td><strong>1.61&times;</strong></td><td>1.42&times;</td><td>1.43&times;</td></tr>
<tr><td>genes, peri/arm density</td><td>1.12&times;</td><td>1.15&times;</td><td>1.15&times;</td></tr>
</table>
<p class="cap">TEs outnumber genes <strong>2.4 : 1</strong> in the pericentromere (8,657 vs 3,622) &mdash; the reverse of genome-wide.</p>
<div class="key"><strong>Conclusion.</strong> The &ldquo;DSBs are enriched in TEs&rdquo; statement is largely a <strong>compartment
effect</strong>, not an intrinsic property of transposons: per feature a TE is no more break-prone than a gene genome-wide,
and <strong>on the arms TEs are the <em>least</em>-fragile features</strong> (0.69&times;). The idea that <strong>pericentromeres
are fragile &ldquo;because of&rdquo; TEs is partly borne out</strong> &mdash; the pericentromeric DSB excess is carried
disproportionately by TEs (they gain ~1.4&ndash;1.6&times; break density from the compartment vs ~1.13&times; for genes, and
dominate the sequence 2.4:1) &mdash; but the driver is the <strong>heterochromatic pericentromeric context</strong> acting on
that TE-rich DNA, not a per-element TE fragility (within the pericentromere TEs are no more broken than the genes beside them).
<em>Script:</em> <code>te_vs_gene_enrichment.py</code>. Caveat: pericentromeric mappability is lower, so peri magnitudes rest
on the uniquely-mappable subset; direction is robust across BA1/BA2.</div>
<p>Because the TE/gene ratio is a <em>within-sample</em> relative measure, it is largely robust to the adapter/batch
confound that limits the raw CENH3-OX-vs-WT landscape comparison (&sect;12) &mdash; so we also test it on the
CENH3-overexpression libraries.</p>
<div class="fig">{f("te_vs_gene_by_genotype.png","82%")}
<div class="cap"><strong>CENH3 overexpression does not change transposon-vs-gene DSB enrichment.</strong> TE/gene ratio
(&divide;gDNA, bootstrap-CI) per library, WT (blue) vs CENH3-OX (red), in each compartment. The OX libraries fall within
the WT spread throughout &mdash; and the batch-matched WT (<code>BA1_WT</code>, same delivery) sits among the OX points in
every compartment. (The high outlier is <code>old_BA1_BA2</code>, the low-complexity library.) Consistent with the overall
result that CENH3-OX leaves the DSB landscape unchanged (&sect;12).</div></div>

<h2 id="s8">8. Reproducibility &amp; hotspots</h2>
<div class="fig">{f("sample_correlation.png","80%")}
<div class="cap">10 kb-bin correlation (NOR excluded). BA1&harr;BA2 strongly concordant; old_BA1_BA2 weaker (low complexity).</div></div>
<div class="fig">{f("hotspots_genomewide.png","85%")}
<div class="cap">Genome-wide robust hotspots coloured by class. 6,322 hotspots: 5,856 are 45S rDNA; of the 466 non-rDNA,
<strong>~90% are pericentromeric</strong> (notably a Chr2 ~8.8&ndash;9.3 Mb cluster), 34 centromere, 15 arm.</div></div>
<div class="fig">{f("hotspots_venn.png","80%")}</div>
<div class="warn"><strong>Flat-Poisson vs gDNA-baseline hotspots.</strong> The standard hotspots above use a flat
genome-mean null (liberal). Re-calling against the <strong>WGA gDNA per-window expectation</strong> (DSB-specific,
mappability+copy-number corrected) gives 7,658 robust hotspots of which <strong>7,524 (98%) are 45S rDNA</strong> and only
77 pericentromere / 38 centromere / 19 arm. So calling vs gDNA <em>increases</em> rDNA (genuinely 12× enriched) and
<em>collapses</em> the pericentromeric hotspots &mdash; those were coverage-driven (more break events than the flat average,
but not more than their own gDNA coverage predicts). <strong>The rDNA is essentially the only DSB-specific hotspot locus.</strong></div>
<div class="fig">{f("hotspots_classification.png","60%")}
<div class="cap">3-way hotspot overlap + compartment classification. rDNA dominates raw hotspots (expected, per §4).</div></div>

<h3>8a. Where hotspots fall &mdash; and how that reconciles with the density test (§7d)</h3>
<p>The per-feature <em>density</em> test (§7d) found TEs are <em>not</em> more break-dense than genes (median TE &asymp; gene,
and TEs are <em>depleted</em> on the arms). Yet hotspot <em>counts</em> tell the opposite-looking story &mdash; and both are
correct, because <strong>hotspots are the upper tail of the density distribution, not its middle.</strong></p>
<div class="fig">{f("hotspot_distribution.png","99%")}
<div class="cap"><strong>Distribution of the 466 genuine (non-rDNA) robust hotspots.</strong> <b>A:</b> ~90% are pericentromeric
(417 peri, 34 centromere, 15 arm). <b>B:</b> occupancy fold (%hotspots &divide; %genome): <strong>TEs 1.67&times; enriched</strong>,
genes 0.10&times;, lncRNA 0.24&times;, satellites 0.09&times; (satellites sit in unmappable centromere cores, so hotspots cannot
be called on them). <b>C:</b> the tail metric &mdash; <strong>0.92% of TEs host a hotspot vs only 0.06% of genes (~15&times;)</strong>.</div></div>
<div class="key"><strong>Reconciliation.</strong> A <em>typical</em> transposon is no more break-prone than a gene (§7d), but the
<em>extreme</em> break sites fall overwhelmingly on a small subset of (mostly pericentromeric) TEs. So &ldquo;DSBs are enriched
in TEs&rdquo; is a statement about the <strong>tail / occupancy</strong> (a minority of TEs are strong hotspots), not about the
<strong>median</strong> TE. Density measures the middle of the distribution; hotspot counts and breaks/Mb measure the tail &mdash;
which is why they disagree. <em>Script:</em> <code>hotspot_distribution.py</code>.</div>

<h2 id="s9">9. Centromere satellite (CEN180) phasing — inconclusive</h2>
<p>We tested whether DSBs are phased to the 178 bp CEN180 monomer two ways: (a) a whole-genome metaplot over the
66,131 genomic CEN180 copies, and (b) mapping the sBLISS reads against the 356 bp <strong>CEN180 dimer</strong>
(2&times;178 bp consensus, <code>minimap2 -ax sr</code>) to read off a 5&prime;-cut-site profile across one repeat unit.</p>
<div class="fig">{f("cen178_metaplot_bliss.png","88%")}
<div class="cap"><strong>(a) Whole-genome CEN180 metaplot</strong> (anchor = repeat start, &plusmn;1068 bp &asymp; 12 monomers,
per-window detrended). The DSB signal is <strong>~178 bp-periodic</strong> across the satellite (autocorrelation at 178 bp
&asymp; +0.73, all three libraries) &mdash; DSB density tracks the monomer, non-sinusoidally.</div></div>
<div class="fig">{f("cen180_phasing_main_bliss.png","82%")}
<div class="cap"><strong>(b) sBLISS 5&prime;-end profile on the CEN180 dimer</strong> vs MNase nucleosome occupancy (grey). All three
libraries agree on a structured profile across the 356 bp unit &mdash; but see the normalisation below before reading
biology into it.</div></div>
<div class="fig">{f("cen180_phasing_main_bliss_wganorm.png","82%")}
<div class="cap"><strong>Why the dimer is confounded.</strong> (A) raw DSB 5&prime; density, the WGA (random) gDNA baseline
(dotted) and MNase; (B) DSB &divide; gDNA. The single consensus dimer imposes a <strong>sequence-similarity coverage
profile</strong>: even WGA <em>naked</em> DNA tracks nucleosome occupancy (r = +0.94) in the monomer interior, so the
gDNA &ldquo;baseline&rdquo; already carries the signal it is meant to remove and the normalised DSB&ndash;nucleosome
correlation <strong>flips sign with the window</strong>.</div></div>
<div class="warn"><strong>Verdict: inconclusive.</strong> There is a real, reproducible ~178 bp DSB periodicity over the
centromeric satellite (a), but the dimer profile (b) cannot be cleanly de-confounded: no robust DSB&ndash;nucleosome
phase can be claimed, and short-read sBLISS (~2&ndash;4% of the satellite uniquely mappable) cannot resolve
sub-nucleosome positioning in CEN180. Shown here for completeness, not as a positive result.</div>

<h2 id="s10">10. BiCroLab blissNPanalysis cross-checks</h2>
<p>Mirroring <a href="https://github.com/BiCroLab/blissNPanalysis">BiCroLab/blissNPanalysis</a> (human + mouse), adapted to TAIR12.</p>
<div class="fig">{f("bicrolab_composition.png")}
<div class="cap">Promoter/gene/intergenic + biotype (BA1 ~22%/63%/15%; of genic DSBs ~66% are rRNA = rDNA).</div></div>
<div class="fig">{f("bicrolab_circos.png","60%")}<div class="cap">Circular genome density (1 Mb, z-capped mean+3SD).</div></div>
<div class="fig">{f("bicrolab_pearson_heatmap_100kb.png","42%")} {f("bicrolab_pearson_heatmap_2kb.png","42%")}
<div class="cap">Pearson correlation, 100 kb (human) and 2 kb (mouse) windows.</div></div>
<div class="fig">{f("bicrolab_2kb_distribution.png","55%")}<div class="cap">DSBs per 2 kb window (ECDF) — long tail of hot windows.</div></div>

<h2 id="s11">11. What the hotspots are — external validation &amp; epigenomic context</h2>
<p>To ask <em>what kind of loci</em> the ssBLISS hotspots are, we remapped published Col-0 datasets to <strong>our exact TAIR12</strong>
(identical pipeline; multireads kept, MAPQ&ge;0) and intersected them with the hotspots: WT &gamma;H2A.X (this study + published SRR8434381),
the HR/resection markers <strong>RAD51 / RPA1A</strong> (SRR24938940/935), <strong>R-loops</strong> (ssDRIP, SRR8097477),
<strong>G-quadruplexes</strong> (G4-seq OQS maps, GSE110582, lifted TAIR10&rarr;TAIR12), and <strong>DNA methylation</strong> (Col-CC mCG/mCHG/mCHH).</p>

<h3>11a. Cross-validation with DSB / repair markers</h3>
<div class="fig">{f("bliss_vs_dsb_markers.png","85%")}
<div class="cap">Arm 10 kb Spearman of each marker vs BLISS DSBs. All weak, but ordered by proximity to the break:
<strong>RAD51 (+0.11) &asymp; published &gamma;H2A.X&divide;H3 (+0.10) &gt; RPA1A (+0.07) &gt; our &gamma;H2A.X (~0)</strong>.</div></div>
<div class="fig">{f("rad51_at_bliss_hotspots.png","98%")}
<div class="cap"><strong>Focal test at the hotspots:</strong> RAD51&divide;input and RPA1A&divide;input are <strong>significantly enriched at ssBLISS hotspots</strong>
(RAD51 p=0.0019, RPA1A p=0.0027; BLISS positive control peaks sharply) &mdash; the homologous-recombination machinery marks them,
orthogonally validating the hotspots as genuine repair-engaged DSB sites.</div></div>
<div class="fig">{f("h2ax_log2_q10.png","85%")}
<div class="cap"><strong>&gamma;H2A.X (WT) log2(IP/INPUT), MAPQ&ge;10.</strong> Flat across the unique-mappable arms with no focal peaks (centromere/NOR greyed = mappability blind spot).
In unstressed WT, &gamma;H2A.X is diffuse, so it does <em>not</em> co-localise with BLISS hotspots (p=0.08, ns) &mdash; the expected biology.
Raw &gamma;H2A.X coverage simply tracks H3/nucleosomes (&rho;+0.90), which is why the H3 control was essential.</div></div>

<h3>11b. R-loops and G-quadruplexes</h3>
<div class="fig">{f("bliss_vs_repair_rloop_compartment.png","85%")}
<div class="cap">By compartment (multireads kept): <strong>R-loops (ssDRIP) are the strongest BLISS correlate everywhere</strong>
(pericentromere +0.31, rDNA +0.25, centromere +0.24). The (peri)centromere is where physiological DSBs, R-loops and the HR machinery co-enrich.</div></div>
<div class="key"><strong>G-quadruplexes (G4-seq OQS).</strong> Using the experimental G4 maps (not coverage): non-rDNA ssBLISS hotspots are
<strong>~4&times; enriched for G-quadruplexes</strong> (stable-K 4.4&times;, PDS 4.2&times;; 36 vs ~8.6 expected, significant), although the
<em>average</em> G4 site has background DSB density &mdash; so G4s mark a <strong>subset</strong> of strong hotspots (fragile/replication-stress loci),
not a genome-wide DSB predictor.</div>

<h3>11c. Hotspot classification</h3>
<div class="fig">{f("hotspot_taxonomy.png","98%")}
<div class="cap">Each ssBLISS hotspot classified by mark (positive = above the 75th pct of compartment-matched random; 25% = chance).
<strong>Arm hotspots</strong> = transcription-associated (60% R-loop&#8314;, 68% genic) and HR-repaired (RAD51/RPA1A&#8314;).
<strong>Pericentromeric hotspots</strong> = R-loop&#8314; (86%) and &gamma;H2A.X&#8314; (53%) but RAD51-depleted &mdash; transcription/R-loop breaks <em>not</em> HR-repaired.
rDNA hotspots = unmarked within the array. (Table: <code>hotspots/hotspot_taxonomy.tsv</code>.)</div></div>

<h3>11d. DNA methylation &mdash; DSB-bearing TEs are hypermethylated</h3>
<div class="fig">{f("TE_hotspot_methylation.png","90%")}
<div class="cap">TEs carrying a DSB hotspot are <strong>hyper-methylated</strong> (mCG 70 vs 58%, mCHG 38.5 vs 26.4% p=5e-33).
Compartment-controlled, the signal is real on the <strong>arms</strong>: hotspot-TEs mCG <strong>52.5% vs 29%</strong> (p=2e-6); pericentromeric TEs are
saturated (no difference). So on the arms it is the <strong>silenced/heterochromatic (methylated) TEs</strong> that acquire DSBs &mdash; a route
distinct from the transcription/R-loop class. BLISS&times;methylation 10 kb &rho;: pericentromere +0.26&ndash;0.31, arms ~+0.06.</div></div>

<h3>11e. Genic hotspots sit on internal repeats; lncRNAs</h3>
<div class="fig">{f("lncRNA_and_tandemrepeat_genes.png","90%")}
<div class="cap"><strong>Left:</strong> lncRNA loci (separate from coding genes) are modestly DSB-dense (1.11&times; coding, between genes and TEs).
<strong>Right:</strong> DSB-hotspot genes are <strong>~26&times;</strong> more likely to contain an internal tandem repeat than a random gene (26% vs 1%);
top-DSB genes OR 3.6 (p=6e-11). GO of the repeat genes: <strong>cell-wall structural proteins</strong> (extensins/AGPs, p=1.7e-27) + chromatin/heterochromatin regulators.</div></div>
<div class="fig">{f("dotplot_AT4G22495.png","46%")} {f("tandem_repeat_lengths.png","52%")}
<div class="cap"><strong>Left:</strong> self dot-plot of AT4G22495 &mdash; a "gene hotspot" that is really an AthREP48 satellite array (parallel diagonals).
<strong>Right:</strong> tandem-repeat anatomy: 83% are the 178 bp centromeric satellite; arrays from ~0.5 kb to the 1.6 Mb centromere cores; gene-internal coding repeats are short (28&ndash;97 bp).</div></div>

<div class="key"><strong>Hotspot calling, flat-Poisson vs gDNA-baseline.</strong> Re-calling against the local WGA-gDNA expectation collapses the
<strong>non-rDNA hotspots ~74% (466 &rarr; ~123)</strong> &mdash; most flat-Poisson pericentromeric calls were coverage-driven &mdash; while rDNA calls
<em>rise</em> (genuinely ~12&times; over gDNA). The gDNA-baseline set is the conservative "real DSB-specific" list.</div>

<div class="note"><strong>Synthesis.</strong> ssBLISS DSB hotspots are not random or function-defined; they are <strong>fragile structural loci</strong> arising by (at least) two routes:
(i) <strong>transcription&ndash;replication conflict</strong> on the arms (R-loop&#8314;, G4&#8314;, genic, HR-repaired by RAD51/RPA1A), and
(ii) <strong>heterochromatin fragility</strong> at methylated TEs / pericentromeric repeats. The rDNA is the one dominant gDNA-confirmed hotspot.
&gamma;H2A.X (unstressed WT) is diffuse and does not pinpoint them.</div>

<h2 id="s12">12. Methods &amp; caveats</h2>
<ul>
<li><strong>Regions (curated Col-CC):</strong> centromere/pericentromere remapped CP116280.1&rarr;Chr1…; 45S NOR Chr2:0&ndash;5.41/Chr4:0&ndash;3.85 Mb; 5S Chr3/Chr4; arms = complement.</li>
<li><strong>Mappability:</strong> GenMap k=127 E=2; MAPQ&ge;1 (minimap2) / MAPQ&ge;10 (bowtie2) for unique-read views.</li>
<li><strong>gDNA:</strong> WGA/random (E-MTAB-6257 ERR2215864/865), mapped identically (minimap2 + bowtie2), used as the DNA-content/mappability baseline.</li>
<li><strong>Numerator = distinct break sites</strong> unless stated; per-1k-reads and per-mapped-Mb are proportional.</li>
<li><strong>External datasets (&sect;11):</strong> &gamma;H2A.X (SRR8434381), RAD51/RPA1A/INPUT (SRR24938940/935/942), ssDRIP R-loops (SRR8097477) remapped to TAIR12 (bowtie2/minimap2 -ax sr, dedup, CPM); G4-seq OQS peaks (GSE110582) lifted TAIR10&rarr;TAIR12 by sequence-remapping (~99.7%); methylation = Col-CC mCG/mCHG/mCHH 10 kb. <em>Caveat:</em> the RAD51/RPA1A reads came from a stringent-MAPQ pipeline so their repeat coverage was re-derived MAPQ&ge;0; G4-seq <em>coverage</em> is not G4 signal (only the OQS peak calls are used).</li>
<li><strong>Caveats:</strong> repeats (centromere, rDNA) unquantifiable per-unique-read (2&ndash;4% mappable); 8 bp UMI saturates at rDNA; satellite sub-nucleosome phasing out of reach (§9); expression-stratified gene metagenes skipped (no tissue-matched Col-0 TPM); arm/centromere hotspot-class counts are small (n=40&ndash;76) so effect sizes are indicative.</li>
</ul>
</body></html>"""
with open(OUT,"w") as fh: fh.write(html)
print("Wrote",OUT,f"({os.path.getsize(OUT)//1024} KB)")
