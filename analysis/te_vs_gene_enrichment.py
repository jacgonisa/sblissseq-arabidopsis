#!/usr/bin/env python3
"""
Test whether transposons (TEs) are enriched for physiological DSBs RELATIVE TO genes.

The report already gives fold numbers vs the genome mean (genes ~0.85-0.90x, TEs
~1.2-1.6x) but never a formal significance test of TE vs gene. This does that,
controlling for the two known confounds:
  1. coverage / mappability  -> divide each feature's DSB density by WGA gDNA
     (naked-DNA) mean coverage over the same feature.
  2. compartment             -> TEs concentrate in the pericentromere, which is
     more break-dense; so we test genome-wide (nuclear, rDNA excluded) AND
     arm-restricted (features fully inside chromosome arms) separately. The
     arm-restricted + gDNA-normalised test is the clean headline.

Per feature: DSB density = (sum of break events) / (length/1000)  [breaks per kb].
gDNA-normalised density = DSB density / (mean WGA coverage over feature + eps).
Test: Mann-Whitney U (Wilcoxon rank-sum), TE vs gene, two-sided; effect size =
Cliff's delta (= 2*AUC - 1); plus a bootstrap CI on the median TE/gene ratio.

Outputs:
  results_TAIR12/analysis/te_vs_gene_enrichment.csv
  results_TAIR12/analysis/figures/te_vs_gene_enrichment.png

Env overrides: BLISS_ROOT, BLISS_REGDIR, BLISS_BWDIR, BLISS_BREAKDIR, BLISS_FIGDIR.
"""
import os, subprocess, tempfile
import numpy as np, pandas as pd, pyBigWig
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT   = os.environ.get("BLISS_ROOT",    "/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026")
ANA    = f"{ROOT}/results_TAIR12/analysis"
REG    = os.environ.get("BLISS_REGDIR",   f"{ANA}/regions")
BREAK  = os.environ.get("BLISS_BREAKDIR", f"{ROOT}/results_TAIR12/breaks")
BWWGA  = os.environ.get("BLISS_WGA",      f"{ANA}/coverage_cpm/wga.cov.cpm.bw")
FIG    = os.environ.get("BLISS_FIGDIR",   f"{ANA}/figures")
SAMPLES = ["BA1", "BA2", "old_BA1_BA2", "BA1_WT", "BA2_OX_40", "BA2_OX_80_1", "BA2_OX_80_2"]
GENO    = {"BA1":"WT","BA2":"WT","old_BA1_BA2":"WT","BA1_WT":"WT",
           "BA2_OX_40":"CENH3-OX","BA2_OX_80_1":"CENH3-OX","BA2_OX_80_2":"CENH3-OX"}
VIOLIN  = ["BA1", "BA2", "old_BA1_BA2"]   # the reliable WT set shown in the per-feature violins
NUCLEAR = {"Chr1","Chr2","Chr3","Chr4","Chr5"}
EPS = 1e-6
RNG = np.random.default_rng(0)

def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout

def prep_feature_bed(src, tmp, arm_bed, peri_bed, nor_bed):
    """Return {'all','arm','peri'} — nuclear, rDNA-excluded feature BEDs
       (cols chrom,start,end,name); 'arm'/'peri' = features lying ENTIRELY inside a
       chromosome-arm / pericentromere interval (compartment-stratified)."""
    base = f"{tmp}/{os.path.basename(src)}.clean.bed"
    # keep nuclear chroms, cut to 4 cols, sort; drop features overlapping the 45S NOR
    sh(f"awk 'BEGIN{{OFS=\"\\t\"}} $1~/^Chr[1-5]$/{{print $1,$2,$3,$4}}' {src} | sort -k1,1 -k2,2n "
       f"| bedtools intersect -a - -b {nor_bed} -v > {base}")
    arm  = f"{tmp}/{os.path.basename(src)}.arm.bed"
    peri = f"{tmp}/{os.path.basename(src)}.peri.bed"
    sh(f"bedtools intersect -a {base} -b {arm_bed}  -u -f 1.0 > {arm}")   # entirely within an arm
    sh(f"bedtools intersect -a {base} -b {peri_bed} -u -f 1.0 > {peri}")  # entirely within pericentromere
    return {"all": base, "arm": arm, "peri": peri}

def break_sum(feat_bed, break_bg, tmp):
    """Per-feature sum of break events via bedtools map."""
    out = f"{tmp}/map.txt"
    sh(f"bedtools map -a {feat_bed} -b {break_bg} -c 4 -o sum -null 0 > {out}")
    df = pd.read_csv(out, sep="\t", header=None,
                     names=["chrom","start","end","name","breaks"])
    df["breaks"] = pd.to_numeric(df["breaks"], errors="coerce").fillna(0.0)
    return df

def gdna_mean(df, bw):
    vals = np.empty(len(df))
    for i,(c,s,e) in enumerate(zip(df.chrom, df.start, df.end)):
        try:
            v = bw.stats(c, int(s), int(e), type="mean")[0]
        except Exception:
            v = None
        vals[i] = 0.0 if v is None else float(v)
    return vals

def cliffs_delta(a, b):
    # P(a>b) - P(a<b) via the U statistic (a = TE, b = gene)
    U, _ = mannwhitneyu(a, b, alternative="two-sided")
    return 2.0 * U / (len(a) * len(b)) - 1.0

def boot_ratio_ci(te, ge, n=1000):
    r = np.empty(n)
    for k in range(n):
        rt = np.median(RNG.choice(te, size=len(te), replace=True))
        rg = np.median(RNG.choice(ge, size=len(ge), replace=True))
        r[k] = rt / rg if rg > 0 else np.nan
    r = r[np.isfinite(r)]
    return np.percentile(r, 2.5), np.percentile(r, 97.5)

def analyse(dens_col, gdf, tdf, sample, comparison, norm):
    ge = gdf[dens_col].to_numpy(); te = tdf[dens_col].to_numpy()
    ge = ge[np.isfinite(ge)]; te = te[np.isfinite(te)]
    U, p = mannwhitneyu(te, ge, alternative="two-sided")
    mg, mt = np.median(ge), np.median(te)
    ratio = mt / mg if mg > 0 else np.nan
    lo, hi = boot_ratio_ci(te, ge)
    return dict(sample=sample, comparison=comparison, norm=norm,
                n_gene=len(ge), n_TE=len(te),
                median_gene=mg, median_TE=mt, ratio_TE_over_gene=ratio,
                ratio_CI_low=lo, ratio_CI_high=hi,
                cliffs_delta=cliffs_delta(te, ge), U=U, p_value=p)

def main():
    os.makedirs(FIG, exist_ok=True)
    bw = pyBigWig.open(BWWGA)
    rows, dist = [], {}   # dist[(sample,comparison)] = (gene_df, te_df) with density cols
    COMPS = [("genome_noNOR","all"), ("arms_only","arm"), ("pericentromere","peri")]
    with tempfile.TemporaryDirectory() as tmp:
        gsets = prep_feature_bed(f"{REG}/genes.bed6", tmp, f"{REG}/arms.bed", f"{REG}/pericentromere.bed", f"{REG}/nor_45s.bed")
        tsets = prep_feature_bed(f"{REG}/TEs.bed6",  tmp, f"{REG}/arms.bed", f"{REG}/pericentromere.bed", f"{REG}/nor_45s.bed")
        for sample in SAMPLES:
            bg = f"{BREAK}/{sample}.break.bedgraph"
            for comp, key in COMPS:
                g = break_sum(gsets[key], bg, tmp); t = break_sum(tsets[key], bg, tmp)
                for df in (g, t):
                    kb = (df.end - df.start) / 1000.0
                    df["dens_raw"] = df.breaks / kb
                    df["gdna"] = gdna_mean(df, bw)
                    df["dens_gdna"] = df.dens_raw / (df.gdna + EPS)
                dist[(sample, comp)] = (g, t)
                for norm, col in [("raw","dens_raw"), ("gdna_norm","dens_gdna")]:
                    rows.append(analyse(col, g, t, sample, comp, norm))
    bw.close()

    out = pd.DataFrame(rows)
    csv = f"{ANA}/te_vs_gene_enrichment.csv"
    out.to_csv(csv, index=False, float_format="%.5g")
    print("wrote", csv)
    with pd.option_context("display.width",200,"display.max_columns",20):
        print(out[["sample","comparison","norm","n_gene","n_TE","median_gene",
                   "median_TE","ratio_TE_over_gene","cliffs_delta","p_value"]].to_string(index=False))

    # ---- cross-compartment view: is the pericentromere fragile *because of* TEs? ----
    # compare each feature class ACROSS compartments (median gDNA-norm density).
    print("\ncross-compartment median gDNA-norm DSB density (per feature):")
    print(f"{'sample':12} {'gene_arm':>9} {'gene_peri':>10} {'TE_arm':>8} {'TE_peri':>8} "
          f"{'peri/arm(TE)':>13} {'peri/arm(gene)':>15}")
    for sample in SAMPLES:
        ga = np.median(dist[(sample,'arms_only')][0].dens_gdna)
        gp = np.median(dist[(sample,'pericentromere')][0].dens_gdna)
        ta = np.median(dist[(sample,'arms_only')][1].dens_gdna)
        tp = np.median(dist[(sample,'pericentromere')][1].dens_gdna)
        print(f"{sample:12} {ga:9.3f} {gp:10.3f} {ta:8.3f} {tp:8.3f} {tp/ta:13.2f} {gp/ga:15.2f}")

    # ---- figure: gDNA-normalised density, gene vs TE, per comparison x sample ----
    comps = [("genome_noNOR","Genome-wide (nuclear, rDNA excluded)"),
             ("arms_only","Chromosome arms only"),
             ("pericentromere","Pericentromere only")]
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 11.4), sharey=True)
    for r,(comp,ctitle) in enumerate(comps):
        for c,sample in enumerate(VIOLIN):
            ax = axes[r][c]; g,t = dist[(sample,comp)]
            data = [np.log10(g.dens_gdna.to_numpy()+1e-3), np.log10(t.dens_gdna.to_numpy()+1e-3)]
            parts = ax.violinplot(data, showmedians=True, widths=0.85)
            for pc,col in zip(parts['bodies'], ["#8c564b","#9467bd"]):
                pc.set_facecolor(col); pc.set_alpha(0.6)
            for kk in ("cmedians","cbars","cmins","cmaxes"):
                if kk in parts: parts[kk].set_color("#333")
            row = out[(out["sample"]==sample)&(out["comparison"]==comp)&(out["norm"]=="gdna_norm")].iloc[0]
            p = row.p_value; star = "ns" if p>=0.05 else ("*" if p>=1e-2 else ("**" if p>=1e-3 else "***"))
            arrow = "TE>gene" if row.ratio_TE_over_gene>1 else "TE<gene"
            ax.set_title(f"{sample}  (n {row.n_gene}/{row.n_TE})\nTE/gene={row.ratio_TE_over_gene:.2f} "
                         f"δ={row.cliffs_delta:+.2f} {star}  {arrow}",
                         fontsize=9, fontweight="bold")
            ax.set_xticks([1,2]); ax.set_xticklabels(["genes","TEs"])
            if c==0: ax.set_ylabel(f"{ctitle}\n\nlog10 gDNA-norm DSB density", fontsize=9)
            ax.grid(axis="y", alpha=0.25)
    fig.suptitle("DSB enrichment: transposons vs genes, by compartment (gDNA-normalised, per feature)\n"
                 "Mann-Whitney U;  δ = Cliff's delta (TE vs gene);  *** p<1e-3",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.965])
    png = f"{FIG}/te_vs_gene_enrichment.png"
    fig.savefig(png, dpi=140, bbox_inches="tight")
    print("wrote", png)

    # ---- genotype comparison: does CENH3-OX change the TE/gene ratio? ----
    # The TE/gene ratio is a WITHIN-sample relative measure, so it is largely robust
    # to the adapter/batch confound that plagues the raw OX-vs-WT landscape comparison.
    gd = out[out["norm"]=="gdna_norm"].copy()
    gd["geno"] = gd["sample"].map(GENO)
    comp_order = ["genome_noNOR","arms_only","pericentromere"]
    comp_lbl   = {"genome_noNOR":"Genome-wide","arms_only":"Arms","pericentromere":"Pericentromere"}
    fig2, ax = plt.subplots(figsize=(10.5, 5.6))
    for _, row in gd.iterrows():
        xi = comp_order.index(row["comparison"])
        wt = row["geno"]=="WT"
        xj = xi + (-0.16 if wt else 0.16)
        ax.errorbar(xj, row["ratio_TE_over_gene"],
                    yerr=[[row["ratio_TE_over_gene"]-row["ratio_CI_low"]],
                          [row["ratio_CI_high"]-row["ratio_TE_over_gene"]]],
                    fmt="o", ms=7, capsize=3, lw=1.2,
                    color=("#1f6feb" if wt else "#d62728"),
                    mfc=("#1f6feb" if wt else "#d62728"), alpha=0.9, zorder=3)
        ax.annotate(row["sample"].replace("BA2_OX_","OX").replace("BA1_WT","WT2").replace("old_BA1_BA2","old"),
                    (xj, row["ratio_CI_high"]), fontsize=6.5, ha="center", va="bottom", rotation=90, color="#444")
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xticks(range(len(comp_order))); ax.set_xticklabels([comp_lbl[c] for c in comp_order])
    ax.set_ylabel("TE / gene DSB-density ratio  (÷gDNA, per feature)")
    ax.set_title("Does CENH3 overexpression change transposon-vs-gene DSB enrichment?\n"
                 "WT (blue) vs CENH3-OX (red) — within-sample ratio, robust to batch/adapter", fontsize=11, fontweight="bold")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker="o",ls="",color="#1f6feb",label="WT (4 libraries)"),
                       Line2D([],[],marker="o",ls="",color="#d62728",label="CENH3-OX (3 libraries)")],
              frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.25); fig2.tight_layout()
    png2 = f"{FIG}/te_vs_gene_by_genotype.png"
    fig2.savefig(png2, dpi=140, bbox_inches="tight"); print("wrote", png2)

if __name__ == "__main__":
    main()
