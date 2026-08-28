#!/usr/bin/env bash
# γH2AX broad-domain peak calling vs INPUT (consistent with BLISS gDNA-baseline method).
# Waits for INPUT_WT.dedup.bam, then: IP/INPUT log2 bigwig + binned Poisson domains.
set -euo pipefail
export PATH="$HOME/miniforge3/envs/bliss/bin:$PATH"
H=/mnt/ssd-8tb/NOVOGENE/H2AX_June2026/results_TAIR12
B=/mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/results_TAIR12
REG=$B/analysis/regions; OUT=$H/peaks; mkdir -p $OUT
IP=$H/bam/IP_WT.dedup.bam; IN=$H/bam/INPUT_WT.dedup.bam
echo "[wait] for $IN ..."; while [ ! -s "$IN" ] || [ ! -s "${IN}.bai" ]; do sleep 30; done; echo "[ok] INPUT ready"
# (1) INPUT-normalised log2 enrichment track (matched units, like the BLISS gDNA track)
bamCompare -b1 $IP -b2 $IN --operation log2 --scaleFactorsMethod None --normalizeUsing CPM \
  --binSize 1000 --pseudocount 1 --skipZeroOverZero -p 8 \
  -o $H/bigwig/IP_over_INPUT.log2.bw
# (2) broad-domain calling: 1kb windows, λ from depth-scaled INPUT, Poisson, merge<=5kb gaps
python - <<'PY'
import pyBigWig, numpy as np, subprocess, scipy.stats as st
H="/mnt/ssd-8tb/NOVOGENE/H2AX_June2026/results_TAIR12"
chrs=["Chr1","Chr2","Chr3","Chr4","Chr5"]; W=1000
import pysam
ip=pysam.AlignmentFile(f"{H}/bam/IP_WT.dedup.bam"); inp=pysam.AlignmentFile(f"{H}/bam/INPUT_WT.dedup.bam")
ipN=ip.mapped; inN=inp.mapped; sf=ipN/inN
rows=[]
for c in chrs:
    L=ip.get_reference_length(c); n=L//W
    ipc=np.array(ip.count_coverage(c,0,n*W,quality_threshold=0)).sum(0).reshape(n,W).sum(1)/150.0
    inc=np.array(inp.count_coverage(c,0,n*W,quality_threshold=0)).sum(0).reshape(n,W).sum(1)/150.0
    lam=np.maximum(inc*sf,0.5)
    p=st.poisson.sf(ipc-1,lam)
    for i in np.where((p<1e-4)&(ipc>=lam*2)&(ipc>=20))[0]:
        rows.append((c,i*W,(i+1)*W,float(ipc[i]),float(ipc[i]/lam[i]),float(p[i])))
import csv
with open(f"{H}/peaks/h2ax_enriched_1kb.bed","w") as f:
    for r in rows: f.write("%s\t%d\t%d\tw\t%.1f\t.\t%.2f\t%.2e\n"%r)
print("enriched 1kb windows:",len(rows))
PY
# merge into broad domains (<=5kb gap)
sort -k1,1 -k2,2n $H/peaks/h2ax_enriched_1kb.bed | bedtools merge -d 5000 -c 5,7 -o max,min \
  > $H/peaks/h2ax_domains.bed
echo "broad domains: $(wc -l < $H/peaks/h2ax_domains.bed)"
