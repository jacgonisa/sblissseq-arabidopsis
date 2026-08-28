#!/usr/bin/env bash
set -euo pipefail

# sBLISS-seq analysis pipeline
# Run from: /mnt/ssd-8tb/NOVOGENE/BLISSseq_June2026/
# conda activate bliss  (or use: conda run -n bliss ...)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR/pipeline"
THREADS=8

cd "$PIPELINE_DIR"

conda run -n bliss snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores "$THREADS" \
    --rerun-incomplete \
    --printshellcmds \
    "$@"
