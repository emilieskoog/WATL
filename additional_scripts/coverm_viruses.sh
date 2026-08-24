#!/bin/bash
set -euo pipefail

BAMDIR="per_sample_bams"
OUT="coverm_viruses_rpkm_mcf0.75_per_sample.tsv"

coverm contig \
    -t 50 \
    -b "$BAMDIR"/*.sorted.bam \
    -m rpkm \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 85 \
    --min-covered-fraction 0.75 \
    > "$OUT"