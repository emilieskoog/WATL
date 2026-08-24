#!/usr/bin/env bash
set -euo pipefail

BAMDIR="/scratch/eskoog/pI_paper_AMG_story/plasmid_derep/per_sample_bams_plasmid"
OUTDIR="/scratch/eskoog/pI_paper_AMG_story/plasmid_derep/coverm"

mkdir -p "$OUTDIR"

coverm contig \
    --bam-files "$BAMDIR"/*.sorted.bam \
    --methods rpkm covered_fraction \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 85 \
    --min-covered-fraction 0.50 \
    --threads 32 \
    --output-file "$OUTDIR/coverm_plasmid_multi_ALL.tsv"