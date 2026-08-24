#!/bin/bash
# run_host_dram.sh — annotate all WATL MAGs with DRAM (KOfam+Pfam+dbCAN, no UniRef)
set -euo pipefail

MAGDIR="/scratch/eskoog/WATL_wet_dry_MAGs"
OUT="/scratch/eskoog/pI_paper_AMG_story/host_dram"
THREADS=50
NSHARDS=5

mkdir -p "$OUT"
cd "$OUT"

if [ ! -d shards ]; then
  mkdir -p shards
  i=0
  for f in "$MAGDIR"/*.fa; do
    s=$(( i % NSHARDS )); mkdir -p "shards/shard_${s}"
    ln -sf "$f" "shards/shard_${s}/"
    i=$(( i + 1 ))
  done
  echo "split $i MAGs into $NSHARDS shards"
fi

for s in $(seq 0 $(( NSHARDS - 1 ))); do
  if [ -f "shard_${s}_out/annotations.tsv" ]; then
    echo "shard_${s} already done, skipping"; continue
  fi
  rm -rf "shard_${s}_out"          # DRAM refuses to write into an existing dir
  echo "=== annotating shard_${s} $(date) ==="
  DRAM.py annotate \
    -i "shards/shard_${s}/*.fa" \
    -o "shard_${s}_out" \
    --min_contig_size 2500 \
    --threads "$THREADS" \
    2>&1 | tee "dram_shard_${s}.log"
done

echo "=== merging $(date) ==="
DRAM.py merge_annotations -i "shard_*_out" -o merged
echo "DONE -> $OUT/merged/annotations.tsv"
wc -l "$OUT/merged/annotations.tsv"
