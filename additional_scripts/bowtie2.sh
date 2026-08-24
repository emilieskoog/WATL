
#!/usr/bin/env bash

set -euo pipefail

# ===== inputs =====
REF="/scratch/eskoog/WATL_wet_dry_viral_analysis/checkV_output/all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2.fna"

# WET 2022 (curated brine list + reads dir)
BRINE_LIST="/scratch/eskoog/WATL_metagenomes_viral_work/brine_metagenomes.txt"
WET_DIR="/data_store/shared_work/metagenomes_watl/fastp_output"        # <sample>_R1_qc.fastq.gz / _R2_qc.fastq.gz

# DRY 2024 (already brine-only)
DRY_DIR="/scratch/eskoog/WATL24_metagnomes/only_brine_filt_files"      # <sample>_filt_R1.fq.gz / _filt_R2.fq.gz

OUTDIR="per_sample_bams"
THREADS=16

IDX="${REF%.*}_bt2idx"
mkdir -p "$OUTDIR"
[[ -r $REF ]] || { echo "cannot read REF: $REF" >&2; exit 1; }

# ---- build index once from the NEW catalog, reuse ----
if [[ -e "${IDX}.1.bt2" || -e "${IDX}.1.bt2l" ]]; then
    echo ">> reusing index $IDX" >&2
else
    echo ">> building bowtie2 index from $REF" >&2
    bowtie2-build --threads "$THREADS" "$REF" "$IDX"
fi

# ---- collect (sample, R1, R2) ----
SAMPLES=(); R1S=(); R2S=(); SKIPPED=()

# WET 2022: from the curated brine list
[[ -r $BRINE_LIST ]] || { echo "cannot read brine list: $BRINE_LIST" >&2; exit 1; }
nwet=0
while read -r f || [[ -n $f ]]; do
    f="${f%$'\r'}"; [[ -z $f ]] && continue          # drop trailing CR / blank lines
    r1="${WET_DIR}/${f}_R1_qc.fastq.gz"
    r2="${WET_DIR}/${f}_R2_qc.fastq.gz"
    if [[ -r $r1 && -r $r2 ]]; then
        SAMPLES+=("$f"); R1S+=("$r1"); R2S+=("$r2"); nwet=$((nwet+1))
    else
        SKIPPED+=("$f")
    fi
done < "$BRINE_LIST"
echo ">> WET 2022: $nwet sample(s) found from $(basename "$BRINE_LIST")" >&2

# DRY 2024: everything in the only-brine dir
[[ -d $DRY_DIR ]] || { echo "dry dir not found: $DRY_DIR" >&2; exit 1; }
ndry=0
while IFS= read -r -d '' r1; do
    s=$(basename "$r1" _filt_R1.fq.gz)
    r2="${r1/_filt_R1.fq.gz/_filt_R2.fq.gz}"
    [[ -r $r2 ]] || { echo "!! dry: missing mate for $r1 (expected $r2)" >&2; exit 1; }
    SAMPLES+=("$s"); R1S+=("$r1"); R2S+=("$r2"); ndry=$((ndry+1))
done < <(find "$DRY_DIR" -maxdepth 1 -type f -iname '*_filt_R1.fq.gz' -print0 | sort -z)
echo ">> DRY 2024: $ndry sample(s) found" >&2

if (( ${#SKIPPED[@]} > 0 )); then
    echo ">> WARNING: ${#SKIPPED[@]} wet sample(s) listed but no reads in $WET_DIR (skipped):" >&2
    printf '     %s\n' "${SKIPPED[@]}" >&2
fi
(( ${#SAMPLES[@]} > 0 )) || { echo "no read pairs found" >&2; exit 1; }
echo ">> ${#SAMPLES[@]} sample(s) total to map:" >&2
printf '   %s\n' "${SAMPLES[@]}" >&2

# ---- map each sample -> its own sorted, indexed BAM (resumable) ----
for i in "${!SAMPLES[@]}"; do
    s="${SAMPLES[$i]}"
    bam="${OUTDIR}/${s}.sorted.bam"
    if [[ -s "$bam" && -s "${bam}.bai" ]]; then
        echo ">> [$((i+1))/${#SAMPLES[@]}] $s already mapped, skipping" >&2
        continue
    fi
    echo ">> [$((i+1))/${#SAMPLES[@]}] mapping $s -> $bam" >&2
    bowtie2 -x "$IDX" -1 "${R1S[$i]}" -2 "${R2S[$i]}" -p "$THREADS" --no-unal \
        | samtools sort -@ "$THREADS" -o "$bam" -
    samtools index "$bam"
done
echo ">> done. per-sample BAMs in $OUTDIR/ (one per sample -> one RPKM column in CoverM)" >&2