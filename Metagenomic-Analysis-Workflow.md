

# 🧬 Metagenome Analysis Workflow Co-assembly

This document outlines a complete metagenome analysis pipeline for Illumina paired-end data collected for the Tank Experiment. It includes:

- Quality filtering
- Assembly
- Coverage mapping
- Binning
- Genome quality assessment
- Abundance estimation

---

## 1. 🔍 Quality Filtering with `fastp`
>This step removes low-quality reads and trims adapters using fastp, a fast and versatile quality control tool for Illumina sequencing data. Filtering is performed on both forward and reverse reads to improve the accuracy and reliability of downstream steps like assembly and binning. The filtered reads are retained in paired-end format, and summary reports are generated for quality assessment.


```bash
#!/bin/bash
set -e

INPUT_DIR="/data_store/seq_data/igm_20250502"

for R1 in ${INPUT_DIR}/*_SE_*_R1_001.fastq.gz; do
  R2="${R1/_R1_/_R2_}"
  if [[ ! -f "$R2" ]]; then
    echo "Missing R2 pair for $R1 — skipping." >&2
    continue
  fi

  BASENAME=$(basename "$R1" _R1_001.fastq.gz)

  fastp -i "$R1" -I "$R2" \
        -o "${BASENAME}_filt_R1.fq.gz" -O "${BASENAME}_filt_R2.fq.gz" \
        -e 30 -w 1 \
        -j "fastp_${BASENAME}.json" -h "fastp_${BASENAME}.html" \
        --failed_out "${BASENAME}_failed.fq.gz"
done



```


## 2. Concatenate all  reads for co-assembly

```
cat *_filt_R1.fq.gz > all_filt_R1.fq.gz
cat *_filt_R2.fq.gz > all_filt_R2.fq.gz
```


## 3. 🧱 Co-Assembly with `Megahit`
>Note:
This step assembles the quality-filtered paired-end reads into longer contiguous sequences (contigs) using MetaSPAdes, a de Bruijn graph-based assembler optimized for metagenomic data. Because individual samples had low DNA input and relatively high read duplication, we performed a co-assembly using MetaSPAdes to improve genome recovery. Co-assembly increases the effective sequencing depth and read diversity available for assembly, especially for low-abundance genomes that may not be recoverable from individual samples. By pooling quality-filtered and deduplicated reads across all samples, we enabled the reconstruction of more complete and contiguous metagenome-assembled genomes (MAGs) from shared microbial populations.




```
#!/bin/bash

megahit -1 all_filt_R1.fq.gz -2 all_filt_R2.fq.gz -o megahit_coassembly_out -t 64 --min-contig-len 1000
```  


## 4. 📏 Assembly Quality Assessment with `MetaQUAST`
>Note:
This step uses MetaQUAST to evaluate the quality of each individual assembly. It calculates key assembly metrics such as N50, total length, number of contigs, and GC content. These metrics help assess whether the assembly is of sufficient quality for downstream analyses. While reference genomes are not used here, MetaQUAST still provides valuable reference-free statistics for comparative assessment across samples.
>* **N50**
*Definition: The length of the shortest contig at 50% of the total assembly length.
Interpretation: Higher N50 indicates longer contigs and generally better assembly continuity.
Caveat: Can be misleading if there are a few very long contigs but many fragmented ones.*

>* **Total Length**
*Definition: The sum of the lengths of all contigs in the assembly.
Interpretation: Reflects how much of the genome(s) was recovered.
Context: Useful for comparing completeness across samples.*

>* **Number of Contigs**
*Definition: The total number of contigs above a certain length threshold (often 500 bp).
Interpretation: Fewer contigs usually suggest a more contiguous (better) assembly. 
Context: High numbers may indicate a fragmented assembly or complex community.**

```bash
#!/bin/bash
metaquast.py megahit_coassembly_out/final.contigs.fa \
  --no-plots -o metaquast_output_megahit -t 64
```



## 5. 🗂 Index Combined Assembly with `BWA`
>Note:
Indexing the concatenated FASTA file with bwa index prepares it for efficient read alignment in the next step. This is necessary for mapping reads from all samples back to the combined assembly so that we can calculate per-contig coverage across samples which is a necessary input for binning tools like MetaBAT2.

```bash
#!/bin/bash

bwa index megahit_coassembly_out/final.contigs.fa
```

## 6. 🧾 Cross-Mapping & Coverage Calculation
>Note:
In this step, we align the filtered reads from each sample to the combined, deduplicated assembly using bwa mem. This is known as cross-mapping, where each sample's reads are mapped not just to their own contigs but to the entire pooled reference. We then use samtools to sort the resulting BAM files and jgi_summarize_bam_contig_depths to compute coverage statistics for each contig across all samples.

>This cross-sample coverage profile is a key input for binning, as tools like MetaBAT2 use both sequence composition and differential coverage across samples to group contigs into genome bins.
```bash
#!/bin/bash

set -e
CPUS=64
REFERENCE="megahit_coassembly_out/final.contigs.fa"
TOMAPSORT=""

for R1 in *_filt_R1.fq.gz; do
    SAMPLE=$(basename "$R1" _filt_R1.fq.gz)
    R2="${SAMPLE}_filt_R2.fq.gz"

    bwa mem -t "$CPUS" "$REFERENCE" "$R1" "$R2" | \
        samtools view -Sb - > "${SAMPLE}_map_megahit_coassembly.bam"

    samtools sort -@ "$CPUS" -o "${SAMPLE}_map_sorted_megahit_coassembly.bam" "${SAMPLE}_map_megahit_coassembly.bam"
    TOMAPSORT="${TOMAPSORT} ${SAMPLE}_map_sorted_megahit_coassembly.bam"
    rm -f "${SAMPLE}_map_megahit_coassembly.bam"
done

jgi_summarize_bam_contig_depths \
  --outputDepth crossmap_contig_depth_megahit_coassembly.txt $TOMAPSORT

```

## 7. 🧱 Binning with `MetaBAT2`
>Note:
This step uses MetaBAT2 to group assembled contigs into metagenome-assembled genomes (MAGs). It leverages both tetranucleotide sequence composition and differential coverage profiles (from the cross-mapped samples) to cluster contigs that likely originate from the same genome. Binning enables the reconstruction of draft microbial genomes directly from metagenomic data, which can then be assessed for quality, taxonomy, and functional potential. Only contigs ≥1500 bp are considered to improve binning accuracy and reduce false positives.

```bash
#!/bin/bash

metabat2 \
  -i megahit_coassembly_out/final.contigs.fa \
  -a crossmap_contig_depth_megahit_coassembly.txt \
  -o bins_dir_megahit_coassembly/bin \
  --seed 1234 -m 1500 -v
```

## 8. ✅ Genome Quality Check with `CheckM2`
> Note:  
> This step uses `CheckM2` to evaluate the quality of each MAG based on lineage-specific marker genes. It estimates two key metrics:  
> 
> - **Completeness**: The percentage of expected marker genes that are present in the bin  
> - **Contamination**: The percentage of duplicated marker genes, which may indicate the bin contains sequences from multiple organisms  
> 
> These values are used to assess whether each MAG is suitable for downstream analysis. Bins with ≥50% completeness and ≤10% contamination are typically considered **medium- to high-quality**.

```bash
#!/bin/bash

checkm2 predict -i bins_dir_megahit_coassembly/ -x fa -o checkm2_output_megahit_coassembly -t 64
```

## 9. 🧹 Filter High-Quality Bins
> Note:  
> This step filters the results from `CheckM2` to retain only **high- or medium-quality genome bins**, defined here as those with:
>
> - **Completeness ≥ 50%**
> - **Contamination ≤ 10%**
>
> The filtered list of bin filenames is saved to `high_quality_bins.txt`, which is then used to extract those bins for downstream analysis (e.g., taxonomic classification, functional annotation, abundance estimation).
```bash
#!/bin/bash

awk 'NR > 1 && $2 >= 50 && $3 <= 10 { print $1 }' checkm2_output_megahit_coassembly/quality_report.tsv > high_quality_bins_megahit_coassembly.txt

sed -i 's/$/.fa/' high_quality_bins_megahit_coassembly.txt
```

## 10. 📁 Extract High-Quality Bins
> **Note:**  
> This step copies only the genome bins that passed the quality filter (as listed in `high_quality_bins.txt`) into a dedicated directory called `high_quality_bins/`. This simplifies downstream processing by isolating bins that are suitable for taxonomic classification, annotation, and abundance estimation.

```bash
#!/bin/bash

mkdir -p high_quality_bins_megahit_coassembly

while read BIN; do
    cp bins_dir_megahit_coassembly/"$BIN" high_quality_bins_megahit_coassembly/
done < high_quality_bins_megahit_coassembly.txt

```
## 11. Dereplicate High-quality bins
> Note:  Dereplicating MAGs with dRep removes redundant genomes by clustering similar bins based on average nucleotide identity (ANI), ensuring that each genome in the final set is unique and representative. This step improves downstream analyses by reducing bias from overrepresented strains and retaining the highest-quality version of each genome.
```
#!/bin/bash

dRep dereplicate drep_output_megahit_coassembly \
  -g high_quality_bins_megahit_coassembly/*.fa \
  -p 64 \
  --checkM_method taxonomy_wf \
  --completeness 50 \
  --contamination 10 \
  --S_algorithm fastANI \
  --clusterAlg average \
  --cov_thresh 0.10 \
  -sa 0.99
  ```
## 12. 📁 Contig ID Cleanup for CoverM Compatibility

After binning and extracting high-quality MAGs, we prepared the resulting genomes for abundance estimation with CoverM. However, CoverM requires that the contig IDs in the genome FASTA files **exactly match** those present in the BAM alignment files. However, **MetaBAT2 appends additional metadata** to each contig header, such as length, coverage, and sample depths, which means that CoverM will find no matches.

For example, a contig header from a MetaBAT2 bin might look like:

```
>sample1_contig0001	length=5000	cov=12.3
```

Whereas the BAM file expects:

```
@SQ	SN:sample1_contig0001	...
```

To fix this, we stripped everything after the first whitespace (space or tab) in each FASTA header using the following command:

```
sed -i 's/[ \t].*//' drep_output_megahit_coassembly/dereplicated_genomes/*.fa
```

This cleanup ensures that the contig IDs in the bin FASTAs match those in the BAM files, allowing CoverM to correctly assign reads to bins and calculate abundance metrics such as `relative_abundance`, `rpkm`, and `count`.

## 13. 📊 Abundance Estimation with `CoverM`
> Note:  
> This step uses `CoverM` to quantify the abundance of each high-quality genome bin across all samples. It calculates:
>
> - **Relative abundance**: Proportion of reads mapping to each bin per sample
> - **RPKM**: Reads per kilobase per million mapped reads, accounting for bin length and sequencing depth
> - **Count**: Raw number of reads mapping to each bin
>
> These metrics are useful for community composition analyses, normalization, and linking genome bins to environmental variables.


```bash
#!/bin/bash

coverm genome \
  --bam-files *_map_sorted_megahit_coassembly.bam \
  --genome-fasta-directory drep_output_megahit_coassembly/dereplicated_genomes\
  --genome-fasta-extension fa \
  --methods relative_abundance rpkm count \
  --min-covered-fraction 0 \
  --threads 64 \
  --output-file high_quality_derep_bin_abundance_megahit_coassembly.tsv
  
```

## 14. 🧬 Taxonomic Classification with `GTDB-Tk`

> **Note:**  
> This step uses `GTDB-Tk` to assign standardized taxonomy to high-quality genome bins based on genome-wide phylogenetic placement. GTDB-Tk provides species- and genus-level classifications that are consistent across studies, and it outputs a `classification.tsv` table with results for each bin.

```bash
#!/bin/bash

gtdbtk classify_wf \
  --genome_dir drep_output_megahit_coassembly/dereplicated_genomes \
  --out_dir gtdbtk_output_derep_HQ_bins \
  --extension fa \
  --cpus 64
```
---

## 📌 Notes

- Requires: `fastp`, `metaspades.py`, `metaquast.py`, `bwa`, `samtools`, `jgi_summarize_bam_contig_depths`, `metabat2`, `checkm2`, `coverm`
- Adjust paths/filenames to match your data organization.

```
#!/bin/bash

genomad end-to-end --cleanup /scratch/eskoog/WATL24_metagnomes/megahit_coassembly_out/final.contigs.fa genomad_WATL24_megahit_coassembly_output /home/eskoog/genomad_db
```
