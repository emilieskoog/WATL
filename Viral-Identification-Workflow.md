---
title: Viral Identification

---

## Viral Identification

This document outlines a complete viral genome recovery and curation pipeline for metagenomic assemblies generated from the WATL datasets. It includes:

- Viral detection using independent classifiers (geNomad and VirSorter2)
- Provirus identification and extraction
- Dereplication of viral contigs
- Genome completeness and contamination assessment with CheckV
- Quality filtering based on length and gene count
- Generation of a curated, non-redundant viral genome set
- Preparation of viral genomes for functional annotation (DRAM-V) and isoelectric point analysis

### Run both geNomad and VirSorter2 on WATL assemblies

## 1a. 🔍 Identify viruses with `Genomad`


```bash

#!/bin/bash

genomad end-to-end --cleanup /scratch/eskoog/WATL24_metagnomes/megahit_coassembly_out/final.contigs.fa genomad_WATL24_megahit_coassembly_output /home/eskoog/genomad_db

```
## 1b. Combine all viruses identified by `Genomad`

```
cat /scratch/eskoog/WATL24_metagnomes/genomad_WATL24_megahit_coassembly_output/final.contigs_find_proviruses/final.contigs_provirus.fna \
    /scratch/eskoog/WATL24_metagnomes/genomad_WATL24_megahit_coassembly_output/final.contigs_summary/final.contigs_virus.fna \
    > /scratch/eskoog/WATL24_metagnomes/genomad_WATL24_megahit_coassembly_output/all_genomad_viruses_and_proviruses.fna
    
        
```
## 1c. 🧱 Rename duplicate entries before  `CheckV`

```
awk '/^>/ {count[$0]++; if (count[$0]>1) {$0 = $0"_"count[$0]}; print; next} 1'   all_genomad_viruses_and_proviruses.fna > all_genomad_viruses_and_proviruses_dedup.fna
```
## 1d. 🔍 Identify viruses with `Virsorter2`

## 1e. Tag contigs before merging

> Before combining FASTA files from geNomad and VirSorter2, prefix each contig ID to retain source info:

```
# Add prefix to geNomad contigs for WATL24
sed 's/^>/\>genomad_/' /scratch/eskoog/WATL24_metagnomes/genomad_WATL24_megahit_coassembly_output/all_genomad_viruses_and_proviruses.fna > /scratch/eskoog/WATL24_metagnomes/genomad_WATL24_megahit_coassembly_output/all_genomad_viruses_and_proviruses_tagged.fna

# Add prefix to Virsorter2 contigs for WATL24

sed 's/^>/\>virsorter2_/' /scratch/eskoog/WATL24_metagnomes/virsorter2_output_from_kbase/WATL_24_final-viral-combined.fa > /scratch/eskoog/WATL24_metagnomes/virsorter2_output_from_kbase/WATL_24_final-viral-combined_virsorter2_tagged.fna


```


## 2. Concatenate all .fna files (with unique IDs).

```
mkdir /scratch/eskoog/WATL24_metagnomes/combined_genomad_virsorter2_viral_sequences

cat /scratch/eskoog/WATL24_metagnomes/genomad_WATL24_megahit_coassembly_output/all_genomad_viruses_and_proviruses_tagged.fna /scratch/eskoog/WATL24_metagnomes/virsorter2_output_from_kbase/WATL_24_final-viral-combined_virsorter2_tagged.fna > /scratch/eskoog/WATL24_metagnomes/combined_genomad_virsorter2_viral_sequences/WATL24_genomad_virsorter2_combined_viruses.fna

mkdir /scratch/eskoog/Tank_Exp/combined_genomad_virsorter2_viral_sequences


```


## 3. Dereplicate using `CD-HIT-EST` at 95% identity (common threshold for viral clustering)

> The -d 0 keeps full headers (important!).
-i: input FASTA file
-o: output FASTA (representative sequences)
-c: identity threshold (e.g., 0.95 = 95%)
-n: word size (must match -c; e.g., 10 for 0.95)
-T: number of threads
-M: memory limit in MB

```
cd-hit-est -i /scratch/eskoog/WATL24_metagnomes/combined_genomad_virsorter2_viral_sequences/WATL24_genomad_virsorter2_combined_viruses.fna -o /scratch/eskoog/WATL24_metagnomes/combined_genomad_virsorter2_viral_sequences/WATL24_genomad_virsorter2_combined_dereplicated_viruses.fna -c 0.97 -n 10 -d 0 -T 64 -M 1000

```

## 4. Check viral quality with `CheckV`

```
checkv end_to_end /scratch/eskoog/WATL24_metagnomes/combined_genomad_virsorter2_viral_sequences/WATL24_genomad_virsorter2_combined_dereplicated_viruses.fna checkV_output/ -t 60

```

## 5. Filter for quality viral sequences

### ✅ Filter 1: Typical viruses 
>Criteria:
>* checkv_quality ∈ Medium-quality, High-quality, Complete
>* contig_length ≥ 10000
>* gene_count ≥ 5


```

 awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || 
(($8 == "Medium-quality" || $8 == "High-quality" || $8 == "Complete") && 
 $2 >= 10000 && $5 >= 5)' quality_summary.tsv > filtered_typical.tsv


```

### ✅ Filter 2: Atypical viruses
>Criteria:
>* checkv_quality ∈ Low-quality, Medium-quality, High-quality, Complete
>* 3500 < contig_length < 10000
>* gene_count ≥ 5


```
 
 awk -F'\t' 'BEGIN {OFS="\t"} NR==1 || 
(($8 == "Low-quality" || $8 == "Medium-quality" || $8 == "High-quality" || $8 == "Complete") && 
 $2 > 3500 && $2 < 10000 && $5 >= 5)' quality_summary.tsv > filtered_atypical.tsv
 ```

### Merge quality-filtered typical and atypical viruses

```
head -n 1 filtered_typical.tsv > genomad_virsorter2_dereplicated_checkv_quality_filtered_typical_atypical_viruses.tsv
tail -n +2 filtered_typical.tsv >> genomad_virsorter2_dereplicated_checkv_quality_filtered_typical_atypical_viruses.tsv
tail -n +2 filtered_atypical.tsv >> genomad_virsorter2_dereplicated_checkv_quality_filtered_typical_atypical_viruses.tsv

```

## 6. Filter for quality viral sequences

### Merge the proviruses.fna and viruses.fna files outputted by checkV

```
cat proviruses.fna viruses.fna > all_viruses_after_checkv_not_quality_filtered.fna
```


### Get contig IDs of quality filtered viruses


```
cut -f1 genomad_virsorter2_dereplicated_checkv_quality_filtered_typical_atypical_viruses.tsv > filtered_ids.txt
```

### Extract contigs from file containing all viral sequences using the list contig IDs

#### 🔹 Create a version of the FASTA with cleaned headers
```
sed -E 's/^>(.*)_1 .*/>\1/' all_viruses_after_checkv_not_quality_filtered.fna > all_viruses_after_checkv_not_quality_filtered_cleaned_headers.fna
```
#### 🔹 Match exactly with seqkit

```
seqkit grep -f filtered_ids.txt all_viruses_after_checkv_not_quality_filtered_cleaned_headers.fna \
  -o all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2.fna
```

#### 🔹 Confirm recovery of all viral sequences

The following two should be the same/ off by 1 count (the first line of filtered_ids.txt is a header and not a viral sequence)
```
grep -c "^>" all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2.fna
```
```
wc -l filtered_ids.txt
```



