# Identifying and Extracting Viral Structural Proteins Workflow

Viruses were previously identified using VirSorter2 and geNomad, then combined, dereplicated, and quality-filtered to generate a curated viral genome set as noted [here](https://github.com/emilieskoog/WATL/blob/main/Viral-Identification-Workflow.md). We now outline the workflow used to assign viral taxonomy, annotate viral genes, and extract viral structural proteins required for the subsequent viral proteome isoelectric point analysis.

## 1. 🔍 Identify viral taxonomy with `Genomad` (again)
> This time we are using the combine viruses from genomad and Virsorter2
> 
```
#!/bin/bash

genomad end-to-end --cleanup /scratch/eskoog/WATL24_metagnomes/checkV_output/all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2.fna genomad_WATL24_megahit_coassembly_already_derep_quality_filtered_viruses_output /home/eskoog/genomad_db
```

---

## 2.  Identify of many of your sequences were classified

```
wc -l /scratch/eskoog/WATL24_metagnomes/checkV_output/genomad_WATL24_megahit_coassembly_already_derep_quality_filtered_viruses_output/all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_summary/all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_summary.tsv
```


## 3. Predict protein-coding genes for each virus using `Prodigal-gv`
>Run viruses through Prodigal-gv adapted for viruses. (I also moved my file continaing all quality-filtered, dereplicated viruses re-taxonomically identified by Genomad to my home directory: /scratch/eskoog/WATL24_metagnomes/)
```
#!/bin/bash

# Define the directory containing the .fa file
DIR="/scratch/eskoog/WATL24_metagnomes/"

# Loop through each .fa file in the directory
for file in "$DIR"*.fna; do
    # Extract the filename without the extension
    base_name=$(basename "$file" .fna)

    # Run Prodigal with the specified options
    parallel-prodigal-gv.py -t 50 -i "$file" -a "${DIR}${base_name}.faa"
done
```

 ## 4. Annotate viral proteins with VOGdb and BLASTdb
>You will need to install VOGdb and NCBI RefSeq db

>Instaling NCBI RefSeq viral protein sequence database:
```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz

gunzip viral.1.protein.faa.gz

#Make a blastable database out of the viral_protein_db
makeblastdb -in viral.1.protein.faa -dbtype prot -out viral_protein_db
```

### BLASTdb
```
blastp \
-query .faa \
-db /home/eskoog/old_fram_dir/viral_protein_db \
-out all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_blast_db.tsv \
-evalue 1e-5 \
-qcov_hsp_perc 50 \
-outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp stitle" \
2>&1 | tee freshwater_NCBI_blastp_log.txt
```

### VOGdb
```
blastp \
-query .faa \
-db /home/eskoog/old_fram_dir/vogdb \
-out all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_db.tsv \
-evalue 1e-5 \
-qcov_hsp_perc 50 \
-outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp stitle" \
2>&1 | tee freshwater_VOG_log.txt
```
### Filter out hits with less than 30% identity 
```
awk -F'\t' '$3 >= 30' \
all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_blast_db.tsv \
> all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_blast_db_pident30.tsv

awk -F'\t' '$3 >= 30' \
all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_db.tsv \
> all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_vog_db_pident30.tsv

```

## 5. Get list of all identified genes by cleaning up the file and identifying unique genes of interest to look for

### a. Remove square brackets and their contents at the end of the last column
```
sed -E 's/\[.*\]$//' all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_blast_db_pident30.tsv >
all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_blast_db_2.tsv

sed -E 's/\[.*\]$//' all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_vog_db_pident30.tsv
 > all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_db_2.tsv
```

### b. remove the ID number from the last column that corresponds to the first column
```
awk 'BEGIN{FS=OFS="\t"} {gsub($2, "", $NF); print}' all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_blast_db_2.tsv > all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_blast_db_3.tsv

awk 'BEGIN{FS=OFS="\t"} {gsub($2, "", $NF); print}' all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_db_2.tsv > all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_db_3.tsv
```

### c. Get unique genes from tsv file
```
awk -F'\t' '{print $NF}' all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_blast_db_3.tsv | sort | uniq > unique_blast_gene_hits.txt

awk -F'\t' '{print $NF}' all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_db_3.tsv | sort | uniq > unique_vog_gene_hits.txt
```

### d. Manually pick out the genes relevant to viral structure including phage capsid, tail, portal, neck, and baseplate proteins 

>I have a file that I use called merged_list_of_viral_structural_genes_VOG_NCBI_RefSeq.txt

### e. Produce a file containing only the protein IDs whose best (lowest e-value) BLAST hit annotation matches an entry in merged_list_of_viral_structural_genes_VOG_NCBI_RefSeq.txt

```
sort -t$'\t' -k1,1 -k11,11g all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_blast_db_3.tsv > NCBI_blast_sorted.tsv

sort -t$'\t' -k1,1 -k11,11g all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_db_3.tsv > VOG_blast_sorted.tsv

awk -F'\t' '
NR==FNR {terms[$0]=1; next}
{
    pid=$1
    if (pid in seen) next
    seen[pid]=1

    ann=$NF
    gsub(/^[ \t]+|[ \t]+$/, "", ann)
    if (!(ann in terms)) next
    if (ann ~ /terminase|endolysin|nucleotide kinase|phosphofructokinase|DNA helicase|DNA methyltransferase|DNA polymerase|DNA repair exonuclease|HNH endonuclease|exonuclease|head maturation protease|head morphogenesis|tail fiber assembly|tail assembly|WhiB|tail sheath stabilizer|FhaB|peptidase|ferredoxin|lysozyme|putative coat protein|putative lysozyme|virion structural protein and packaging/) next
    print pid
}' \
merged_list_of_viral_structural_genes_VOG_NCBI_RefSeq.txt \
NCBI_blast_sorted.tsv \
> NCBI_structural_besthit.tsv

awk -F'\t' '
NR==FNR {terms[$0]=1; next}
{
    pid=$1
    if (pid in seen) next
    seen[pid]=1

    ann=$NF
    gsub(/^[ \t]+|[ \t]+$/, "", ann)
    if (!(ann in terms)) next
    if (ann ~ /terminase|endolysin|nucleotide kinase|phosphofructokinase|DNA helicase|DNA methyltransferase|DNA polymerase|DNA repair exonuclease|HNH endonuclease|exonuclease|head maturation protease|head morphogenesis|tail fiber assembly|tail assembly|WhiB|tail sheath stabilizer|FhaB|peptidase|ferredoxin|lysozyme|putative coat protein|putative lysozyme|virion structural protein and packaging/) next
    print pid
}' \
merged_list_of_viral_structural_genes_VOG_NCBI_RefSeq.txt \
VOG_blast_sorted.tsv \
> VOG_structural_besthit.tsv

```

### f. Combine the NCBI and VOG best-hit structural protein lists into a single deduplicated list

```
cat NCBI_structural_besthit.tsv VOG_structural_besthit.tsv | sort -u > all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_NCBI_RefSeq_VOG_3_pident30_viral_structural_genes1_unique_contigIDs_besthit.tsv

```

### g. Count the resulting unique protein IDs

```
wc -l all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_NCBI_RefSeq_VOG_3_pident30_viral_structural_genes1_unique_contigIDs_besthit.tsv
```

### e. Extract viral structural proteins of interest from .faa file using known (annotated) structural protein IDs
```
#!/bin/bash

# Define input file with contig IDs
contig_file="all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_NCBI_RefSeq_VOG_3_pident30_viral_structural_genes1_unique_contigIDs_besthit.tsv"

# Define output directory
output_dir="all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_blast_db_viral_structural_genes1_contigIDs_merged_viral_structural_genes_dir_besthit"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

echo "Reading contig IDs from $contig_file..."

# Load contig IDs into an associative array (trim spaces)
declare -A contig_ids
while read -r id; do
    id_trimmed=$(echo "$id" | tr -d ' ')  # Remove spaces
    contig_ids["$id_trimmed"]=1
done < "$contig_file"

echo "Processing .faa files..."
# Process each .faa file
for faa_file in all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2.faa; do
    echo "Processing: $faa_file"

    # Define new output filename
    new_faa_file="$output_dir/viral_structural_genes_${faa_file#quality_filtered_brine_viral_contigs}"

    # Ensure new file is empty before writing
    > "$new_faa_file"

    awk -v out="$new_faa_file" '
        BEGIN {
            # Load contig IDs into an array
            while ((getline < "'"$contig_file"'") > 0) {
                gsub(/^[ \t]+|[ \t]+$/, "", $1)  # Trim spaces
                contig[$1] = 1
            }
            close("'"$contig_file"'")
        }
        /^>/ {
            # Extract contig ID (everything before first #)
            match($0, /^>([^#]+)/, arr)
            id = arr[1]
            gsub(/^[ \t]+|[ \t]+$/, "", id)  # Trim spaces from extracted ID

            # Check if the contig ID matches the list
            keep = (id in contig) ? 1 : 0
        }
        keep { print > out }  # Print header and full sequence to output
    ' "$faa_file"

    # Remove empty files
    if [ ! -s "$new_faa_file" ]; then
        rm "$new_faa_file"
    else
        echo "Saved filtered sequences to: $new_faa_file"
    fi
done

echo "All files processed. Check $output_dir"
```

### f. Split .faa file into files of unique viruses with structural proteins
```
#!/bin/bash

input="all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_blast_db_viral_structural_genes1_contigIDs_merged_viral_structural_genes_dir_besthit/viral_structural_genes_all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2.faa"
output_dir="WATL_wet_dry_viruses_split_by_contig_faa_besthit"
mkdir -p "$output_dir"

awk -v outdir="$output_dir" '
    /^>/ {
        # Extract contig ID by removing final _digit(s) at end of ID
        match($0, /^>([^ ]+)/, a)         # Get full ID before first space
        full_id = a[1]
        contig = full_id

        # Remove the final _digit part (gene ID) ONLY if it exists after contig ID
        sub(/_[0-9]+$/, "", contig)

        # Make filename safe by replacing pipe
        safe_id = contig
        gsub(/\|/, "_", safe_id)
        outfile = outdir "/" safe_id ".faa"
    }
    {
        print >> outfile
    }
' "$input"
```


> You can then confirm that in directory `split_by_contig_faa`
```
find . -name "*.faa" -exec grep -h "^>" {} + | wc -l
```
> is the same as `viral_structural_genes_all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus.faa` file 
```
grep ">" -c viral_structural_genes_all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus.faa
```
> in directory `/scratch/eskoog/WATL24_metagnomes/all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2_virus_vog_blast_db_viral_structural_genes1_contigIDs_merged_viral_structural_genes_dir
`
