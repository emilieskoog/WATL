#!/bin/bash
virsorter run \
--seqname-suffix-off \
--viral-gene-enrich-off \
--provirus-off \
--prep-for-dramv \
-i all_checkv_quality_filtered_typical_atypical_dereplicated_viruses_from_genomad_virsorter2.fna \
-w vs2_for_dramv \
--include-groups dsDNAphage,ssDNA,NCLDV,lavidaviridae \
--min-length 1000 \
--min-score 0 \
-j 50 \
all