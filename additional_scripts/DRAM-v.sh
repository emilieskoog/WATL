DRAM-v.py annotate \
  -i vs2_for_dramv/for-dramv/final-viral-combined-for-dramv.fa \
  -v vs2_for_dramv/for-dramv/viral-affi-contigs-for-dramv.tab \
  -o dramv_annotate --threads 20
DRAM-v.py distill -i dramv_annotate/annotations.tsv -o dramv_distill