#!/bin/bash
set -xeuo pipefail
mkdir -p outputs/4_varcall
bcftools mpileup -f rawdata/ecoli_rel606.fa outputs/3_mergebam/merged.bam | \
    bcftools call --ploidy 1 -c | \
    bcftools view -e "QUAL<30 || ALT != '.'" -o outputs/4_varcall/variants.vcf.gz -Oz

