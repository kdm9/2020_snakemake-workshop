#!/bin/bash
set -xeuo pipefail
mkdir -p outputs/3_varcall
bcftools mpileup -f rawdata/ecoli_rel606.fa $(for srr in SRR097977 SRR098026; do ls "outputs/2_align/${srr}.bam"; done) | \
    bcftools call --ploidy 1 -c | \
    bcftools view -e "QUAL<30 || ALT != '.'" -o outputs/3_varcall/variants.vcf.gz -Oz

