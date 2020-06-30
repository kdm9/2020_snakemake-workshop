#!/bin/bash
set -xeuo pipefail
mkdir -p outputs/3_mergebam
samtools merge  -f outputs/3_mergebam/merged.bam $(for srr in SRR097977 SRR098026; do ls "outputs/2_align/${srr}.bam"; done)
samtools index outputs/3_mergebam/merged.bam 
