#!/bin/bash
set -xeuo pipefail
mkdir -p outputs/2_align/
for srr in SRR097977 SRR098026
do
	# Align to reference
    bwa mem -R "@RG\tID:${srr}\tSM:${srr}" rawdata/ecoli_rel606.fa "outputs/1_qc/${srr}_qcd.fastq" | \
        samtools view -Su -o - | \
        samtools sort -o "outputs/2_align/${srr}.bam" && \
        samtools index "outputs/2_align/${srr}.bam"
done
