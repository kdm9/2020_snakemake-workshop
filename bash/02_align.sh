#!/bin/bash
set -xeuo pipefail
mkdir -p outputs/2_align/
if [ ! -f rawdata/ecoli_rel606.fa ]
then
    curl -L -o rawdata/ecoli_rel606.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
   gunzip rawdata/ecoli_rel606.fa.gz 
fi
bwa index rawdata/ecoli_rel606.fa
for srr in SRR097977 SRR098026
do
	# Align to reference
    bwa mem -R "@RG\tID:${srr}\tSM:${srr}" rawdata/ecoli_rel606.fa "outputs/1_qc/${srr}_qcd.fastq" | \
        samtools view -Su -o - | \
        samtools sort -o "outputs/2_align/${srr}.bam" && \
        samtools index "outputs/2_align/${srr}.bam"
done
