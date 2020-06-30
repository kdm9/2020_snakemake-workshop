#!/bin/bash
set -xeuo pipefail
mkdir -p outputs/1_qc/
for srr in SRR097977 SRR098026
do
	# run scythe, an adaptor trimmer, on each raw read file
	scythe \
		-a "rawdata/illumina_adapters.fa" \
		-q sanger \
		-o "outputs/1_qc/${srr}_qcd.fastq" \
		"rawdata/reads/${srr}.fastq"
done
