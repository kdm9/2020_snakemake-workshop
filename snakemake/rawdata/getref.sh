curl -L -o ecoli_rel606.fa.gz http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip ecoli_rel606.fa.gz 
bwa index ecoli_rel606.fa
