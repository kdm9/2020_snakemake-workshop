module purge
export TMPDIR=${PBS_JOBFS:-/tmp}
source /g/data/xe2/gadi/profile.sh
useconda
conda activate paneuc

module load seqhax samtools bcftools bwa freebayes gatk
