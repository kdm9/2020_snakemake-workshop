#!/bin/bash
source gadi/gadimod.sh

set -ueo pipefail
logdir=gadi/log
mkdir -p $logdir
mkdir -p data/log/
export TMPDIR=${PBS_JOBFS:-$TMPDIR}
TARGET=${TARGET:-all}

QSUB="qsub -q {cluster.queue} -l ncpus={threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem} -N {cluster.name} -l storage=scratch/xe2+gdata/xe2"
QSUB="$QSUB -l wd -j oe -o $logdir -P {cluster.project}"

snakemake                                                          \
    -j 1000                                                        \
    --max-jobs-per-second 2                                        \
    --cluster-config gadi/cluster.yaml                             \
    --local-cores ${PBS_NCPUS:-1}                                  \
    --js gadi/jobscript.sh                                         \
    --nolock                                                       \
    --rerun-incomplete                                             \
    --keep-going                                                   \
    --cluster "$QSUB"                                              \
    "$TARGET"                                                      

