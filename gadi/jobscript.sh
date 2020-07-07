#!/bin/bash
# properties = {properties}

source /g/data1/xe2/gadi/profile.sh
source raijin/gadimod.sh

export TMPDIR=$PBS_JOBFS


set -ueo pipefail
{exec_job}
