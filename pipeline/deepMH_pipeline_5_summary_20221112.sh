#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
# Job name
#$ -N summary
# Number of cpu cores required
#$ -pe smp 32
# RAM requirement per cpu core
#$ -l h_vmem=4G
# require computer to handle
#$ -l h=c003
# Email
#$ -q processing.q

export LC_ALL=C
export MALLOC_ARENA_MAX=4

if [ "$#" -ne 2 ]; then
 echo "Usage: SAMPLEID" >&2
 exit 1
fi

dir_call=${1}
dir_acc=${2}
ncpu=32

date
echo "==START=="

echo ${dir_call}, ${dir_acc}, ${ncpu}

Rscript ./deepMH_analysis_5_summary_v1.0.0.R \
  -c ${dir_call} \
  -a ${dir_acc} \
  -t ${ncpu}

echo "==END=="
date
