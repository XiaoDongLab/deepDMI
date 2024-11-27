#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
# Job name
#$ -N gpu1_epi
# Number of cpu cores required
#$ -pe smp 1
# RAM requirement per cpu core
#$ -l h_vmem=80G
# require computer to handle
# -l h=x002
#$ -q gpu1.q

export LC_ALL=C
export MALLOC_ARENA_MAX=4
source ~/apps/python-3.8-tensorflow-venv/bin/activate
export CUDA_VISIBLE_DEVICES=0; ### use 0 for gpu queue 1; and 1 for gpu queue 2

if [ "$#" -ne 1 ]; then
 echo "Usage: SAMPLEID" >&2
 exit 1
fi

sn=${1}

date
echo "==START=="

hostname

python deepMH_analysis_2_modeling_v1.0.0.py -i ${1}.autosome -b 100

echo "==END=="
date
