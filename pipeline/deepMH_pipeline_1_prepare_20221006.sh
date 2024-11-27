#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
# Job name
#$ -N pre
# Number of cpu cores required
#$ -pe smp 24
# RAM requirement per cpu core
#$ -l h_vmem=20G
# require computer to handle
# -l h=x001
# Email
#$ -q processing.q

export LC_ALL=C
export MALLOC_ARENA_MAX=4
source ~/apps/python-3.8-tensorflow-venv/bin/activate

if [ "$#" -ne 1 ]; then
 echo "Usage: SAMPLEID" >&2
 exit 1
fi

sn=${1}
ref=~/references/homo_sapiens/GRCh38.p13_BS/Homo_sapiens.GRCh38.dna.primary_assembly.fa

date
echo "==START=="

echo ${sn}

awk '$1<=22' ../alignment/4-2-pool/${sn}.cpg.pooled.bed \
  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4/($4+$5)}' \
  > ${sn}.autosome.bed

python3 deepMH_analysis_1_prepare_v1.0.0.py -i ${sn}.autosome -r ${ref} -b 100 -c 10

rm ${sn}.autosome.bed

echo "==END=="
date

