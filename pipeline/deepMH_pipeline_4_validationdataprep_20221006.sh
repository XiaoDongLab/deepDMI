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
# -q all.q

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

mkdir -p ./tmp
chmod 755 ./tmp

echo ${sn}

awk '$1<=22' ../alignment/4-2-pool/${sn}.cpg.pooled.bed \
  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4/($4+$5)}' \
  | sort -R -T ./tmp > ${sn}.autosome.tmp1.bed

len=$(wc -l ${sn}.autosome.tmp1.bed | awk '{print $1}')

# pick only 2% data as artifical problems
len1=$(expr $len / 50)
len2=$(expr $len - $len1)

head -n ${len1} ${sn}.autosome.tmp1.bed \
  | awk '$4 > 0.75' | awk '{print $1 "\t" $2 "\t" $3 "\t" 0}' \
  > ${sn}.autosome.tmp2.bed

head -n ${len1} ${sn}.autosome.tmp1.bed \
  | awk '$4 < 0.25' | awk '{print $1 "\t" $2 "\t" $3 "\t" 1}' \
  >> ${sn}.autosome.tmp2.bed

sort -T ./tmp -k1,1d -k2,2n ${sn}.autosome.tmp2.bed > ${sn}.autosome.artificialonly.bed

head -n ${len1} ${sn}.autosome.tmp1.bed \
  | awk '$4 >= 0.25 && $4 <= 0.75' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' \
  >> ${sn}.autosome.tmp2.bed

tail -n ${len2} ${sn}.autosome.tmp1.bed >> ${sn}.autosome.tmp2.bed

sort -T ./tmp -k1,1d -k2,2n ${sn}.autosome.tmp2.bed > ${sn}.autosome.artificial.bed

rm ${sn}.autosome.tmp?.bed

python3 deepMH_analysis_1_prepare_v1.0.0.py -i ${sn}.autosome.artificial -r ${ref} -b 100 -c 10

echo "==END=="
date

