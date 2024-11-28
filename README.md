# deepDMI
A deep-learning framework for detecting DNA Methylation Instability from single-cell bisulfite sequencing data

Version: 1.0.1

Updated date: 2024.11.27

Citation: Lei Zhang, Josh Bartz, Yiwei Zhao, Linshan Laux, Shamsed Mahmud, Moonsook Lee, Xavier Revelo, Alexander Y. Maslov, Albert-László Barabási, Laura J. Niedernhofer, Paul D. Robbins, Jan Vijg and Xiao Dong. Concurrent single-cell methylomic and transcriptomic analysis shows increased epimutation burden and intra-cell type transcriptional noise in human hepatocytes during aging and in fatty liver disease. In submission.

#####
## Author and License

Author: Xiao Dong

Email: dong0265@umn.edu (X.D.)

Licensed under the GNU Affero General Public License version 3 or later.

#####
## Dependencies
• python 3.11
• tensorflow 2.18.0
• tensorflow-probability 0.25.0
• tf-keras 2.18.0

#####
## The deep neural network architecture of deepDMI

![alt text](https://github.com/XiaoDongLab/deepDMI/blob/main/figures/fig1.png)

(A) The deep neural network (DNN) architecture of deepDMI. The DNN predicts the DNA methylation (DNAme) status of a CpG, referred to as the "target CpG," using three layers of input features: (i) the [-100, +100] flanking DNA sequence; (ii) DNAme status of 9 flanking CpGs within the 201-bp region; and (iii) distances of the 9 flanking CpGs from the target CpG. Feature (i) is processed through an LSTM layer followed by a fully connected layer, while features (ii) and (iii) are processed through a GRU layer, also followed by a fully connected layer. These two fully connected layers are integrated through three additional fully connected layers, culminating in a single neuron that outputs the predicted DNAme level.

(B) Epimutation detection. The DNN is trained using all CpGs of a single cell and is applied to every CpG within that cell. A significant difference between the observed DNAme status and the predicted status of a CpG indicates that its DNAme status deviates from the genome-wide pattern within the cell, and it is subsequently classified as an epimutation. 

#####
## Usage

### Step 1. Prepare input file

• DeepDMI is supposed to take a single input file of DNAme patterns across the genome of a single cell. The input file should be formated as a tab-spaced table as the following WITHOUT its header:

| chromosome  | first position of a CpG | second position of a CpG | methylation level [0-1] |
| ----------- | ---------- | ---------- |  ---------- |
| 1 | 10484 | 10485 |  0 |
| 1 | 10489 | 10490 |  1 |
| 1 | 10493 | 10494 |  0.5 |

• Then run deepDMI_1_prepare.py to reformat the input file
```shell
version=v1.0.1; # the version of deepDMI

sn=YourSampleID; # your input file should be named as: ${sn}.cpg.pooled.bed

ref=YourReferenceGenome.fasta

# keep only autosomes
awk '$1<=22' ${sn}.cpg.pooled.bed \
  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4/($4+$5)}' \
  > ${sn}.autosome.bed

python3 deepDMI_1_prepare_${version}.py -i ${sn}.autosome -r ${ref} -b 100 -c 10
```
• The above generates a file named "${sn}.autosome.pkl", which is used as the input in step 2.


### Step 2. Train and apply deepDMI model

• In step 2, train and apply deepDMI model on the data above.
```shell
python deepDMI_2_model_${version}.py -i ${sn}.autosome -b 100
```

### Step 3. Introduce artifical epimutations to the input data

• In step 2, create 2% of artifical epimutations to the input data, and reformat the data for usage in step 4.
```shell
mkdir -p ./tmp
chmod 755 ./tmp

echo ${sn}

awk '$1<=22' ${sn}.cpg.pooled.bed \
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

python3 deepDMI_1_prepare_${version}.py -i ${sn}.autosome.artificial -r ${ref} -b 100 -c 10
```

• The above generates a file named "${sn}.autosome.artificial.pkl", which is used as the input in step 4.

### Step 4. Train and apply deepDMI model on the input with artifical epimutations

• In step 4, train and apply deepDMI model on the data above.
```shell
python3 deepDMI_2_model_${version}.py -i ${sn}.autosome.artificial -b 100
```


### Step 5. Call epimutations based on step 2 output and estimate accuracy based on step 4 output

• In step 5, utilize the R functions in the following scripts to call epimutations and estimate accuracy:
https://github.com/XiaoDongLab/deepDMI/blob/main/code/deepDMI_3_v1.0.1.R


#####
## Release Notes
• v1.0.0, 2024.11.27, Noticed capatability issue with current version of NumPy. Working on a new version.

• v1.0.0, 2022.11.22, 1st release version.
