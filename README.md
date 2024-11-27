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
## The deep learning network

![alt text](https://github.com/XiaoDongLab/deepDMI/blob/main/figures/fig1.png)

(A) The deep neural network (DNN) architecture of deepDMI. The DNN predicts the DNA methylation (DNAme) status of a CpG, referred to as the "target CpG," using three layers of input features: (i) the [-100, +100] flanking DNA sequence; (ii) DNAme status of 9 flanking CpGs within the 201-bp region; and (iii) distances of the 9 flanking CpGs from the target CpG. Feature (i) is processed through an LSTM layer followed by a fully connected layer, while features (ii) and (iii) are processed through a GRU layer, also followed by a fully connected layer. These two fully connected layers are integrated through three additional fully connected layers, culminating in a single neuron that outputs the predicted DNAme level.

(B) Epimutation detection. The DNN is trained using all CpGs of a single cell and is applied to every CpG within that cell. A significant difference between the observed DNAme status and the predicted status of a CpG indicates that its DNAme status deviates from the genome-wide pattern within the cell, and it is subsequently classified as an epimutation. 

#####
## Release Notes
• v1.0.0, 2024.11.27, Noticed capatability issue with current version of NumPy. Working on a new version.

• v1.0.0, 2022.11.22, 1st release version.
