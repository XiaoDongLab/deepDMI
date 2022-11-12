
homedir="/data/gv_h/lzhang/projects/2022-scepi/2021-scepi-mm"
# homedir="/data/gv_h/lzhang/projects/2022-scepi/2022-liver-batch1/bs"
# homedir="/data/gv_h/lzhang/projects/2022-scepi/2022-liver-batch2/bs"
numCores=60

# dir_acc=paste(homedir, "/epimut_call_b100c10_accuracy", sep="")
dir_acc=paste(homedir, "/epimut_call_accuracy", sep="")
dir_call=paste(homedir, "/epimut_call", sep="")

require(parallel)
require(MASS)
require(ggplot2)

setwd(dir_acc)

sns_dat <- read.table("../samplelist.txt", header=F, sep="\t")
colnames(sns_dat) <- c("sampleid", "subjectid", "type", "age", "sex", "class")

sns = as.vector(sns_dat$sampleid)

### 1. function: calling epimutation from a dataset ###
# call_epi = function(dat, cri_diffmeth=0.5, cri_pval=0.01){
call_epi = function(dat){
	observed=dat[,3]; predicted=dat[,4]
	mse=1/nrow(dat) * sum((observed - predicted)^2)

	# Different to microarray, due to missing data, site specific mse may not known, so, just use the genome-wide sample-specific mse to estimate sd.
	cdf1=pnorm(q=observed, mean=predicted,sd=mse^(1/2))
	cdf2=1-cdf1
	cdf=cbind(cdf1, cdf2)
	pval=apply(cdf, 1, min)
	pval.adjust = p.adjust(pval, method = "BH")

	# methdiff, defined as predicted - observed; if <0, it means unexpected hypermeth; if >0, unexpected hypometh
	methdiff = predicted - observed
	dat <- cbind(dat, methdiff, pval)
	colnames(dat) = c("chr", "pos", "observed_meth", "predicted_meth", "methdiff", "pval")
	return(dat)
}

### 2. function: determine precision & recall ###
acc = function(sn, cri_diffmeth=0.5, cri_pval=0.01){
	# default criteria:
	# 1. diff_meth > 0.5
	# 2. p < 0.01
	### 2. epimutations from original data ###
	dat_ori=read.table(paste(dir_call, "/", sn, ".autosome.predict2", sep=""), header=F)
	dat_ori_all = call_epi(dat_ori)
	write.table(dat_ori_all, paste(dir_call, "/", sn, ".autosome.results.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
	dat_ori_epimut = dat_ori_all[abs(dat_ori_all$methdiff)>cri_diffmeth & dat_ori_all$pval<cri_pval, ]

	### 3. epimutations from artifical data ###
	dat_acc=read.table(paste(dir_acc, "/", sn, ".autosome.artificial.predict2", sep=""), header=F)
	dat_acc_all = call_epi(dat_acc)
	dat_acc_epimut = dat_acc_all[abs(dat_acc_all$methdiff)>cri_diffmeth & dat_acc_all$pval<cri_pval, ]

	### 4. load artifically generated epimutations ###
	dat_art=read.table(paste(dir_acc, "/", sn, ".autosome.artificialonly.bed", sep=""), header=F)

	### 5. sensitivity ### 
	# 5.1. find dat_art in dat_ori, but not in dat_ori_epimut
	# dat_ori was a subset of all input because of my filter; dat_art was created before the filter
	tmp = merge(dat_art[,c(1,3)], dat_ori[, 1:2], by=1:2, sort=F)
	tmp = merge(tmp, dat_ori_epimut, by=1:2, all=T, sort=F)
	tmp = tmp[is.na(tmp[, 3]), ]
	# 5.2. check if the above artifical ones are observed in acc data
	tmp1 = merge(dat_acc_epimut, tmp[, c(1,2)], by=1:2)
	dat_recall = nrow(tmp1)/nrow(tmp)

	### 6. precision ### 
	# 6.1. FPs: find those in dat_acc_epimut, but not in either dat_ori_epimut or dat_art
	tmp2=merge(dat_ori_epimut, dat_art[, c(1,3)], by=1:2, all=T)[, 1:2]
	tmp3=merge(dat_acc_epimut[, 1:2], cbind(tmp2, "all"), by=1:2, all.x=T, sort=F)
	tmp4=tmp3[is.na(tmp3[,3]), ]
	# 6.2. precision
	dat_precision = nrow(tmp1)/(nrow(tmp1) + nrow(tmp4))
	dat_out = c(nrow(dat_ori_epimut)/nrow(dat_ori), nrow(dat_ori_epimut)/nrow(dat_ori) / dat_recall * dat_precision, dat_precision, dat_recall, 2*(dat_precision*dat_recall)/(dat_precision+dat_recall))
	names(dat_out)=c("epimut_freq_raw", "epimut_freq_corrected", "precision", "recall", "f1")

	return(dat_out)
}

calculate = function(i){
	print(sns[i])
	dat_summary=vector()
	x1 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.05), 4)
	x2 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.01), 5)
	x3 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.001), 6)
	x4 = c(sns[i], acc(sns[i], cri_diffmeth=0.75, cri_pval=0.001), 7)
	x5 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.0005), 8)
	x6 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.0001), 9)
	x7 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.00001), 10)
	x8 = c(sns[i], acc(sns[i], cri_diffmeth=0.75, cri_pval=0.00001), 11)
	x9 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.000001), 12)
	x10 = c(sns[i], acc(sns[i], cri_diffmeth=0.5, cri_pval=0.0000001), 13)
	dat_summary=rbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
	return(dat_summary)
}

result=vector();
result <- mclapply(1:length(sns), calculate, mc.cores = numCores)

tmp=vector()
for(i in 1:length(result)){
	tmp = rbind(tmp, result[[i]])
}

colnames(tmp) = c('sampleid', 'epimut_freq_raw', "epimut_freq_corrected", 'precision', 'recall', 'f1', 'cri_no')

write.table(tmp, paste(dir_call, "/summary_epimutfreq.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

tmp <- read.table(paste(dir_call, "/summary_epimutfreq.txt", sep=""), header=T)
for(i in 4:13){
	tmp1=tmp[tmp$cri_no==i, ]
	tmp1 = merge(tmp1, sns_dat, by=1)
	tmp1$subjectid = as.factor(tmp1$subjectid)
	tmp1$type = as.factor(tmp1$type)
	tmp1$age= as.factor(tmp1$age)
	tmp1$class= as.factor(tmp1$class)
	tmp2 = tmp1[tmp1$type=="Cell", ]; # limiting only to single cells, excluding bulk

	pdf(paste("summary_accuracy_precision_cri", i,".pdf", sep=''))
	p <- ggplot(tmp2, aes(x=subjectid, y=precision)) + 
	        ylim(0, 1) + 
	        geom_boxplot(outlier.shape = NA) +
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			theme_classic()
	print(p)
	dev.off()

	pdf(paste("summary_accuracy_recall_cri", i,".pdf", sep=''))
	p <- ggplot(tmp2, aes(x=subjectid, y=recall)) + 
	        ylim(0, 1) + 
	        geom_boxplot(outlier.shape = NA) + 
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			theme_classic()
	print(p)
	dev.off()

	pdf(paste("summary_accuracy_f1_cri", i,".pdf", sep=''))
	p <- ggplot(tmp2, aes(x=subjectid, y=f1)) + 
	        ylim(0, 1) + 
	        geom_boxplot(outlier.shape = NA) +
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			theme_classic()
	print(p)
	dev.off()

	pdf(paste("epimut_freq_cri", i,".pdf", sep=''))
	p <- ggplot(tmp2, aes(x=subjectid, y=epimut_freq_corrected)) + 
	        ylim(0, ceiling(max(tmp2$epimut_freq_corrected)*100) / 100) + 
	        geom_boxplot(outlier.shape = NA) + 
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			theme_classic()
	print(p)
	dev.off()

	pdf(paste("epimut_freqraw_cri", i,".pdf", sep=''))
	p <- ggplot(tmp2, aes(x=subjectid, y=epimut_freq_raw)) + 
	        ylim(0, ceiling(max(tmp2$epimut_freq_corrected)*100) / 100) + 
	        geom_boxplot(outlier.shape = NA) + 
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			theme_classic()
	print(p)
	dev.off()
}

