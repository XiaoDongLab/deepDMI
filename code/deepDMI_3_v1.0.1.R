# R functions to call epimutations and determine accuracy based on artifical epimutations:
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

### 2. function: mse ###
call_mse = function(dat){
	observed=dat[,3]; predicted=dat[,4]
	mse=1/nrow(dat) * sum((observed - predicted)^2)
	return(mse)
}

### 3. function: determine precision & recall ###
acc = function(sn, dir_call=".", dir_acc=".", cri_diffmeth=0.5, cri_pval=0.01){
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
	
	### 7. mse ###
	dat_mse=call_mse(dat_ori)
	
	### out ###
	dat_out = c(dat_mse, nrow(dat_ori_epimut)/nrow(dat_ori), nrow(dat_ori_epimut)/nrow(dat_ori) / dat_recall * dat_precision, dat_precision, dat_recall, 2*(dat_precision*dat_recall)/(dat_precision+dat_recall))
	names(dat_out)=c("mse", "epimut_freq_raw", "epimut_freq_corrected", "precision", "recall", "f1")

	return(dat_out)
}

