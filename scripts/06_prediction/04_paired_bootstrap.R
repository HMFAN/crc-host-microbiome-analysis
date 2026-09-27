#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(digest)})

root <- "03_Results_9/RepeatedCV_mGWAS57_Pfam"
td <- file.path(root,"Tables"); qd <- file.path(root,"QA")
dir.create(td,recursive=TRUE,showWarnings=FALSE); dir.create(qd,recursive=TRUE,showWarnings=FALSE)
infile <- file.path(td,"06_sample_OOF_predictions.tsv")
stopifnot(file.exists(infile))
d <- fread(infile)
base_cols <- paste("mGWAS57",1:10,sep="__")
comb_cols <- paste("mGWAS57_plus_Pfam33",1:10,sep="__")
stopifnot(all(c("IID","Cohort","CRC",base_cols,comb_cols) %in% names(d)),nrow(d)==1441L)
d[, p_mGWAS57 := rowMeans(as.matrix(.SD)), .SDcols=base_cols]
d[, p_combined := rowMeans(as.matrix(.SD)), .SDcols=comb_cols]
stopifnot(all(is.finite(d$p_mGWAS57)),all(is.finite(d$p_combined)))

auc_fast <- function(y,p) {
  r <- rank(p,ties.method="average")
  (sum(r[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))
}
ap_fast <- function(y,p) {
  o <- order(p,decreasing=TRUE); yy <- y[o]
  sum(cumsum(yy)[yy==1]/which(yy==1))/sum(yy)
}
stats <- function(ii) {
  y <- d$CRC[ii]; p0 <- d$p_mGWAS57[ii]; p1 <- d$p_combined[ii]
  c(delta_AUC=auc_fast(y,p1)-auc_fast(y,p0),
    delta_PR_AUC=ap_fast(y,p1)-ap_fast(y,p0),
    delta_Brier=mean((p1-y)^2)-mean((p0-y)^2))
}

strata <- split(seq_len(nrow(d)),interaction(d$Cohort,d$CRC,drop=TRUE,lex.order=TRUE))
B <- as.integer(Sys.getenv("N_BOOT","10000")); seed <- 20260921L
set.seed(seed)
observed <- stats(seq_len(nrow(d)))
boot <- matrix(NA_real_,nrow=B,ncol=length(observed),dimnames=list(NULL,names(observed)))
for (b in seq_len(B)) {
  ii <- unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE)
  boot[b,] <- stats(ii)
  if (b %% 1000L == 0L) message("completed bootstrap ",b,"/",B)
}

summ <- rbindlist(lapply(seq_along(observed),function(j) {
  z <- boot[,j]; lo <- (1+sum(z<=0))/(B+1); hi <- (1+sum(z>=0))/(B+1)
  data.table(metric=names(observed)[j],definition=c(
    "AUC(combined) - AUC(mGWAS57)",
    "PR-AUC(combined) - PR-AUC(mGWAS57)",
    "Brier(combined) - Brier(mGWAS57); negative favours combined")[j],
    observed_delta=observed[j],CI_low=quantile(z,.025,names=FALSE),CI_high=quantile(z,.975,names=FALSE),
    two_sided_P=min(1,2*min(lo,hi)),bootstrap_mean=mean(z),bootstrap_sd=sd(z),B=B,seed=seed)
}))
boot_dt <- as.data.table(boot); boot_dt[,bootstrap:=seq_len(.N)]; setcolorder(boot_dt,"bootstrap")
sample_probs <- d[,.(IID,Cohort,CRC,p_mGWAS57,p_combined)]

fwrite(summ,file.path(td,"09_paired_stratified_bootstrap_summary.tsv"),sep="\t")
fwrite(boot_dt,file.path(td,"10_paired_stratified_bootstrap_replicates.tsv.gz"),sep="\t")
fwrite(sample_probs,file.path(td,"11_ensemble_probabilities_for_bootstrap.tsv"),sep="\t")
fwrite(data.table(input=infile,sha256=digest(infile,file=TRUE,algo="sha256"),N=nrow(d),
  strata=length(strata),bootstrap_replicates=B,seed=seed,resampling="paired within cohort x CRC strata",
  CI="percentile",P="two-sided bootstrap sign P with plus-one correction"),
  file.path(qd,"paired_bootstrap_QA.tsv"),sep="\t")
print(summ)
