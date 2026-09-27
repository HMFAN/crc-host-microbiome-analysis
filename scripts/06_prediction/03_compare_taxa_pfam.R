#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(glmnet); library(digest)})

out <- "03_Results_9/RepeatedCV_mGWAS57_Pfam"
td <- file.path(out, "Tables"); qd <- file.path(out, "QA")
dir.create(td, recursive=TRUE, showWarnings=FALSE)
dir.create(qd, recursive=TRUE, showWarnings=FALSE)

inputs <- c(
  metadata="02_Analysis/BioSampleDedup_Metadata_Delivery_20260907_v1/03_covariates_NO_OLD_PCS_v2.tsv",
  abundance="02_Analysis/BioSampleDedup_Reanalysis_20260906_v1/02_Core_rSGB_v2/02_all_rSGB_MMUPHin_adjusted.tsv",
  taxonomy="03_Results_9/Figure3/Tables/Inputs/Fig3A_all190_CRC_effects.tsv",
  selected57="02_Analysis/Heri_mGWAS_879_20260909_v1/linear_results/cohort_robustness/candidate_review_v2/retained_57_traits.tsv",
  pfam_tpm="03_Results_9/Figure5/Linked_Fig5CD_Analysis/Tables_v2/Pfam38_abundance_TPM.tsv",
  pfam_map="03_Results_9/Figure5/Linked_Fig5CD_Analysis/Tables/Pfam_gene_catalog_mapping_summary.tsv"
)
stopifnot(all(file.exists(inputs)))

m <- fread(inputs["metadata"])
ab <- fread(inputs["abundance"])
tx <- fread(inputs["taxonomy"])
sel <- fread(inputs["selected57"])$trait
pm <- fread(inputs["pfam_map"])
pf <- fread(inputs["pfam_tpm"])
stopifnot(nrow(m)==1441L, uniqueN(m$BioSample)==1441L,
          sum(m$condition=="CRC")==715L, sum(m$condition=="CTR")==726L,
          uniqueN(m$Cohort)==11L, length(sel)==57L, uniqueN(sel)==57L)
m[, `:=`(CRC=as.integer(condition=="CRC"), Cohort=as.character(Cohort))]

ab <- ab[match(tx$feature, feature)]
stopifnot(!anyNA(ab$feature), all(m$IID %in% names(ab)))
xt <- log10(t(as.matrix(ab[, m$IID, with=FALSE])) + 1e-8)
colnames(xt) <- tx$trait; xt <- xt[, sel, drop=FALSE]

setnames(pf, 1L, "Pfam")
pf_ids <- sub(".quant", "", names(pf)[-1L], fixed=TRUE)
stopifnot(all(m$IID %in% pf_ids), nrow(pf)==33L, uniqueN(pf$Pfam)==33L)
xp <- t(as.matrix(pf[, -1L])); rownames(xp) <- pf_ids; colnames(xp) <- pf$Pfam
xp <- log10(xp[m$IID, , drop=FALSE] + 1)
pfam33 <- pm[n_gene_catalog_genes > 0, pfam]
pfam11_all <- pm[in_Fig5B_FDR11 == TRUE, pfam]
pfam11 <- pm[in_Fig5B_FDR11 == TRUE & n_gene_catalog_genes > 0, pfam]
pfam11_zero <- pm[in_Fig5B_FDR11 == TRUE & n_gene_catalog_genes == 0, pfam]
stopifnot(length(pfam33)==33L, length(pfam11_all)==11L, length(pfam11)==9L,
          all(pfam33 %in% colnames(xp)), all(pfam11 %in% colnames(xp)))
xp33 <- xp[, pfam33, drop=FALSE]; xp11 <- xp[, pfam11, drop=FALSE]
stopifnot(all(is.finite(xt)), all(is.finite(xp33)), all(apply(xt,2,sd)>0), all(apply(xp33,2,sd)>0))

sets <- list(mGWAS57=xt, Pfam33=xp33,
             mGWAS57_plus_Pfam33=cbind(xt,xp33),
             Pfam11=xp11, mGWAS57_plus_Pfam11=cbind(xt,xp11))

make_folds <- function(cohort,y,k,seed) {
  set.seed(seed); s <- interaction(cohort,y,drop=TRUE,lex.order=TRUE)
  f <- integer(length(y))
  for (z in levels(s)) {ii <- sample(which(s==z)); f[ii] <- rep(seq_len(k),length.out=length(ii))}
  f
}
auc_fast <- function(y,p) {r<-rank(p,ties.method="average"); (sum(r[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))}
ap_fast <- function(y,p) {o<-order(p,decreasing=TRUE); yy<-y[o]; sum(cumsum(yy)[yy==1]/which(yy==1))/sum(yy)}
sens_at_spec <- function(y,p,target) {
  thr<-c(Inf,sort(unique(p),decreasing=TRUE),-Inf)
  tp<-vapply(thr,function(t)mean(p[y==1]>=t),numeric(1)); sp<-vapply(thr,function(t)mean(p[y==0]<t),numeric(1))
  ok<-sp>=target; j<-which.max(ifelse(ok,tp,-Inf))
  c(sensitivity=tp[j],achieved_specificity=sp[j],threshold=thr[j])
}
metrics <- function(y,p) {
  a<-sens_at_spec(y,p,.90); b<-sens_at_spec(y,p,.95)
  data.table(ROC_AUC=auc_fast(y,p), PR_AUC_average_precision=ap_fast(y,p), Brier=mean((p-y)^2),
    sensitivity_at_90_specificity=a["sensitivity"], achieved_specificity_90=a["achieved_specificity"], threshold_90=a["threshold"],
    sensitivity_at_95_specificity=b["sensitivity"], achieved_specificity_95=b["achieved_specificity"], threshold_95=b["threshold"])
}

pred <- list(); met <- list(); fits <- list(); fold_columns <- list()
for (rr in 1:10) {
  outer <- make_folds(m$Cohort,m$CRC,5,20260920+rr)
  fold_columns[[rr]] <- outer
  for (ss in names(sets)) {
    x <- sets[[ss]]; pr <- rep(NA_real_,nrow(m))
    for (ff in 1:5) {
      tr <- outer!=ff; te <- !tr
      mu <- colMeans(x[tr,,drop=FALSE]); sdv <- apply(x[tr,,drop=FALSE],2,sd); sdv[sdv==0] <- 1
      ztr <- sweep(sweep(x[tr,,drop=FALSE],2,mu,"-"),2,sdv,"/")
      zte <- sweep(sweep(x[te,,drop=FALSE],2,mu,"-"),2,sdv,"/")
      inner <- make_folds(m$Cohort[tr],m$CRC[tr],5,20261920+100*rr+ff)
      cv <- cv.glmnet(ztr,m$CRC[tr],family="binomial",alpha=.5,standardize=FALSE,
                      intercept=TRUE,foldid=inner,type.measure="deviance",nlambda=60)
      pr[te] <- as.numeric(predict(cv,newx=zte,s="lambda.min",type="response"))
      fits[[length(fits)+1L]] <- data.table(model=ss,repetition=rr,outer_fold=ff,
        lambda_min=cv$lambda.min,lambda_1se=cv$lambda.1se,
        training_n=sum(tr),test_n=sum(te),nonzero=sum(as.matrix(coef(cv,s="lambda.min"))[-1,1]!=0))
    }
    stopifnot(all(is.finite(pr)))
    pred[[paste(ss,rr,sep="__")]] <- pr
    met[[length(met)+1L]] <- cbind(model=ss,repetition=rr,metrics(m$CRC,pr))
  }
  message("completed repetition ",rr,"/10")
}

met <- rbindlist(met); fits <- rbindlist(fits)
summary <- melt(met,id.vars=c("model","repetition"),variable.name="metric",value.name="value")[,
  .(mean=mean(value),sd=sd(value),min=min(value),max=max(value)),by=.(model,metric)]

delta <- dcast(met, repetition ~ model, value.var="ROC_AUC")
for (suffix in c("Pfam33","Pfam11")) delta[, paste0("delta_AUC_combined_vs_mGWAS57_",suffix) :=
  get(paste0("mGWAS57_plus_",suffix))-mGWAS57]
delta_long <- melt(delta,id.vars="repetition",measure.vars=patterns("delta_AUC"),variable.name="comparison",value.name="delta_AUC")
delta_summary <- delta_long[, {z<-delta_AUC; se<-sd(z)/sqrt(.N); q<-qt(.975,.N-1L); .(mean_delta_AUC=mean(z),sd_delta_AUC=sd(z),CI_low=mean(z)-q*se,CI_high=mean(z)+q*se,n_repeats=.N)},by=comparison]

pred_dt <- cbind(m[,.(IID,BioSample,Cohort,condition,CRC)],
  as.data.table(setNames(fold_columns,paste0("fold_repeat",1:10))),
  as.data.table(pred))
ensemble <- rbindlist(lapply(names(sets),function(ss) {
  cols <- paste(ss,1:10,sep="__"); p <- rowMeans(as.matrix(pred_dt[,..cols])); cbind(model=ss,aggregation="mean of 10 OOF probabilities",metrics(m$CRC,p))
}))

feature_manifest <- rbindlist(list(
  data.table(feature_set="mGWAS57",feature=colnames(xt),quantifiable=TRUE),
  pm[,.(feature_set="Pfam33_frozen38",feature=pfam,quantifiable=n_gene_catalog_genes>0,in_Pfam11=in_Fig5B_FDR11,n_gene_catalog_genes)],
  pm[in_Fig5B_FDR11==TRUE,.(feature_set="Pfam11",feature=pfam,quantifiable=n_gene_catalog_genes>0,in_Pfam11=TRUE,n_gene_catalog_genes)]
),fill=TRUE)

fwrite(met,file.path(td,"01_metrics_by_repeat.tsv"),sep="\t")
fwrite(summary,file.path(td,"02_metrics_summary.tsv"),sep="\t")
fwrite(delta_long,file.path(td,"03_paired_delta_AUC_by_repeat.tsv"),sep="\t")
fwrite(delta_summary,file.path(td,"04_paired_delta_AUC_summary.tsv"),sep="\t")
fwrite(ensemble,file.path(td,"05_ensemble_metrics.tsv"),sep="\t")
fwrite(pred_dt,file.path(td,"06_sample_OOF_predictions.tsv"),sep="\t")
fwrite(fits,file.path(td,"07_outer_fold_models.tsv"),sep="\t")
fwrite(feature_manifest,file.path(td,"08_fixed_feature_manifest.tsv"),sep="\t")
fwrite(data.table(role=names(inputs),path=unname(inputs),sha256=vapply(inputs,digest,"",file=TRUE,algo="sha256")),file.path(qd,"input_manifest.tsv"),sep="\t")
fwrite(data.table(N=nrow(m),CRC=sum(m$CRC),control=sum(m$CRC==0),cohorts=uniqueN(m$Cohort),
  taxa=57,pfam33=33,pfam11_frozen=11,pfam11_quantifiable=9,
  pfam11_zero_coverage=paste(pfam11_zero,collapse=";"),repeats=10,outer_folds=5,alpha=.5),file.path(qd,"run_summary.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(qd,"sessionInfo.txt"))
print(summary); print(delta_summary); print(ensemble)
