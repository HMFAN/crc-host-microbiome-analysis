#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(glmnet);library(digest)})
set.seed(20260920)

out <- "03_Results_9/Fig6b_RepeatedCV_SHAP_mGWAS57"
td <- file.path(out,"Tables"); qd <- file.path(out,"QA")
dir.create(td,recursive=TRUE,showWarnings=FALSE);dir.create(qd,recursive=TRUE,showWarnings=FALSE)
inputs <- c(
 metadata="02_Analysis/BioSampleDedup_Metadata_Delivery_20260907_v1/03_covariates_NO_OLD_PCS_v2.tsv",
 abundance="02_Analysis/BioSampleDedup_Reanalysis_20260906_v1/02_Core_rSGB_v2/02_all_rSGB_MMUPHin_adjusted.tsv",
 taxonomy="03_Results_9/Figure3/Tables/Inputs/Fig3A_all190_CRC_effects.tsv",
 selected57="02_Analysis/Heri_mGWAS_879_20260909_v1/linear_results/cohort_robustness/candidate_review_v2/retained_57_traits.tsv",
 fig6a_predictions="03_Results_9/RepeatedCV_mGWAS57/Tables/04_sample_OOF_predictions.tsv",
 fig6a_models="03_Results_9/RepeatedCV_mGWAS57/Tables/05_outer_fold_models.tsv"
)
stopifnot(all(file.exists(inputs)))
m <- fread(inputs["metadata"]); ab <- fread(inputs["abundance"]); tx <- fread(inputs["taxonomy"])
sel <- fread(inputs["selected57"])$trait; prior <- fread(inputs["fig6a_predictions"]); prior_models <- fread(inputs["fig6a_models"])
stopifnot(nrow(m)==1441,uniqueN(m$BioSample)==1441,sum(m$condition=="CRC")==715,sum(m$condition=="CTR")==726,
 uniqueN(m$Cohort)==11,length(sel)==57,uniqueN(sel)==57,all(sel%in%tx$trait),nrow(prior_models)==50)
m[,`:=`(CRC=as.integer(condition=="CRC"),Cohort=as.character(Cohort))]
stopifnot(identical(m$IID,prior$IID),identical(m$BioSample,prior$BioSample))
ab <- ab[match(tx$feature,feature)];stopifnot(!anyNA(ab$feature),all(m$IID%in%names(ab)))
x <- log10(t(as.matrix(ab[,m$IID,with=FALSE]))+1e-8);colnames(x)<-tx$trait;x<-x[,sel,drop=FALSE]
stopifnot(all(is.finite(x)),all(apply(x,2,sd)>0))

make_folds <- function(cohort,y,k,seed){set.seed(seed);s<-interaction(cohort,y,drop=TRUE,lex.order=TRUE);f<-integer(length(y));for(z in levels(s)){ii<-sample(which(s==z));f[ii]<-rep(seq_len(k),length.out=length(ii))};f}

n <- nrow(x); p <- ncol(x); R <- 10L
shap_sum <- matrix(0,n,p,dimnames=list(m$IID,sel)); z_sum <- matrix(0,n,p,dimnames=list(m$IID,sel))
model_rows <- vector("list",R*5L); coef_rows <- vector("list",R*5L); add_err <- numeric(R*5L); prob_err <- numeric(R*5L)
kk <- 0L
for(rep_id in 1:R){
 outer <- make_folds(m$Cohort,m$CRC,5,20260920+rep_id)
 stopifnot(identical(as.integer(outer),as.integer(prior[[paste0("fold_repeat",rep_id)]])))
 for(ff in 1:5){
  kk <- kk+1L; tr <- outer!=ff; te <- !tr
  mu <- colMeans(x[tr,,drop=FALSE]); ss <- apply(x[tr,,drop=FALSE],2,sd); ss[ss==0]<-1
  ztr <- sweep(sweep(x[tr,,drop=FALSE],2,mu,"-"),2,ss,"/")
  zte <- sweep(sweep(x[te,,drop=FALSE],2,mu,"-"),2,ss,"/")
  inner <- make_folds(m$Cohort[tr],m$CRC[tr],5,20261920+100*rep_id+ff)
  cv <- cv.glmnet(ztr,m$CRC[tr],family="binomial",alpha=.5,standardize=FALSE,intercept=TRUE,
    foldid=inner,type.measure="deviance",nlambda=60)
  cf <- as.matrix(coef(cv,s="lambda.min"))[,1]; beta <- cf[sel]; base <- unname(cf["(Intercept)"])
  sv <- sweep(zte,2,beta,"*")
  eta <- base+rowSums(sv); phat <- plogis(eta)
  shap_sum[te,] <- shap_sum[te,,drop=FALSE]+sv; z_sum[te,] <- z_sum[te,,drop=FALSE]+zte
  refp <- prior[[paste0("probability_repeat",rep_id)]][te]
  prob_err[kk] <- max(abs(phat-refp)); add_err[kk] <- max(abs(eta-as.numeric(predict(cv,newx=zte,s="lambda.min",type="link"))))
  model_rows[[kk]] <- data.table(repetition=rep_id,outer_fold=ff,training_n=sum(tr),test_n=sum(te),lambda_min=cv$lambda.min,
    lambda_1se=cv$lambda.1se,nonzero_taxa=sum(beta!=0),intercept=base,max_probability_difference_vs_Fig6a=prob_err[kk],max_SHAP_additivity_error=add_err[kk])
  coef_rows[[kk]] <- data.table(repetition=rep_id,outer_fold=ff,trait=sel,coefficient=beta,training_mean=mu[sel],training_sd=ss[sel])
 }
}
shap <- shap_sum/R; zavg <- z_sum/R
stopifnot(all(is.finite(shap)),all(is.finite(zavg)),max(prob_err)<1e-12,max(add_err)<1e-10)
modeldt <- rbindlist(model_rows); coefdt <- rbindlist(coef_rows)

lut <- unique(tx[trait%in%sel,.(trait,feature,species,phylum,class,order,genus,beta_Fig3=beta)])
stopifnot(nrow(lut)==57,all(sel%in%lut$trait))
overall <- data.table(trait=sel,mean_abs_SHAP=colMeans(abs(shap)),mean_signed_SHAP=colMeans(shap),
 mean_SHAP_CRC=colMeans(shap[m$CRC==1,,drop=FALSE]),mean_SHAP_control=colMeans(shap[m$CRC==0,,drop=FALSE]),
 abundance_SHAP_Spearman=vapply(seq_len(p),function(j)cor(zavg[,j],shap[,j],method="spearman"),numeric(1)))
overall[,`:=`(CRC_minus_control=mean_SHAP_CRC-mean_SHAP_control,
 abundance_SHAP_direction=fifelse(abundance_SHAP_Spearman>0,"higher abundance -> more positive CRC log-odds SHAP",
  fifelse(abundance_SHAP_Spearman<0,"higher abundance -> more negative CRC log-odds SHAP","no monotonic direction")))]
overall <- merge(overall,lut,by="trait",all.x=TRUE,sort=FALSE)
setorder(overall,-mean_abs_SHAP,trait);overall[,rank_mean_abs_SHAP:=seq_len(.N)]
overall[,lineage:=fifelse(order%chin%c("Bacteroidales","Oscillospirales","Lachnospirales"),order,"Other")]
top20 <- overall[1:20]
bc <- overall[species=="Bacteroides_caccae" | trait=="core161"]
stopifnot(nrow(bc)==1L)

samplemeta <- m[,.(IID,BioSample,Cohort,condition,CRC)]
shap_wide <- cbind(samplemeta,as.data.table(shap)); z_wide <- cbind(samplemeta,as.data.table(zavg))
plotdt <- rbindlist(lapply(top20$trait,function(v){j<-match(v,sel); z<-lut[trait==v]
 data.table(IID=m$IID,BioSample=m$BioSample,Cohort=m$Cohort,condition=m$condition,CRC=m$CRC,trait=v,
  species=z$species,order=z$order,lineage=fifelse(z$order%chin%c("Bacteroidales","Oscillospirales","Lachnospirales"),z$order,"Other"),
  mean_OOF_SHAP=shap[,j],mean_OOF_standardized_abundance=zavg[,j],rank_mean_abs_SHAP=top20[trait==v,rank_mean_abs_SHAP])
}))

fwrite(overall,file.path(td,"01_taxon_SHAP_summary.tsv"),sep="\t")
fwrite(top20,file.path(td,"02_top20_taxa_lineage.tsv"),sep="\t")
fwrite(bc,file.path(td,"03_Bacteroides_caccae_SHAP.tsv"),sep="\t")
fwrite(shap_wide,file.path(td,"04_sample_mean_OOF_SHAP_matrix.tsv"),sep="\t")
fwrite(z_wide,file.path(td,"05_sample_mean_OOF_standardized_abundance.tsv"),sep="\t")
fwrite(plotdt,file.path(td,"06_Fig6b_plot_data.tsv"),sep="\t")
fwrite(modeldt,file.path(td,"07_restored_Fig6a_models.tsv"),sep="\t")
fwrite(coefdt,file.path(td,"08_restored_model_coefficients_scaling.tsv"),sep="\t")
fwrite(data.table(role=names(inputs),path=unname(inputs),sha256=vapply(inputs,digest,"",file=TRUE,algo="sha256")),file.path(qd,"input_manifest.tsv"),sep="\t")
fwrite(data.table(N=n,CRC=sum(m$CRC),control=sum(!m$CRC),cohorts=uniqueN(m$Cohort),taxa=p,repeats=R,outer_models=nrow(modeldt),
 SHAP_scale="CRC log-odds",SHAP_method="exact interventional linear SHAP with outer-training background mean",
 max_probability_difference_vs_Fig6a=max(prob_err),max_SHAP_additivity_error=max(add_err),Bcaccae_rank=bc$rank_mean_abs_SHAP),file.path(qd,"run_summary.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(qd,"sessionInfo.txt"))
print(overall[,.(rank_mean_abs_SHAP,trait,species,lineage,mean_abs_SHAP,mean_SHAP_CRC,mean_SHAP_control,CRC_minus_control,abundance_SHAP_Spearman)][1:20])
print(bc)
