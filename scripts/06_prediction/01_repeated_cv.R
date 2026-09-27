#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(glmnet);library(parallel);library(digest)})
set.seed(20260920)
out<-"03_Results_9/RepeatedCV_mGWAS57";td<-file.path(out,"Tables");qd<-file.path(out,"QA")
dir.create(td,recursive=TRUE,showWarnings=FALSE);dir.create(qd,recursive=TRUE,showWarnings=FALSE)
inputs<-c(metadata="02_Analysis/BioSampleDedup_Metadata_Delivery_20260907_v1/03_covariates_NO_OLD_PCS_v2.tsv",
 abundance="02_Analysis/BioSampleDedup_Reanalysis_20260906_v1/02_Core_rSGB_v2/02_all_rSGB_MMUPHin_adjusted.tsv",
 taxonomy="03_Results_9/Figure3/Tables/Inputs/Fig3A_all190_CRC_effects.tsv",
 selected57="02_Analysis/Heri_mGWAS_879_20260909_v1/linear_results/cohort_robustness/candidate_review_v2/retained_57_traits.tsv")
stopifnot(all(file.exists(inputs)))
m<-fread(inputs["metadata"]);ab<-fread(inputs["abundance"]);tx<-fread(inputs["taxonomy"]);sel<-fread(inputs["selected57"])$trait
stopifnot(nrow(m)==1441,uniqueN(m$BioSample)==1441,sum(m$condition=="CRC")==715,sum(m$condition=="CTR")==726,
 uniqueN(m$Cohort)==11,length(sel)==57,uniqueN(sel)==57,all(sel%in%tx$trait))
m[,`:=`(CRC=as.integer(condition=="CRC"),Cohort=as.character(Cohort))]
ab<-ab[match(tx$feature,feature)];stopifnot(!anyNA(ab$feature),all(m$IID%in%names(ab)))
x<-log10(t(as.matrix(ab[,m$IID,with=FALSE]))+1e-8);colnames(x)<-tx$trait;x<-x[,sel,drop=FALSE]
stopifnot(all(is.finite(x)),all(apply(x,2,sd)>0))

make_folds<-function(cohort,y,k,seed){set.seed(seed);s<-interaction(cohort,y,drop=TRUE,lex.order=TRUE);f<-integer(length(y));for(z in levels(s)){ii<-sample(which(s==z));f[ii]<-rep(seq_len(k),length.out=length(ii))};f}
auc_fast<-function(y,p){r<-rank(p,ties.method="average");(sum(r[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0))}
ap_fast<-function(y,p){o<-order(p,decreasing=TRUE);yy<-y[o];sum(cumsum(yy)[yy==1]/which(yy==1))/sum(yy)}
sens_at_spec<-function(y,p,target){thr<-c(Inf,sort(unique(p),decreasing=TRUE),-Inf);tp<-vapply(thr,function(t)mean(p[y==1]>=t),numeric(1));sp<-vapply(thr,function(t)mean(p[y==0]<t),numeric(1));ok<-sp>=target;j<-which.max(ifelse(ok,tp,-Inf));c(sensitivity=tp[j],achieved_specificity=sp[j],threshold=thr[j])}
metrics<-function(y,p){a<-sens_at_spec(y,p,.90);b<-sens_at_spec(y,p,.95);data.table(ROC_AUC=auc_fast(y,p),PR_AUC_average_precision=ap_fast(y,p),Brier=mean((p-y)^2),
 sensitivity_at_90_specificity=a["sensitivity"],achieved_specificity_90=a["achieved_specificity"],threshold_90=a["threshold"],
 sensitivity_at_95_specificity=b["sensitivity"],achieved_specificity_95=b["achieved_specificity"],threshold_95=b["threshold"])}

one_repeat<-function(y,rep_id,save_detail=FALSE){
 outer<-make_folds(m$Cohort,y,5,20260920+rep_id);pr<-rep(NA_real_,nrow(x));mods<-vector("list",5)
 for(ff in 1:5){tr<-outer!=ff;te<-!tr;mu<-colMeans(x[tr,,drop=FALSE]);ss<-apply(x[tr,,drop=FALSE],2,sd);ss[ss==0]<-1
  ztr<-sweep(sweep(x[tr,,drop=FALSE],2,mu,"-"),2,ss,"/");zte<-sweep(sweep(x[te,,drop=FALSE],2,mu,"-"),2,ss,"/")
  inner<-make_folds(m$Cohort[tr],y[tr],5,20261920+100*rep_id+ff)
  cv<-cv.glmnet(ztr,y[tr],family="binomial",alpha=.5,standardize=FALSE,intercept=TRUE,foldid=inner,type.measure="deviance",nlambda=60)
  pr[te]<-as.numeric(predict(cv,newx=zte,s="lambda.min",type="response"))
  if(save_detail)mods[[ff]]<-data.table(repetition=rep_id,outer_fold=ff,training_n=sum(tr),test_n=sum(te),training_crc=sum(y[tr]),test_crc=sum(y[te]),lambda_min=cv$lambda.min,lambda_1se=cv$lambda.1se,nonzero_taxa=sum(as.matrix(coef(cv,s="lambda.min"))[-1,1]!=0))
 }
 stopifnot(all(is.finite(pr)));list(pred=pr,fold=outer,metric=cbind(repetition=rep_id,metrics(y,pr)),models=if(save_detail)rbindlist(mods) else NULL)
}

obs<-lapply(1:10,function(r)one_repeat(m$CRC,r,TRUE));predmat<-do.call(cbind,lapply(obs,`[[`,"pred"));colnames(predmat)<-paste0("probability_repeat",1:10)
foldmat<-do.call(cbind,lapply(obs,`[[`,"fold"));colnames(foldmat)<-paste0("fold_repeat",1:10)
repmet<-rbindlist(lapply(obs,`[[`,"metric"));modeldt<-rbindlist(lapply(obs,`[[`,"models"))
avgpred<-rowMeans(predmat);ensemble<-cbind(aggregation="mean of 10 OOF probabilities",metrics(m$CRC,avgpred))
summarydt<-melt(repmet,id.vars="repetition",variable.name="metric",value.name="value")[,.(mean=mean(value),sd=sd(value),min=min(value),max=max(value)),by=metric]

B<-as.integer(Sys.getenv("N_PERM","200"));requested_workers<-as.integer(Sys.getenv("N_WORKERS","8"));cores<-suppressWarnings(detectCores());if(is.na(cores))cores<-requested_workers;workers<-min(requested_workers,cores,B)
perm_one<-function(bb){set.seed(20270000+bb);yp<-unlist(lapply(split(seq_len(nrow(m)),m$Cohort),function(ii)sample(m$CRC[ii])),use.names=FALSE)
 ii_order<-unlist(split(seq_len(nrow(m)),m$Cohort),use.names=FALSE);yy<-integer(nrow(m));yy[ii_order]<-yp
 aa<-vapply(1:10,function(r)one_repeat(yy,r,FALSE)$metric$ROC_AUC,numeric(1));c(permutation=bb,mean_repeat_AUC=mean(aa),sd_repeat_AUC=sd(aa))}
if(.Platform$OS.type=="unix" && workers>1){perm<-mclapply(seq_len(B),perm_one,mc.cores=workers,mc.set.seed=FALSE)}else perm<-lapply(seq_len(B),perm_one)
permdt<-rbindlist(lapply(perm,as.list));observed_stat<-mean(repmet$ROC_AUC);ge<-sum(permdt$mean_repeat_AUC>=observed_stat)
ptest<-data.table(statistic="mean ROC-AUC across 10 repeated OOF runs",observed=observed_stat,n_permutations=B,n_permuted_ge_observed=ge,
 empirical_one_sided_P=(1+ge)/(B+1),null_mean=mean(permdt$mean_repeat_AUC),null_sd=sd(permdt$mean_repeat_AUC),
 null_median=median(permdt$mean_repeat_AUC),null_95_CI_low=quantile(permdt$mean_repeat_AUC,.025,names=FALSE),
 null_95_CI_high=quantile(permdt$mean_repeat_AUC,.975,names=FALSE),
 observed_percentile_in_null=100*mean(permdt$mean_repeat_AUC<=observed_stat),
 null_min=min(permdt$mean_repeat_AUC),null_max=max(permdt$mean_repeat_AUC))

predout<-cbind(m[,.(IID,BioSample,Cohort,condition,CRC)],as.data.table(foldmat),as.data.table(predmat));predout[,mean_OOF_probability:=avgpred]
fwrite(repmet,file.path(td,"01_metrics_by_repeat.tsv"),sep="\t");fwrite(summarydt,file.path(td,"02_metrics_summary.tsv"),sep="\t")
fwrite(ensemble,file.path(td,"03_mean_OOF_probability_metrics.tsv"),sep="\t");fwrite(predout,file.path(td,"04_sample_OOF_predictions.tsv"),sep="\t")
fwrite(modeldt,file.path(td,"05_outer_fold_models.tsv"),sep="\t");fwrite(permdt,file.path(td,"06_permutation_null_AUC.tsv"),sep="\t");fwrite(ptest,file.path(td,"07_permutation_test.tsv"),sep="\t")
fwrite(data.table(role=names(inputs),path=unname(inputs),sha256=vapply(inputs,digest,"",file=TRUE,algo="sha256")),file.path(qd,"input_manifest.tsv"),sep="\t")
fwrite(data.table(N=nrow(m),CRC=sum(m$CRC),control=sum(!m$CRC),cohorts=uniqueN(m$Cohort),taxa=ncol(x),repeats=10,outer_folds=5,alpha=.5,permutations=B,workers=workers),file.path(qd,"run_summary.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(qd,"sessionInfo.txt"));print(repmet);print(summarydt);print(ensemble);print(ptest)
