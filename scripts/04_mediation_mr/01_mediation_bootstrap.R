options(stringsAsFactors=FALSE);suppressPackageStartupMessages(library(data.table));setDTthreads(1)
o<-'03_Results_9/Figure4/Tables';od<-file.path(o,'Mediation/AllCandidatesBootstrap10000');dir.create(od,recursive=TRUE,showWarnings=FALSE);b<-'02_Analysis/Heri_mGWAS_879_20260909_v1';l<-fread(file.path(b,'linear_results/cohort_robustness/candidate_review_v2/retained_71_signals.tsv'));m<-fread(file.path(b,'inputs/metadata.tsv'));p<-fread(file.path(b,'inputs/phenotypes.tsv'));pc<-fread(file.path(o,'Inputs/02_Variants/DS_PCA.eigenvec'));g<-fread('02_Analysis/CRC_Disease_Stratified_Effects/results/lead_dosages.raw');p<-p[match(m$IID,IID)];pc<-pc[match(m$IID,IID)];g<-g[match(m$IID,IID)];pn<-paste0('PC',1:10)
cov<-data.frame(cohort=factor(m$Cohort),pc[,..pn]);full<-data.frame(cov,Age=suppressWarnings(as.numeric(m$Age)),BMI=suppressWarnings(as.numeric(m$BMI)),Sex=factor(ifelse(tolower(m$Sex)%in%c('m','male'),'M',ifelse(tolower(m$Sex)%in%c('f','female'),'F',NA))));keep<-complete.cases(full)&full$Age>0&full$BMI>0;Y<-as.integer(m$condition=='CRC');stopifnot(sum(keep)==605)
source('scripts/04_mediation_mr/mediation_calc.R')

stopifnot(nrow(l)==71,uniqueN(l$trait)==57,nrow(m)==879,!anyDuplicated(m$IID))
for(z in list(p,pc,g))stopifnot(identical(as.character(m$IID),as.character(z$IID)))
jobs<-CJ(model=c('Primary','Interaction','Full_covariates'),i=seq_len(nrow(l)),sorted=FALSE)
writeLines('RUNNING',file.path(od,'STATUS.txt'))
fwrite(jobs[,.(job=.I,model,i,trait=l$trait[i],SNP=l$SNP[i],seed=85000+.I)],file.path(od,'jobs.tsv'),sep='\t')
run<-function(j){
 fn<-file.path(od,sprintf('job_%03d.rds',j));if(file.exists(fn))return(readRDS(fn)$result)
 mode<-jobs$model[j];idx<-jobs$i[j];set.seed(85000+j)
 ii<-if(mode=='Full_covariates')which(keep) else seq_len(nrow(m))
 X<-model.matrix(~.,droplevels(if(mode=='Full_covariates')full[ii,] else cov[ii,]))
 d<-round(g[[paste0(l$SNP[idx],'_',l$A1[idx])]][ii]*16384)/16384;M<-p[[l$trait[idx]]][ii];y<-Y[ii]
 stopifnot(all(is.finite(d)),all(is.finite(M)))
 t0<-if(mean(d)>1.5)1 else 0;int<-mode=='Interaction';point<-calc(X,d,M,y,int,t0)
 strata<-split(seq_along(ii),interaction(m$Cohort[ii],y,drop=TRUE))
 draws<-matrix(NA_real_,10000,3,dimnames=list(NULL,c('ACME','ADE','total')))
 for(k in seq_len(10000)){
 ix<-unlist(lapply(strata,function(v)sample(v,length(v),replace=TRUE)),use.names=FALSE)
 draws[k,]<-tryCatch(calc(X[ix,,drop=FALSE],d[ix],M[ix],y[ix],int,t0),error=function(e)c(NA,NA,NA))
 }
 rr<-rbindlist(lapply(names(point),function(ef){v<-draws[,ef];v<-v[is.finite(v)];ci<-quantile(v,c(.025,.975));pv<-min(1,2*min((sum(v<=0)+1)/(length(v)+1),(sum(v>=0)+1)/(length(v)+1)))
 data.table(model=mode,trait=l$trait[idx],species=l$species[idx],SNP=l$SNP[idx],N=length(ii),contrast_from=t0,contrast_to=t0+1,effect=ef,estimate=point[ef],CI_low=ci[1],CI_high=ci[2],bootstrap_success=length(v),nonpositive=sum(v<=0),nonnegative=sum(v>=0),P=pv)}))
 saveRDS(list(result=rr,point=point,draws=draws,seed=85000+j),paste0(fn,'.tmp'));file.rename(paste0(fn,'.tmp'),fn)
 cat('DONE',j,mode,l$trait[idx],format(Sys.time()),'\n');flush.console();rr
}
rr<-parallel::mclapply(seq_len(nrow(jobs)),run,mc.cores=6,mc.preschedule=FALSE)
stopifnot(!any(vapply(rr,inherits,logical(1),'try-error')))
r<-rbindlist(rr);stopifnot(nrow(r)==639)
r[,`:=`(BH=p.adjust(P,'BH'),Holm=p.adjust(P,'holm'),Bonferroni=p.adjust(P,'bonferroni')) ,by=.(model,effect)]
r[,`:=`(BH_all213=p.adjust(P,'BH'),Holm_all213=p.adjust(P,'holm'),Bonferroni_all213=p.adjust(P,'bonferroni')),by=effect]
fwrite(r,file.path(od,'all71_three_models.tsv'),sep='\t')
a<-r[effect=='ACME']
corrections<-c('BH','Holm','Bonferroni','BH_all213','Holm_all213','Bonferroni_all213')
s<-rbindlist(lapply(corrections,function(cc)a[,.(correction=cc,paths=sum(get(cc)<.05),taxa=uniqueN(trait[get(cc)<.05])),by=model]))
w<-dcast(a,trait+species+SNP~model,value.var=c('estimate',corrections,'P'))
for(cc in corrections){cols<-paste0(cc,'_',c('Primary','Interaction','Full_covariates'));w[,(paste0('robust_',cc)):=Reduce(`&`,lapply(.SD,function(x)x<.05)) & sign(estimate_Primary)==sign(estimate_Interaction) & sign(estimate_Primary)==sign(estimate_Full_covariates),.SDcols=cols]}
s<-rbind(s,rbindlist(lapply(corrections,function(cc){v<-w[[paste0('robust_',cc)]];data.table(model='All_three_same_path_same_direction',correction=cc,paths=sum(v),taxa=uniqueN(w$trait[v]))})))
fwrite(s,file.path(od,'retention_counts.tsv'),sep='\t');fwrite(w,file.path(od,'all71_robustness.tsv'),sep='\t')
fwrite(w[robust_BH==TRUE],file.path(od,'robust_BH_paths.tsv'),sep='\t')
fwrite(w[robust_Bonferroni==TRUE],file.path(od,'robust_Bonferroni_paths.tsv'),sep='\t')
writeLines(if(all(r$bootstrap_success==10000))'COMPLETE_PENDING_INDEPENDENT_QA' else 'COMPLETE_WITH_BOOTSTRAP_FAILURES',file.path(od,'STATUS.txt'));print(s)
