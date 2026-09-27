#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(MMUPHin)})
setDTthreads(1)
a<-commandArgs(FALSE); script<-normalizePath(sub('^--file=','',a[grepl('^--file=',a)]))
fig<-dirname(dirname(script));root<-dirname(dirname(fig));out<-file.path(fig,'Tables')
qa<-file.path(fig,'QA/Core879');stopifnot(!dir.exists(qa));dir.create(qa,recursive=TRUE)
writeLines('RUNNING',file.path(qa,'STATUS.txt'))
w<-function(x,n)fwrite(x,file.path(out,paste0('Fig1F_',n,'.tsv')),sep='\t',quote=FALSE,na='NA')
wm<-function(x,n)w(data.table(feature=rownames(x),as.data.table(x)),n)
input<-file.path(root,'02_Analysis/STITCH_879_Exclude6_20260907_v1/inputs')
keep<-readLines(file.path(input,'keep_runs.txt'));meta<-fread(file.path(input,'metadata.tsv'))
meta<-meta[match(keep,IID)]
stopifnot(length(keep)==879,!anyDuplicated(keep),!anyDuplicated(meta$BioSample),identical(meta$IID,keep),sum(meta$condition=='CRC')==434,sum(meta$condition=='CTR')==445)
valid<-!is.na(meta$local_sample_label)&nzchar(meta$local_sample_label)
stopifnot(!anyDuplicated(paste(meta$Cohort[valid],meta$local_sample_label[valid],sep=':')))
canonical<-file.path(root,'01_Data_Final/01_Figure1/Legacy_01_Data/01_Taxon/bacteria_species.txt')
rawpath<-file.path(root,'02_Analysis/BioSampleDedup_20260906_v1/inputs/04_species_abundance_unadjusted.tsv')
orig<-fread(canonical); ids<-sub('.bac.bracken.report','',names(orig)[-1],fixed=TRUE)
stopifnot(!anyDuplicated(ids),all(keep %in% ids))
x<-as.matrix(orig[,match(keep,ids)+1L,with=FALSE]);features<-orig[[1]];rownames(x)<-features;colnames(x)<-keep
raw<-fread(rawpath);stopifnot(identical(raw[[1]],features),all(keep %in% names(raw)),all(unname(x)==unname(as.matrix(raw[,keep,with=FALSE]))))
stopifnot(nrow(x)==2092,!anyDuplicated(features),all(is.finite(x)),all(x>=0),all(colSums(x)>0))
rm(raw,orig);gc()
rel<-sweep(x,2,colSums(x),'/');stopifnot(max(abs(colSums(rel)-1))<1e-12)
md<-as.data.frame(meta);md$Cohort<-factor(md$Cohort);md$Study_condition<-factor(md$condition,levels=c('CTR','CRC'));rownames(md)<-keep
design<-model.matrix(~Cohort+Study_condition,md);stopifnot(qr(design)$rank==ncol(design))
set.seed(20260906)
message('Starting MMUPHin: ',Sys.time())
fit<-adjust_batch(feature_abd=rel,batch='Cohort',covariates='Study_condition',data=md,control=list(verbose=FALSE))
message('MMUPHin finished: ',Sys.time())
adj<-as.matrix(fit$feature_abd_adj);stopifnot(identical(dimnames(adj),dimnames(rel)),all(is.finite(adj)))
neg<-sum(adj<0);adj[adj<0]<-0
stats<-data.table(feature=features,genome_id=sub('.*s__([^:]+).*','\\1',features),species=sub('.*:','',features),n_detected=rowSums(adj>=1e-5),prevalence_adjusted=rowMeans(adj>=1e-5),mean_adjusted=rowMeans(adj),prevalence_raw=rowMeans(rel>=1e-5),mean_raw=rowMeans(rel))
stats[,core:=prevalence_adjusted>=.20 & mean_adjusted>=.001]
stats[,failure_reason:=fifelse(core,'pass',fifelse(prevalence_adjusted<.2 & mean_adjusted<.001,'prevalence_and_mean',fifelse(prevalence_adjusted<.2,'prevalence','mean_abundance')))]
base<-file.path(root,'02_Analysis/BioSampleDedup_Reanalysis_20260906_v1')
oldpaths<-c(old225=file.path(root,'03_Results_8/Figure1/Tables/06_Fig1F/Fig1F_core_species_table.txt'),old218=file.path(base,'02_Core_rSGB_v2/07_trait_map.tsv'),old200=file.path(base,'13_Exclude5Projects_Core933_20260907_v1/07_trait_map.tsv'))
comparisons<-list()
for(label in names(oldpaths)){
 old<-fread(oldpaths[[label]],select='feature')$feature
 stopifnot(!anyDuplicated(old),all(old %in% features))
 stats[,(label):=feature %in% old]
 comparisons[[label]]<-data.table(reference=label,previous_n=length(old),new_n=sum(stats$core),shared=sum(stats$core & stats$feature %in% old),added=sum(stats$core & !stats$feature %in% old),removed=sum(!stats$core & stats$feature %in% old))
}
h18path<-file.path(base,'11_Exclude5Projects_18Taxa_20260907_v1/inputs/04_historical18_trait_map.tsv')
h18<-fread(h18path);stats[,historical18:=feature %in% h18$feature]
aniPath<-file.path(root,'01_Data_Final/01_Figure1/Legacy_01_Data/gtdbtk.Bac120.summary.tsv')
ani<-fread(aniPath,na.strings=c('NA','N/A',''),quote='')
stopifnot(!anyDuplicated(ani$user_genome))
z<-match(stats$genome_id,ani$user_genome)
maximum<-pmax(as.numeric(ani$fastani_ani[z]),as.numeric(ani$closest_placement_ani[z]),na.rm=TRUE);maximum[!is.finite(maximum)]<-NA_real_
stats[,max_ani:=maximum];stats[,classification:=ani$classification[z]]
stats[,novelty_detail:=fifelse(is.na(max_ani),'Unclassified',fifelse(max_ani>=95,'Known','Novel'))]
stats[,new_trait:=NA_character_];stats[core==TRUE,new_trait:=sprintf('r879_%04d',seq_len(.N))]
wm(rel,'all_rSGB_relative_abundance879');wm(adj,'all_rSGB_MMUPHin_adjusted879')
w(stats,'all_rSGB_core_selection879');w(stats[core==TRUE][order(-mean_adjusted)],'core_species879')
wm(adj[stats$core,,drop=FALSE],'core_adjusted_abundance879')
w(rbindlist(comparisons),'core_comparison879')
w(stats[core|old225|old218|old200],'core_membership_changes879')
w(stats[match(h18$feature,feature)],'historical18_core_status879')
w(meta,'sample_manifest879');w(meta[,.N,by=.(Cohort,condition,BioProject)],'cohort_counts879')
w(stats[core==TRUE,.N,by=novelty_detail],'core_novelty_counts879')
bycohort<-rbindlist(lapply(unique(meta$Cohort),function(co){s<-meta[Cohort==co,IID];data.table(Cohort=co,n=length(s),feature=features,prevalence_adjusted=rowMeans(adj[,s,drop=FALSE]>=1e-5),mean_adjusted=rowMeans(adj[,s,drop=FALSE]))}))
w(bycohort,'core_cohort_statistics879')
saveRDS(list(control=fit$control,seed=20260906),file.path(qa,'MMUPHin_control.rds'))
capture.output(sessionInfo(),file=file.path(qa,'sessionInfo.txt'))
paths<-c(canonical,rawpath,file.path(input,'keep_runs.txt'),file.path(input,'metadata.tsv'),oldpaths,h18path,aniPath,script)
hashes<-vapply(paths,function(p)strsplit(system2('/usr/bin/shasum',c('-a','256',shQuote(p)),stdout=TRUE),' ')[[1]][1],character(1))
fwrite(data.table(path=paths,sha256=hashes),file.path(qa,'source_sha256.tsv'),sep='\t')
writeLines(c(paste('negative_values_clamped',neg),paste('sample_labels_missing',sum(!valid)),paste('core',sum(stats$core))),file.path(qa,'checks.txt'))
print(rbindlist(comparisons));print(stats[core==TRUE,.N,by=novelty_detail]);print(stats[genome_id=='ERR1293731.10']);print(stats[historical18==TRUE & core==FALSE,.(species,genome_id,prevalence_adjusted,mean_adjusted,failure_reason)])
writeLines('COMPLETE_PENDING_INDEPENDENT_QC',file.path(qa,'STATUS.txt'))
