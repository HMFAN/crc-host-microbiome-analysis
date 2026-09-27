options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table));setDTthreads(1)
o <- '03_Results_9/Figure4'
i <- file.path(o,'Tables','Inputs','MR_LeadProxy')
t <- file.path(o,'Tables','MR','LeadProxy')
q <- file.path(o,'QA','MR_LeadProxy')
rd <- function(x) {
  part <- strsplit(x, '/', fixed=TRUE)[[1]]
  base <- if (part[1] == 'Inputs') i else if (part[1] == 'Tables') t else stop('Unsupported table scope: ', part[1])
  fread(file.path(base, paste(part[-1], collapse='/')))
}
wr <- function(x,n)fwrite(x,file.path(t,n),sep='\t',na='NA')
l <- rd('Inputs/current_retained_71.tsv');p <- rd('Inputs/proxy_exposure.tsv')
v <- rd('Tables/variant_build_mapping.tsv');stopifnot(!anyDuplicated(v$SNP))
d <- l[,.(trait,species,lead_SNP=SNP,SNP,A1,bx=BETA,sx=SE,exposure_P=P,exposure_AF=AF,exposure_N=N,match_type='Direct',LD_r=1,LD_r2=1,distance=0)]
pr <- p[,.(trait,species,lead_SNP,SNP=proxy_SNP,A1,bx=BETA,sx=SE,exposure_P=P,exposure_AF=AF,exposure_N=N,match_type='Proxy',LD_r,LD_r2,distance)]
a0 <- rbindlist(list(d,pr),use.names=TRUE)
a0[,`:=`(F_stat=(bx/sx)^2,record_id=.I)]
stopifnot(!anyDuplicated(a0[,.(trait,lead_SNP,SNP)]),all(a0$exposure_N==879),all(pr$LD_r2>=0.8-1e-10),all(pr$SNP!=pr$lead_SNP))
a0 <- merge(a0,v,by='SNP',all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(a0$chr38),all(a0$A1==a0$ref38|a0$A1==a0$alt38),all(is.finite(a0$exposure_AF)),all(a0$exposure_AF>0&a0$exposure_AF<1))
comp <- c(A='T',T='A',C='G',G='C')
outcomes <- c('GCST012879','GCST90129505','FinnGen_DF13')
y <- rd('Inputs/outcomes_normalized.tsv')
stopifnot(all(c('outcome','chr','pos','effect_allele','other_allele','beta','se','eaf','rsid','build')%in%names(y)))
y[,`:=`(chr=as.character(chr),pos=as.integer(pos),effect_allele=toupper(effect_allele),other_allele=toupper(other_allele))]
y[,key:=paste(chr,pos,sep=':')]
audit <- list()
for(oo in outcomes){
 yy <- y[outcome==oo];build <- unique(yy$build)
 stopifnot(nrow(yy)>0,length(build)==1,build%in%c(37L,38L))
 for(i in seq_len(nrow(a0))){
  z <- copy(a0[i]);ea<-z$A1;oa<-if(ea==z$alt38)z$ref38 else z$alt38
  query_key <- if(build==38)z$key38 else z$key37
  if(build==37&&z$lift_count==1&&z$strand37=='-'){ea<-unname(comp[ea]);oa<-unname(comp[oa])}
  rr <- yy[yy$key==query_key]
  st<-'no_outcome_record';sg<-NA_real_;al<-NA_character_;j<-NA_integer_;fd<-NA_real_
  if(build==37&&z$lift_count!=1){st<-'nonunique_build_mapping';rr<-rr[0]}
  if(nrow(rr)>0){
   same<-rr$effect_allele==ea&rr$other_allele==oa;swap<-rr$effect_allele==oa&rr$other_allele==ea
   csame<-rr$effect_allele==unname(comp[ea])&rr$other_allele==unname(comp[oa])
   cswap<-rr$effect_allele==unname(comp[oa])&rr$other_allele==unname(comp[ea])
   ix<-which(same|swap|csame|cswap);st<-'allele_mismatch'
   if(length(ix)>1)st<-'ambiguous_multiple_records'
   if(length(ix)==1){
    j<-ix[1];pal<-paste(sort(c(ea,oa)),collapse='')%in%c('AT','CG')
    if(pal){
     af<-rr$eaf[j];xf<-z$exposure_AF
     st<-'palindrome_no_outcome_EAF'
     if(is.finite(af)&&is.finite(xf)&&af>0&&af<1){
      st<-'palindrome_ambiguous_frequency'
      diffs<-c(abs(xf-af),abs(xf-(1-af)))
      if(min(xf,1-xf)<.42&&min(af,1-af)<.42&&min(diffs)<=.10&&abs(diff(diffs))>.10){
       sg<-if(which.min(diffs)==1)1 else -1;al<-if(sg==1)'palindrome_EAF_same' else 'palindrome_EAF_swap'
       fd<-min(diffs);st<-'harmonized'
      }
     }
    }else{
     sg<-if(same[j]|csame[j])1 else -1
     al<-if(same[j])'same' else if(swap[j])'swap' else if(csame[j])'strand_same' else 'strand_swap'
     st<-'harmonized'
    }
   }
  }
  z[,`:=`(outcome=oo,outcome_build=build,outcome_key=query_key,ea_outcome_build=ea,oa_outcome_build=oa,status=st,alignment=al,orientation_sign=sg,EAF_distance=fd,
          outcome_rsid=if(is.na(j))NA_character_ else as.character(rr$rsid[j]),
          outcome_effect_allele=if(is.na(j))NA_character_ else rr$effect_allele[j],
          outcome_other_allele=if(is.na(j))NA_character_ else rr$other_allele[j],
          outcome_EAF=if(is.na(j))NA_real_ else rr$eaf[j],
          by_raw=if(is.na(j))NA_real_ else rr$beta[j],
          by=if(is.na(j))NA_real_ else sg*rr$beta[j],sy=if(is.na(j))NA_real_ else rr$se[j])]
  if(z$status=='harmonized'&&(!is.finite(z$by)||!is.finite(z$sy)||z$sy<=0))z$status<-'invalid_outcome_effect'
  if(z$status=='harmonized'&&(!is.finite(z$bx)||!is.finite(z$sx)||z$sx<=0||z$bx==0))z$status<-'invalid_exposure_effect'
  if(z$status=='harmonized'&&z$F_stat<=10)z$status<-'weak_instrument_F_le_10'
  audit[[length(audit)+1]]<-z
 }
}
a <- rbindlist(audit)
a[,`:=`(wald_beta=by/bx,first_order_SE=sy/abs(bx),delta_SE=sqrt(sy^2/bx^2+by^2*sx^2/bx^4))]
a[,`:=`(first_order_P=2*pnorm(-abs(wald_beta/first_order_SE)),delta_P=2*pnorm(-abs(wald_beta/delta_SE)))]
a[,inference_eligible:=status=='harmonized']
a[inference_eligible==FALSE,`:=`(wald_beta=NA_real_,first_order_SE=NA_real_,delta_SE=NA_real_,first_order_P=NA_real_,delta_P=NA_real_)]
wr(a,'harmonization_all_leads_proxies.tsv')
wr(a[,.(records=.N,unique_variants=uniqueN(SNP)),by=.(outcome,match_type,status)],'matching_summary.tsv')
g <- rd('Inputs/proxy_alt_dosages.tsv');stopifnot(nrow(g)==879,!anyDuplicated(g$IID))
meta <- fread('02_Analysis/Heri_mGWAS_879_20260909_v1/inputs/metadata.tsv')
stopifnot(setequal(g$IID,meta$IID),all(unique(a$SNP)%in%names(g)))
ldfun <- function(x,y){
 if(x==y)return(1)
 if(sub(':.*','',x)!=sub(':.*','',y))return(0)
 r<-cor(g[[x]],g[[y]],use='complete.obs');stopifnot(is.finite(r));r^2
}
pc <- unique(pr[,.(lead_SNP,SNP,LD_r,LD_r2)])
pc[,measured_r:=mapply(function(x,y)cor(g[[x]],g[[y]],use='complete.obs'),lead_SNP,SNP)]
pc[,`:=`(r_error=abs(measured_r-LD_r),r2_error=abs(measured_r^2-LD_r2))]
stopifnot(max(pc$r_error)<1e-5,max(pc$r2_error)<1e-5)
wr(pc,'proxy_LD_independent_check.tsv')
strategies <- c('Lead_only','Proxy_only','Lead_plus_proxy')
selection<-list();ldlog<-list();result<-list();pathaudit<-list()
for(ss in strategies)for(oo in outcomes)for(tt in unique(l$trait)){
 base<-a[trait==tt&outcome==oo];chosen<-base[0]
 for(lead in l[trait==tt,SNP]){
  dz<-base[lead_SNP==lead&match_type=='Direct'&status=='harmonized']
  pz<-base[lead_SNP==lead&match_type=='Proxy'&status=='harmonized']
  setorderv(pz,c('LD_r2','distance','SNP'),c(-1,1,1))
  pick<-if(ss=='Lead_only')dz else if(ss=='Proxy_only')head(pz,1) else if(nrow(dz)>0)dz else head(pz,1)
  pathaudit[[length(pathaudit)+1]]<-data.table(strategy=ss,outcome=oo,trait=tt,lead_SNP=lead,
   direct_status=base[lead_SNP==lead&match_type=='Direct',status],
   proxy_candidates=base[lead_SNP==lead&match_type=='Proxy',.N],eligible_proxies=nrow(pz),
   selected_SNP_before_LD=if(nrow(pick))pick$SNP[1] else NA_character_,
   selection_status=if(nrow(pick))'selected_before_within_trait_LD' else if(ss=='Lead_only')'no_eligible_direct' else if(ss=='Proxy_only')'no_eligible_proxy' else 'no_eligible_direct_or_proxy')
  if(nrow(pick)>0){pick[,lead_selection_P:=l[trait==tt&SNP==lead,P]];chosen<-rbind(chosen,pick,fill=TRUE)}
 }
 if(nrow(chosen)>0){
  setorderv(chosen,c('lead_selection_P','lead_SNP'),c(1,1));chosen[,LD_retained:=FALSE];kept<-integer()
  for(i in seq_len(nrow(chosen))){
   pass<-TRUE
   for(jj in kept){
    r2<-ldfun(chosen$SNP[i],chosen$SNP[jj]);reject<-r2>=.001
    ldlog[[length(ldlog)+1]]<-data.table(strategy=ss,outcome=oo,trait=tt,SNP1=chosen$SNP[i],SNP2=chosen$SNP[jj],r2=r2,rejected=reject)
    if(reject)pass<-FALSE
   }
   if(pass){chosen$LD_retained[i]<-TRUE;kept<-c(kept,i)}
  }
  chosen[,strategy:=ss];selection[[length(selection)+1]]<-chosen
 }
 z<-if(nrow(chosen))chosen[LD_retained==TRUE] else chosen;n<-nrow(z)
 be<-se<-se1<-pv<-p1<-Q<-QP<-NA_real_;meth<-'not_testable';ins<-NA_character_
 if(n==1){be<-z$wald_beta;se<-z$delta_SE;se1<-z$first_order_SE;meth<-'Wald_delta'}
 if(n>1){
  den<-sum(z$bx^2/z$sy^2);be<-sum(z$bx*z$by/z$sy^2)/den
  Q<-sum((z$by-be*z$bx)^2/z$sy^2);QP<-pchisq(Q,n-1,lower.tail=FALSE)
  se<-sqrt(1/den)*max(1,sqrt(Q/(n-1)));se1<-se;meth<-'IVW_multiplicative_RE'
 }
 if(n>0){pv<-2*pnorm(-abs(be/se));p1<-2*pnorm(-abs(be/se1));ins<-if(all(z$match_type=='Direct'))'Direct' else if(all(z$match_type=='Proxy'))'Proxy' else 'Mixed'}
 result[[length(result)+1]]<-data.table(strategy=ss,trait=tt,species=base$species[1],outcome=oo,candidate_leads=uniqueN(base$lead_SNP),nsnp=n,n_proxy=if(n)sum(z$match_type=='Proxy') else 0,
  method=meth,instrument_status=ins,beta=be,SE=se,P=pv,OR=exp(be),CI_low=exp(be-1.96*se),CI_high=exp(be+1.96*se),
  first_order_SE=se1,first_order_P=p1,Q=Q,Q_P=QP,F_min=if(n)min(z$F_stat) else NA_real_,
  instruments=if(n)paste(z$SNP,collapse=';') else NA_character_,
  lead_signals=if(n)paste(z$lead_SNP,collapse=';') else NA_character_,
  status=if(n)'testable' else 'no_eligible_harmonized_instrument')
}
r<-rbindlist(result);sel<-rbindlist(selection,fill=TRUE)
adj<-function(x,n){y<-rep(NA_real_,length(x));i<-which(is.finite(x));y[i]<-p.adjust(x[i],'BH',n=n);y}
r[,FDR57:=adj(P,57),by=.(strategy,outcome)]
r[,FDR171:=adj(P,171),by=strategy]
r[,first_order_FDR57:=adj(first_order_P,57),by=.(strategy,outcome)]
stopifnot(nrow(r)==513,!anyDuplicated(r[,.(strategy,trait,outcome)]))
wr(r,'MR_taxon_results.tsv');wr(sel,'selected_instruments.tsv');wr(rbindlist(ldlog),'within_trait_LD.tsv');wr(rbindlist(pathaudit),'lead_path_selection_audit.tsv')
wr(a[match_type=='Direct'],'lead_Wald_all71_three_outcomes.tsv')
wr(a[match_type=='Proxy'],'proxy_Wald_all_candidates_three_outcomes.tsv')
wr(r[,.(planned=57L,testable=sum(nsnp>0),nominal=sum(P<.05,na.rm=TRUE),FDR57_significant=sum(FDR57<.05,na.rm=TRUE),FDR171_significant=sum(FDR171<.05,na.rm=TRUE)),by=.(strategy,outcome)],'MR_summary.tsv')
wr(r[trait=='core161'],'Bc_MR_summary.tsv')
fwrite(data.table(item=c('unique_taxa','lead_associations','unique_leads','proxy_associations','unique_proxies','outcomes','strategies','result_rows'),value=c(57,71,67,nrow(pr),uniqueN(pr$SNP),3,3,nrow(r))),file.path(q,'counts.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(q,'MR_sessionInfo.txt'))
cat('COMPLETE: 513 planned taxon x outcome x strategy rows\n');print(rd('Tables/MR_summary.tsv'));print(r[trait=='core161',.(strategy,outcome,nsnp,OR,CI_low,CI_high,P,FDR57)])
