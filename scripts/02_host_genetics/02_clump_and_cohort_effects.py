#!/usr/bin/python3
import csv,json,math,statistics,os,time,signal,subprocess,fcntl,traceback,hashlib
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
H=Path(os.environ['CRC_GENETICS_ROOT']);L=H/'linear_mgwas190_v1';O=L/'P1e6_clump_cohort_v1';P=H/'software/plink2';stage='preflight'
def write(p,fields,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def status(state,**kw):
 d=dict(state=state,stage=stage,pid=os.getpid(),updated=time.strftime('%F %T %z'),**kw);p=O/'status.tmp';p.write_text(json.dumps(d,indent=2));p.replace(O/'status.json');print(json.dumps(d),flush=True)
def resource():
 m={l.split(':')[0]:int(l.split()[1]) for l in Path('/proc/meminfo').read_text().splitlines()};return os.getloadavg()[0],m['MemAvailable']/1048576

def run(name,args):
 global stage
 stage=name
 if (O/'logs'/(name+'.done')).exists():return
 while True:
  load,mem=resource()
  if load<85 and mem>96:break
  status('WAITING_RESOURCES',load=load,available_GiB=mem);time.sleep(30)
 status('RUNNING',command=list(map(str,args)))
 with (O/'logs'/(name+'.stdout')).open('w') as f:
  p=subprocess.Popen(list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,start_new_session=True);paused=False
  while p.poll() is None:
   load,mem=resource()
   if not paused and (load>104 or mem<48):
    try:os.killpg(p.pid,signal.SIGSTOP)
    except ProcessLookupError:pass
    paused=True;status('PAUSED_RESOURCES')
   elif paused and load<80 and mem>96:
    try:os.killpg(p.pid,signal.SIGCONT)
    except ProcessLookupError:pass
    paused=False;status('RUNNING')
   time.sleep(5)
  assert p.returncode==0,(name,p.returncode)
 (O/'logs'/(name+'.done')).write_text(time.strftime('%F %T'))

def qsf(a,x):
 if x<=0:return 1.
 if x<a+1:
  ap=a;term=1/a;s=term
  for _ in range(10000):
   ap+=1;term*=x/ap;s+=term
   if abs(term)<abs(s)*1e-14:break
  return max(0.,min(1.,1-s*math.exp(-x+a*math.log(x)-math.lgamma(a))))
 b=x+1-a;c=1e300;d=1/b;h=d
 for i in range(1,10000):
  an=-i*(i-a);b+=2;d=an*d+b
  if abs(d)<1e-300:d=1e-300
  c=b+an/c
  if abs(c)<1e-300:c=1e-300
  d=1/d;delta=d*c;h*=delta
  if abs(delta-1)<1e-14:break
 return max(0.,min(1.,math.exp(-x+a*math.log(x)-math.lgamma(a))*h))

def main():
 global stage
 O.mkdir(exist_ok=True)
 lock=open(O/'lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 for d in ['inputs','clumps','logs','cohorts','summary']:(O/d).mkdir(exist_ok=True)
 status('RUNNING')
 assert json.loads((L/'VALIDATED.json').read_text())['traits']==190
 hits=json.loads((L/'suggestive_hits.json').read_text());hits=[r for r in hits if float(r['P'])<1e-6]
 traits=sorted({r['trait'] for r in hits});assert len(hits)==787 and len(traits)==62 and len({r['ID'] for r in hits})==747
 mapping={r['trait']:r for r in csv.DictReader((H/'inputs/trait_map.tsv').open(),delimiter='\t')}
 meta=list(csv.DictReader((H/'inputs/metadata.tsv').open(),delimiter='\t'));mm={r['IID']:r for r in meta};base=list(csv.DictReader((H/'mgwas/covariates.tsv').open(),delimiter='\t'))
 assert len(mm)==len(base)==879 and {r['IID'] for r in base}==set(mm)
 spec=spec_from_file_location('helpers',Path(__file__).with_name('rank_helper.py'));helpers=module_from_spec(spec);spec.loader.exec_module(helpers)
 manifest=dict(source='completed all190 DS linear regression',selection='nominal P<1e-6, post-results exploratory threshold',traits=62,candidate_pairs=787,unique_SNPs=747,clump=dict(p1=1e-6,p2=.01,r2=.1,kb=1000,LD='unphased DS in all879'),sensitivity='5 within-cohort and 5 leave-one-cohort-out regressions at original lead SNPs',phenotypes='frozen original CLR+RINT, no re-transform within subsets',covariates='PC1-10+CRC; remaining cohort dummies for leave-one-out',sample_filter='no additional sample exclusion; same879',pass_fail='no automatic pass/fail classification; report uniform evidence for all leads',Bc='same rules as all other traits',threads=4)
 (O/'methods.json').write_text(json.dumps(manifest,indent=2))
 leads=[]
 for idx,t in enumerate(traits):
  stage='prepare_clump_'+t;status('PREPARING_CLUMP',completed_traits=idx,total_traits=62)
  assoc=O/'inputs'/(t+'.assoc');primary={r['ID']:r for r in hits if r['trait']==t}
  if not (O/'logs'/('clump_'+t+'.done')).exists():
   with (L/('all190.'+t+'.glm.linear')).open() as f,assoc.open('w') as out:
    w=csv.writer(out,delimiter='\t');w.writerow(['SNP','P'])
    for r in csv.DictReader(f,delimiter='\t'):
     if r['P']!='NA' and float(r['P'])<=.01:w.writerow([r['ID'],r['P']])
  run('clump_'+t,[P,'--pfile',H/'genotypes/known_DS','--clump',assoc,'--clump-id-field','SNP','--clump-p-field','P','--clump-p1','1e-6','--clump-p2','.01','--clump-r2','.1','--clump-kb','1000','--clump-unphased','--threads','4','--memory','16384','--out',O/'clumps'/t])
  with (O/'clumps'/(t+'.clumps')).open() as f:
   for r in csv.DictReader(f,delimiter='\t'):
    a=primary[r['ID']];leads.append(dict(trait=t,species=mapping[t]['species'],SNP=r['ID'],A1=a['A1'],BETA=a['BETA'],SE=a['SE'],P=a['P'],N=a['OBS_CT'],AF=a['A1_FREQ'],clump_total=r['TOTAL']))
 write(O/'summary/lead_pairs.tsv',list(leads[0]),leads)
 write(O/'summary/Bc_leads.tsv',list(leads[0]),[r for r in leads if r['trait']=='core161'])
 snps=sorted({r['SNP'] for r in leads});(O/'inputs/leads.ids').write_text('\n'.join(snps)+'\n')
 run('extract_lead_genotypes',[P,'--pfile',H/'genotypes/known_DS','--extract',O/'inputs/leads.ids','--make-pgen','--threads','4','--memory','16384','--out',O/'inputs/leads'])
 ph=list(csv.DictReader((H/'inputs/phenotypes.tsv').open(),delimiter='\t'));write(O/'inputs/phenotypes.tsv',['FID','IID']+traits,[{k:r[k] for k in ['FID','IID']+traits} for r in ph])
 cohorts=sorted({r['cohort'] for r in base});assert len(cohorts)==5
 allrows=[];qc=[]
 for mode in ['within','leaveout']:
  for c in cohorts:
   name=mode+'_'+c;selected=[r for r in base if (r['cohort']==c if mode=='within' else r['cohort']!=c)]
   levels=sorted({r['cohort'] for r in selected});names=['FID','IID']+['PC'+str(i) for i in range(1,11)]+['CRC_status']+['cohort_'+v for v in levels[1:]];cv=[]
   for r in selected:
    x={k:r[k] for k in names[:12]};x['CRC_status']=int(r['CRC_status']=='CRC');x.update({'cohort_'+v:int(r['cohort']==v) for v in levels[1:]});cv.append(x)
   constant=[k for k in names[2:] if len({float(r[k]) for r in cv})==1];fields=[k for k in names if k not in constant]
   cols=[[1.]*len(cv)]+[[float(r[k]) for r in cv] for k in fields[2:]];rank=helpers.rank(cols);assert rank==len(cols),(name,rank,len(cols))
   if name=='within_CHN3':
    write(O/'inputs'/(name+'.original_covariates.tsv'),fields,[{k:r[k] for k in fields} for r in cv])
    basis=[]
    for col in cols:
     norm=math.sqrt(sum(v*v for v in col));v=[x/norm for x in col]
     for _ in range(2):
      for bvec in basis:
       dot=sum(x*y for x,y in zip(v,bvec));v=[x-dot*y for x,y in zip(v,bvec)]
     norm=math.sqrt(sum(x*x for x in v));assert norm>1e-9
     basis.append([x/norm for x in v])
    assert max(abs(sum(x*y for x,y in zip(a,b))) for i,a in enumerate(basis) for j,b in enumerate(basis) if i!=j)<1e-10
    for col in cols:
     predicted=[0.]*len(col)
     for bvec in basis:
      dot=sum(x*y for x,y in zip(col,bvec));predicted=[x+dot*y for x,y in zip(predicted,bvec)]
     assert max(abs(x-y) for x,y in zip(col,predicted))<1e-8
    cv=[dict(FID=r['FID'],IID=r['IID'],**{'Q'+str(j):basis[j][i] for j in range(1,len(basis))}) for i,r in enumerate(cv)]
    fields=['FID','IID']+['Q'+str(j) for j in range(1,len(basis))]
   write(O/'inputs'/(name+'.covar'),fields,[{k:r[k] for k in fields} for r in cv])
   (O/'inputs'/(name+'.keep')).write_text(''.join(r['IID']+'\t'+r['IID']+'\n' for r in cv))
   qc.append(dict(analysis=name,N=len(cv),rank=rank,reference_cohort=levels[0],constant_covariates_dropped=','.join(constant)))
   run(name,[P,'--pfile',O/'inputs/leads','--keep',O/'inputs'/(name+'.keep'),'--pheno',O/'inputs/phenotypes.tsv','--covar',O/'inputs'/(name+'.covar'),'--covar-variance-standardize','--glm','hide-covar','omit-ref','--threads','4','--memory','16384','--out',O/'cohorts'/name])
   for t in traits:
    target={r['SNP'] for r in leads if r['trait']==t}
    with (O/'cohorts'/(name+'.'+t+'.glm.linear')).open() as f:
     for r in csv.DictReader(f,delimiter='\t'):
      if r['ID'] not in target:continue
      allrows.append(dict(analysis=name,mode=mode,cohort=c,trait=t,SNP=r['ID'],A1=r['A1'],BETA=r['BETA'],SE=r['SE'],P=r['P'],N=r['OBS_CT'],AF=r['A1_FREQ'],ERRCODE=r['ERRCODE']))
 write(O/'summary/covariate_QC.tsv',list(qc[0]),qc)
 write(O/'summary/all_cohort_estimates.tsv',list(allrows[0]),allrows)
 assert len(allrows)==len(leads)*10
 summary=[]
 for r in leads:
  z=[x for x in allrows if x['trait']==r['trait'] and x['SNP']==r['SNP']];assert all(x['A1']==r['A1'] for x in z)
  valid=lambda x:x['P']!='NA' and x['BETA']!='NA' and x['SE']!='NA' and math.isfinite(float(x['SE'])) and float(x['SE'])>0
  within=[x for x in z if x['mode']=='within' and valid(x)];leave=[x for x in z if x['mode']=='leaveout' and valid(x)]
  b=float(r['BETA']);d=dict(r,within_estimable=len(within),within_same_direction=sum(float(x['BETA'])*b>0 for x in within),leaveout_estimable=len(leave),leaveout_same_direction=sum(float(x['BETA'])*b>0 for x in leave),Q='NA',Q_P='NA',I2='NA',meta_BETA='NA',meta_SE='NA',meta_P='NA',leaveout_minBETA=min([float(x['BETA']) for x in leave],default='NA'),leaveout_maxBETA=max([float(x['BETA']) for x in leave],default='NA'),leaveout_maxP=max([float(x['P']) for x in leave],default='NA'))
  if len(within)>=2:
   weights=[1/float(x['SE'])**2 for x in within];betas=[float(x['BETA']) for x in within];total=sum(weights);mb=sum(w*v for w,v in zip(weights,betas))/total;se=1/math.sqrt(total);q=sum(w*(v-mb)**2 for w,v in zip(weights,betas));df=len(within)-1
   d.update(Q=q,Q_P=qsf(df/2,q/2),I2=max(0,(q-df)/q)*100 if q>0 else 0,meta_BETA=mb,meta_SE=se,meta_P=math.erfc(abs(mb/se)/math.sqrt(2)))
  summary.append(d)
 write(O/'summary/robustness_evidence.tsv',list(summary[0]),summary)
 write(O/'summary/Bc_robustness_evidence.tsv',list(summary[0]),[r for r in summary if r['trait']=='core161'])
 stage='complete';status('COMPLETE_AWAITING_REVIEW',lead_pairs=len(leads),unique_lead_SNPs=len(snps),traits=len(traits),sensitivity_runs=10,no_automatic_pass_fail=True)
 (O/'VALIDATED.json').write_text(json.dumps(dict(lead_pairs=len(leads),unique_lead_SNPs=len(snps),traits=len(traits),cohort_estimates=len(allrows),complete=time.strftime('%F %T')),indent=2))
if __name__=='__main__':
 try:main()
 except BaseException as e:
  if O.exists():status('FAILED',error=str(e))
  traceback.print_exc();raise
