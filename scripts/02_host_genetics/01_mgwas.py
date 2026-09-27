#!/usr/bin/python3
import csv,json,os,time,subprocess,signal,traceback,fcntl,math
from pathlib import Path
H=Path(os.environ['CRC_GENETICS_ROOT']);L=H/'linear_mgwas190_v1';L.mkdir(exist_ok=True)
lock=open(L/'lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def status(state,**kw):
 p=L/'status.tmp';p.write_text(json.dumps(dict(state=state,pid=os.getpid(),time=time.strftime('%F %T'),**kw),indent=2));p.replace(L/'status.json')
def run(stage,args):
 def resources():
  m={l.split(':')[0]:int(l.split()[1]) for l in Path('/proc/meminfo').read_text().splitlines()};return os.getloadavg()[0],m['MemAvailable']/1048576
 while True:
  load,mem=resources()
  if load<85 and mem>96:break
  status('WAITING_RESOURCES',stage=stage);time.sleep(30)
 status('RUNNING',stage=stage,command=list(map(str,args)))
 with (L/(stage+'.stdout')).open('w') as f:
  p=subprocess.Popen(list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,start_new_session=True);paused=False
  while p.poll() is None:
   load,mem=resources()
   if not paused and (load>104 or mem<48):
    try:os.killpg(p.pid,signal.SIGSTOP)
    except ProcessLookupError:pass
    paused=True;status('PAUSED_RESOURCES',stage=stage)
   elif paused and load<80 and mem>96:
    try:os.killpg(p.pid,signal.SIGCONT)
    except ProcessLookupError:pass
    paused=False;status('RUNNING',stage=stage)
   time.sleep(5)
  assert p.returncode==0,(stage,p.returncode)
def main():
 assert not (L/'covariates.tsv').exists(),'Versioned run already exists'
 rows=list(csv.DictReader((H/'mgwas/covariates.tsv').open(),delimiter='\t'));cohorts=sorted({r['cohort'] for r in rows})
 fields=['FID','IID']+['PC'+str(i) for i in range(1,11)]+['CRC_status']+['cohort_'+c for c in cohorts[1:]]
 with (L/'covariates.tsv').open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rows:
   d={k:r[k] for k in fields[:12]};d['CRC_status']=int(r['CRC_status']=='CRC');d.update({'cohort_'+c:int(r['cohort']==c) for c in cohorts[1:]});w.writerow(d)
 assert len(rows)==879 and len({r['IID'] for r in rows})==879
 (L/'targets.ids').write_text('6:94953369:C:T\n1:103289110:C:T\n')
 (L/'methods.json').write_text(json.dumps(dict(model='PLINK2 additive linear regression on continuous DS; original GT not substituted',phenotypes='same 190 CLR+RINT traits',samples=879,covariates=fields[2:],cohort_reference=cohorts[0],threads=8,studywide_Bonferroni=.05/(190*1568972),REGENIE='continues separately; this model has no LOCO prediction',target='Bc core161 chr6:94953369:C:T prespecified'),indent=2))
 common=[H/'software/plink2','--pfile',H/'genotypes/known_DS','--pheno',H/'inputs/phenotypes.tsv','--covar',L/'covariates.tsv','--covar-variance-standardize','--glm','hide-covar','omit-ref','--threads','8','--memory','16384']
 run('01_Bc_targets',common+['--pheno-name','core161','--extract',L/'targets.ids','--out',L/'Bc_targets'])
 target=list(csv.DictReader((L/'Bc_targets.core161.glm.linear').open(),delimiter='\t'));assert len(target)==2
 (L/'Bc_targets.json').write_text(json.dumps(target,indent=2))
 run('02_all190',common+['--out',L/'all190'])
 traits=list(csv.DictReader((H/'inputs/trait_map.tsv').open(),delimiter='\t'));summary=[];hits=[]
 for t in traits:
  n=0;valid=0;mn=1.;nom=0;sig=0
  with (L/('all190.'+t['trait']+'.glm.linear')).open() as f:
   for r in csv.DictReader(f,delimiter='\t'):
    n+=1
    if r['P']=='NA':continue
    p=float(r['P']);assert math.isfinite(p);valid+=1;mn=min(mn,p)
    if p<1e-5:nom+=1;hits.append(dict(trait=t['trait'],species=t['species'],**r))
    if p<.05/(190*1568972):sig+=1
  assert n==1568972 and valid>0,(t['trait'],n,valid)
  summary.append(dict(trait=t['trait'],species=t['species'],variants=n,valid=valid,failed=n-valid,minP=mn,suggestive=nom,studywide=sig))
 with (L/'summary.tsv').open('w') as f:
  w=csv.DictWriter(f,list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
 (L/'suggestive_hits.json').write_text(json.dumps(hits))
 status('COMPLETE',traits=190,samples=879)
 (L/'VALIDATED.json').write_text(json.dumps(dict(traits=190,samples=879,variants_per_trait=1568972)))
try:main()
except BaseException as e:status('FAILED',error=str(e));traceback.print_exc();raise
