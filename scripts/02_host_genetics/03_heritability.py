#!/usr/bin/python3
import csv, fcntl, hashlib, json, math, os, signal, subprocess, time, traceback
from pathlib import Path

H = Path(os.environ['CRC_GENETICS_ROOT'])
I, G = H/'inputs', H/'genotypes'
R = H/'heritability_PC1_10_CRC_cohort'

def save(state, **extra):
    d = dict(state=state, pid=os.getpid(), updated=time.strftime('%F %T %z'), **extra)
    p = R/'status.tmp'
    p.write_text(json.dumps(d, indent=2))
    p.replace(R/'status.json')
    print(json.dumps(d), flush=True)

def table(path, header, rows):
    with path.open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        if header:
            w.writerow(header)
        w.writerows(rows)

def rank(columns):
    basis = []
    for col in columns:
        norm = math.sqrt(sum(v*v for v in col))
        assert norm > 0
        x = [v/norm for v in col]
        for _ in range(2):
            for b in basis:
                dot = sum(v*w for v, w in zip(x, b))
                x = [v-dot*w for v, w in zip(x, b)]
        norm = math.sqrt(sum(v*v for v in x))
        if norm > 1e-9:
            basis.append([v/norm for v in x])
    return len(basis)

def resources():
    m = {x.split(':')[0]: int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines()}
    return os.getloadavg()[0], m['MemAvailable']/1048576

def main():
    assert not R.exists(), f'Refusing to overwrite existing output: {R}'
    R.mkdir()
    lock = (R/'controller.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    hashes = json.loads((I/'phenotype_sha256.json').read_text())
    for name, digest in hashes.items():
        assert hashlib.sha256((I/name).read_bytes()).hexdigest() == digest, name
    meta = list(csv.DictReader((I/'metadata.tsv').open(), delimiter='\t'))
    ids = [x['IID'] for x in meta]
    assert len(ids) == len(set(ids)) == 879
    assert all(x['condition'] in ['CRC', 'CTR'] and x['Cohort'] for x in meta)
    traits = list(csv.DictReader((I/'trait_map.tsv').open(), delimiter='\t'))
    assert len(traits) == 190 and len({x['trait'] for x in traits}) == 190
    bc = [x for x in traits if x['genome_id'] == 'ERR1293731.10']
    assert len(bc) == 1
    ph = [x.split() for x in (I/'phenotypes_gcta.tsv').read_text().splitlines()]
    assert [x[1] for x in ph] == ids and all(len(x) == 192 and x[0] == x[1] for x in ph)
    assert all(math.isfinite(float(v)) for x in ph for v in x[2:])

    ev = [x.split() for x in (G/'DS_PC1_10.eigenvec').read_text().splitlines()]
    header = ev[0]
    iid_col = header.index('IID') if 'IID' in header else header.index('#IID')
    pcs = {x[iid_col]: [float(v) for v in x[-10:]] for x in ev[1:]}
    assert set(pcs) == set(ids)
    grm_ids = [x.split()[1] for x in (G/'known_GT_GRM.grm.id').read_text().splitlines()]
    assert grm_ids == ids and len(grm_ids) == 879

    cohorts = {v: i+1 for i, v in enumerate(sorted({x['Cohort'] for x in meta}))}
    categorical = [[int(x['condition'] == 'CRC'), cohorts[x['Cohort']]] for x in meta]
    quantitative = [pcs[x['IID']] for x in meta]
    columns = [[1.0]*879] + [list(c) for c in zip(*quantitative)]
    for col in zip(*categorical):
        levels = sorted(set(col))
        columns += [[float(v == level) for v in col] for level in levels[1:]]
    design_rank = rank(columns)
    assert design_rank == len(columns), (design_rank, len(columns))

    table(R/'all879.keep', None, [[x, x] for x in ids])
    table(R/'categorical.covar', None, [[x, x]+v for x, v in zip(ids, categorical)])
    table(R/'quantitative.qcovar', None, [[x, x]+v for x, v in zip(ids, quantitative)])
    table(R/'covariate_audit.tsv', ['IID', 'CRC_status', 'cohort_code']+[f'PC{i}' for i in range(1,11)],
          [[x]+a+b for x, a, b in zip(ids, categorical, quantitative)])
    preflight = dict(
        traits=190, samples=879, phenotype='current CLR then RINT',
        fixed_effects='PC1-10 + CRC_status + cohort', omitted_covariates='age, BMI, sex',
        covariate_rank=design_rank, covariate_columns=len(columns), cohort_encoding=cohorts,
        GRM='GT known-union 1,568,972 SNPs', PCA='DS PC1-10', FDR_family=190,
        Bc=bc[0], input_hashes=hashes)
    (R/'preflight.json').write_text(json.dumps(preflight, indent=2))
    (R/'results').mkdir()
    (R/'logs').mkdir()

    queue = bc + [x for x in traits if x not in bc]
    active, done, failed = [], [], []
    paused, last = False, None
    os.environ.update(OMP_NUM_THREADS='2', MKL_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
    while queue or active:
        for item in active[:]:
            trait, proc, handle = item
            if proc.poll() is not None:
                handle.close()
                active.remove(item)
                (done if proc.returncode == 0 else failed).append(trait['trait'])
        load, mem = resources()
        if not paused and (load > 104 or mem < 48):
            for _, proc, _ in active:
                try: os.killpg(proc.pid, signal.SIGSTOP)
                except ProcessLookupError: pass
            paused = True
        elif paused and load < 80 and mem > 96:
            for _, proc, _ in active:
                try: os.killpg(proc.pid, signal.SIGCONT)
                except ProcessLookupError: pass
            paused = False
        while not paused and queue and len(active) < 4 and load < 85 and mem > 96:
            trait = queue.pop(0)
            out = R/'results'/trait['trait']
            args = [str(H/'software/gcta64'), '--grm', str(G/'known_GT_GRM'),
                    '--keep', str(R/'all879.keep'), '--pheno', str(I/'phenotypes_gcta.tsv'),
                    '--mpheno', trait['mpheno'], '--covar', str(R/'categorical.covar'),
                    '--qcovar', str(R/'quantitative.qcovar'), '--reml', '--thread-num', '2',
                    '--out', str(out)]
            (R/'logs'/(trait['trait']+'.command.json')).write_text(json.dumps(args))
            handle = (R/'logs'/(trait['trait']+'.stdout')).open('w')
            proc = subprocess.Popen(args, stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
            active.append((trait, proc, handle))
        snap = (len(done), tuple(failed), tuple(x['trait'] for x, _, _ in active), paused)
        if snap != last:
            save('PAUSED_RESOURCES' if paused else ('RUNNING_REML' if active else 'WAITING_RESOURCES'),
                 completed=len(done), failed=failed, active=[x['trait'] for x, _, _ in active],
                 queued=len(queue), load=load, available_GiB=mem)
            last = snap
        if queue or active:
            time.sleep(10)

    rows = []
    for trait in traits:
        data = {}
        try:
            for line in (R/'results'/(trait['trait']+'.hsq')).read_text().splitlines():
                z = line.split()
                if len(z) >= 2: data[z[0]] = z[1:]
            h2, se, pval = float(data['V(G)/Vp'][0]), float(data['V(G)/Vp'][1]), float(data['Pval'][0])
            n = int(data['n'][0])
            assert all(math.isfinite(x) for x in [h2, se, pval]) and 0 <= pval <= 1 and se >= 0 and n == 879
            assert trait['trait'] not in failed
            rows.append(dict(trait, h2=h2, SE=se, P=pval, FDR='NA', n=n, status='OK'))
        except Exception as exc:
            rows.append(dict(trait, h2='NA', SE='NA', P='NA', FDR='NA', n='NA', status='FAILED: '+str(exc)))
    if all(x['status'] == 'OK' for x in rows):
        order = sorted(range(190), key=lambda i: rows[i]['P'])
        previous = 1.0
        for k in range(189, -1, -1):
            j = order[k]
            previous = min(previous, rows[j]['P']*190/(k+1))
            rows[j]['FDR'] = previous
    fields = list(traits[0]) + ['h2', 'SE', 'P', 'FDR', 'n', 'status']
    table(R/'all190_REML.tsv', fields, [[x[k] for k in fields] for x in rows])
    table(R/'Bc_REML.tsv', fields, [[x[k] for k in fields] for x in rows if x['genome_id'] == 'ERR1293731.10'])
    if not all(x['status'] == 'OK' for x in rows):
        raise RuntimeError('At least one REML failed; FDR withheld. See all190_REML.tsv.')
    (R/'VALIDATED.json').write_text(json.dumps(dict(status='COMPLETE', traits=190, n=879,
        FDR_family=190, fixed_effects='PC1-10 + CRC_status + cohort', completed=time.strftime('%F %T %z')), indent=2))
    save('COMPLETE', completed=190, n=879)

if __name__ == '__main__':
    try:
        main()
    except BaseException as exc:
        if R.exists(): save('FAILED', error=str(exc))
        traceback.print_exc()
        raise
