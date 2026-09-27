#!/usr/bin/env python3
import argparse, csv, re, time
from collections import defaultdict
from pathlib import Path
import numpy as np

def norm_gene(x):
    return re.sub(r"_[0-9]+$", "", x.strip())

ap = argparse.ArgumentParser()
ap.add_argument("--mapping", required=True)
ap.add_argument("--tpm", required=True)
ap.add_argument("--outdir", required=True)
args = ap.parse_args()
out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

gene_pfams = defaultdict(set)
gene_classes = defaultdict(set)
gene_classes11 = defaultdict(set)
pfam_is_sig11 = {}
with open(args.mapping) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        g = norm_gene(r["gene_id"])
        gene_pfams[g].add(r["matched_pfam"])
        gene_classes[g].add(r["functional_class"])
        pfam_is_sig11[r["matched_pfam"]] = r["in_Fig5B_FDR11"].lower() == "true"
        if pfam_is_sig11[r["matched_pfam"]]:
            gene_classes11[g].add(r["functional_class"])

t0 = time.time(); scanned = matched = 0
pfam_acc = {}; class38_acc = {}; class11_acc = {}
matched_rows = []
with open(args.tpm) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    samples = header[1:]
    n = len(samples)
    for line in fh:
        scanned += 1
        gene, sep, rest = line.partition("\t")
        if gene not in gene_pfams:
            if scanned % 1000000 == 0:
                print(f"scanned={scanned:,} matched={matched:,} min={(time.time()-t0)/60:.1f}", flush=True)
            continue
        vals = np.fromstring(rest, sep="\t", dtype=np.float64)
        if len(vals) != n:
            continue
        matched += 1
        matched_rows.append((gene, ",".join(sorted(gene_pfams[gene])), ",".join(sorted(gene_classes[gene]))))
        for pf in gene_pfams[gene]:
            if pf not in pfam_acc: pfam_acc[pf] = np.zeros(n)
            pfam_acc[pf] += vals
        for cl in gene_classes[gene]:
            if cl not in class38_acc: class38_acc[cl] = np.zeros(n)
            class38_acc[cl] += vals
        for cl in gene_classes11.get(gene, ()):
            if cl not in class11_acc: class11_acc[cl] = np.zeros(n)
            class11_acc[cl] += vals
        if scanned % 1000000 == 0:
            print(f"scanned={scanned:,} matched={matched:,} min={(time.time()-t0)/60:.1f}", flush=True)

def write_matrix(path, label, acc):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([label] + samples)
        for k in sorted(acc): w.writerow([k] + [f"{x:.10g}" for x in acc[k]])

write_matrix(out/"Pfam38_abundance_TPM.tsv", "Pfam", pfam_acc)
write_matrix(out/"FunctionalClass38_abundance_TPM.tsv", "functional_class", class38_acc)
write_matrix(out/"FunctionalClass_FDR11_abundance_TPM.tsv", "functional_class", class11_acc)
with open(out/"matched_gene_rows.tsv", "w", newline="") as fh:
    w=csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["gene_id","Pfams","functional_classes"]); w.writerows(matched_rows)
with open(out/"aggregation_QA.tsv", "w", newline="") as fh:
    w=csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["metric","value"])
    for k,v in [("TPM_gene_rows_scanned",scanned),("target_gene_rows_matched",matched),("samples",len(samples)),("output_Pfam",len(pfam_acc)),("output_classes38",len(class38_acc)),("output_classes11",len(class11_acc)),("elapsed_minutes",(time.time()-t0)/60)]: w.writerow([k,v])
print(f"DONE scanned={scanned:,} matched={matched:,} minutes={(time.time()-t0)/60:.1f}")
