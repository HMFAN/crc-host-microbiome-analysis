#!/usr/bin/env python3
import csv
from pathlib import Path

ROOT = Path("03_Results_9/Figure5")
OUT = ROOT / "Linked_Fig5CD_Analysis"
TAB = OUT / "Tables"
TAB.mkdir(parents=True, exist_ok=True)

all38_path = ROOT / "Tables/Fig5b/Fig5B_HQ20_all38_Pfam_accessory_enrichment.tsv"
sig11_path = ROOT / "Tables/Fig5b/Fig5B_HQ20_FDR_significant_Pfam.tsv"
anno_path = Path("01_Data/clusterRe95_rep_seq.emapper.annotations")

def read_targets(path):
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["feature"]: r for r in rows}

all38 = read_targets(all38_path)
sig11 = read_targets(sig11_path)
targets = set(all38)

out_path = TAB / "gene_catalog_genes_with_Fig5A38_Pfam.tsv"
counts = {x: 0 for x in targets}
genes = set()
rows_written = 0

with anno_path.open(errors="replace") as src, out_path.open("w", newline="") as dst:
    writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene_id", "matched_pfam", "in_Fig5B_FDR11", "functional_class",
                     "KEGG_ko", "description", "preferred_name", "all_PFAMs"])
    header = None
    for line in src:
        if line.startswith("##"):
            continue
        if line.startswith("#query"):
            header = line[1:].rstrip("\n").split("\t")
            idx = {v: i for i, v in enumerate(header)}
            continue
        if not line.strip() or header is None:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < len(header):
            continue
        pfams = set(x for x in f[idx["PFAMs"]].split(",") if x and x != "-")
        hit = sorted(pfams & targets)
        if not hit:
            continue
        gene = f[idx["query"]]
        genes.add(gene)
        for pf in hit:
            r = all38[pf]
            writer.writerow([
                gene, pf, pf in sig11, r["functional_class"], f[idx["KEGG_ko"]],
                f[idx["Description"]], f[idx["Preferred_name"]], f[idx["PFAMs"]]
            ])
            counts[pf] += 1
            rows_written += 1

with (TAB / "Pfam_gene_catalog_mapping_summary.tsv").open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["pfam", "in_Fig5B_FDR11", "functional_class", "n_gene_catalog_genes"])
    for pf, r in sorted(all38.items(), key=lambda kv: int(kv[1]["feature_order"])):
        writer.writerow([pf, pf in sig11, r["functional_class"], counts[pf]])

with (TAB / "extraction_QA.tsv").open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value"])
    writer.writerow(["target_pfam", len(targets)])
    writer.writerow(["fdr11_pfam", len(sig11)])
    writer.writerow(["unique_mapped_genes", len(genes)])
    writer.writerow(["gene_pfam_rows", rows_written])

print(f"mapped {len(genes)} genes; wrote {rows_written} gene-Pfam rows")
