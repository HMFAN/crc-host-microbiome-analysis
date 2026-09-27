#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import subprocess
from collections import defaultdict
from pathlib import Path

BC_ID = "ERR1293731.10"

def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def fasta_records(path: Path):
    name = None
    seq = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)

def fisher_greater(a, b, c, d):
    n = a + b + c + d
    row1 = a + b
    col1 = a + c
    max_x = min(row1, col1)

    def log_hypergeom(x):
        return (
            math.lgamma(col1 + 1) - math.lgamma(x + 1) - math.lgamma(col1 - x + 1)
            + math.lgamma(n - col1 + 1) - math.lgamma(row1 - x + 1)
            - math.lgamma((n - col1) - (row1 - x) + 1)
            - (math.lgamma(n + 1) - math.lgamma(row1 + 1) - math.lgamma(n - row1 + 1))
        )

    logs = [log_hypergeom(x) for x in range(a, max_x + 1)]
    peak = max(logs)
    return min(1.0, math.exp(peak) * sum(math.exp(v - peak) for v in logs))

def bh(rows):
    indexed = sorted(enumerate(rows), key=lambda x: x[1]["P"])
    previous = 1.0
    m = len(rows)
    for reverse_rank, (idx, row) in enumerate(reversed(indexed), 1):
        rank = m - reverse_rank + 1
        q = min(previous, row["P"] * m / rank)
        rows[idx]["FDR"] = q
        previous = q

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--pfam-hits", required=True)
    parser.add_argument("--feature-key", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    pfam_path = Path(args.pfam_hits)
    feature_path = Path(args.feature_key)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = read_tsv(manifest_path)
    assert len(manifest) == 20
    assert len({r["gtdb_species"] for r in manifest}) == 20
    assert sum(r["genome_id"] == BC_ID for r in manifest) == 1
    assert all(float(r["checkm_completeness"]) >= 90 for r in manifest)
    assert all(float(r["checkm_contamination"]) <= 5 for r in manifest)

    combined = outdir / "HQ20_all_proteins.faa"
    gene_to_genome = {}
    observed_counts = defaultdict(int)
    with combined.open("w") as out:
        for row in manifest:
            gid = row["genome_id"]
            faa = Path(row["normalized_faa"])
            if not faa.exists():
                raise FileNotFoundError(faa)
            for protein_id, seq in fasta_records(faa):
                new_id = protein_id if protein_id.startswith(gid + "|") else f"{gid}|{protein_id}"
                if new_id in gene_to_genome:
                    raise RuntimeError(f"Duplicate protein id: {new_id}")
                gene_to_genome[new_id] = gid
                observed_counts[gid] += 1
                out.write(f">{new_id}\n{seq}\n")

    for row in manifest:
        expected = int(row["protein_count"])
        observed = observed_counts[row["genome_id"]]
        if expected != observed:
            raise RuntimeError(f"Protein count mismatch {row['genome_id']}: {observed} != {expected}")

    cluster_prefix = outdir / "HQ20_cdhit_c70_aS80"
    clstr = Path(str(cluster_prefix) + ".clstr")
    if not clstr.exists() or clstr.stat().st_size == 0:
        subprocess.run([
            "cd-hit", "-i", str(combined), "-o", str(cluster_prefix),
            "-c", "0.70", "-aS", "0.80", "-G", "0", "-g", "1",
            "-T", str(args.threads), "-M", "24000", "-d", "0"
        ], check=True)

    clusters = []
    current = None
    with clstr.open() as handle:
        for line in handle:
            if line.startswith(">Cluster "):
                if current is not None:
                    clusters.append(current)
                current = {"cluster_id": line.strip()[1:].replace(" ", "_"), "members": []}
            else:
                token = line.split(">", 1)[1].split("...", 1)[0]
                current["members"].append(token)
        if current is not None:
            clusters.append(current)

    bc_membership = []
    selected = set()
    status_counts = defaultdict(int)
    for cluster in clusters:
        genomes = sorted({gene_to_genome[g] for g in cluster["members"]})
        bc_genes = [g for g in cluster["members"] if gene_to_genome[g] == BC_ID]
        if not bc_genes:
            continue
        n_genomes = len(genomes)
        if n_genomes == 1:
            status = "Bc_specific_HQ20"
        elif n_genomes <= 3:
            status = "Bc_rare_HQ20"
        else:
            status = "Bc_shared_HQ20"
        status_counts[status] += 1
        for gene in bc_genes:
            bc_membership.append({
                "cluster_id": cluster["cluster_id"], "gene_id": gene,
                "status": status, "cluster_n_genomes": n_genomes,
                "present_genomes": ";".join(genomes)
            })
            if n_genomes <= 3:
                selected.add(gene)

    bc_all = {gene for gene, genome in gene_to_genome.items() if genome == BC_ID}
    assert len(bc_all) == int(next(r["protein_count"] for r in manifest if r["genome_id"] == BC_ID))
    assert selected <= bc_all

    features = read_tsv(feature_path)
    assert len(features) == 38 and len({r["feature"] for r in features}) == 38
    feature_to_bc_genes = defaultdict(set)
    for row in read_tsv(pfam_path):
        gene = row["query_gene_id"]
        if gene in bc_all:
            feature_to_bc_genes[row["pfam_name"]].add(gene)

    background = bc_all - selected
    enrichment = []
    for feat in features:
        name = feat["feature"]
        carriers = feature_to_bc_genes.get(name, set())
        a = len(selected & carriers)
        b = len(selected) - a
        c = len(background & carriers)
        d = len(background) - c
        p = fisher_greater(a, b, c, d)
        odds_ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        row = dict(feat)
        row.update({
            "selected_with_feature": a, "selected_total": len(selected),
            "background_with_feature": c, "background_total": len(background),
            "selected_rate": a / len(selected), "background_rate": c / len(background),
            "accessory_fraction": a / len(carriers) if carriers else 0.0,
            "log2_OR_HA": math.log2(odds_ha), "P": p
        })
        enrichment.append(row)
    bh(enrichment)
    enrichment.sort(key=lambda r: (r["FDR"], -r["selected_with_feature"], r["feature"]))

    write_tsv(outdir / "HQ20_Bc_gene_cluster_membership.tsv", bc_membership,
              ["cluster_id", "gene_id", "status", "cluster_n_genomes", "present_genomes"])
    write_tsv(outdir / "HQ20_Bc_specific_rare_genes.tsv",
              [r for r in bc_membership if r["gene_id"] in selected],
              ["cluster_id", "gene_id", "status", "cluster_n_genomes", "present_genomes"])
    write_tsv(outdir / "HQ20_38Pfam_specific_rare_enrichment.tsv", enrichment, list(enrichment[0]))

    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    summary = [
        {"metric": "analysis_label", "value": "HQ20_only_no_legacy33"},
        {"metric": "manifest_sha256", "value": manifest_hash},
        {"metric": "n_genomes", "value": len(manifest)},
        {"metric": "n_species", "value": len({r['gtdb_species'] for r in manifest})},
        {"metric": "n_all_proteins", "value": len(gene_to_genome)},
        {"metric": "n_clusters", "value": len(clusters)},
        {"metric": "n_Bc_proteins", "value": len(bc_all)},
        {"metric": "n_Bc_specific_clusters", "value": status_counts["Bc_specific_HQ20"]},
        {"metric": "n_Bc_rare_clusters", "value": status_counts["Bc_rare_HQ20"]},
        {"metric": "n_Bc_specific_rare_proteins", "value": len(selected)},
        {"metric": "n_Bc_background_proteins", "value": len(background)},
        {"metric": "n_tested_Pfam", "value": len(enrichment)},
        {"metric": "n_FDR_lt_0.05", "value": sum(r["FDR"] < 0.05 for r in enrichment)},
        {"metric": "cdhit_identity", "value": 0.70},
        {"metric": "cdhit_shorter_coverage", "value": 0.80},
    ]
    write_tsv(outdir / "HQ20_analysis_summary.tsv", summary, ["metric", "value"])
    write_tsv(outdir / "HQ20_manifest_used.tsv", manifest, list(manifest[0]))
    print(f"DONE HQ20: clusters={len(clusters)} selected_proteins={len(selected)} FDR={sum(r['FDR'] < 0.05 for r in enrichment)}")

if __name__ == "__main__":
    main()
