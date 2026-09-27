# Running order

Run from the repository root. Prepare input files and any required output directories before starting. Input and output data paths retain the original analysis layout; script paths follow the folders below.

1. **Core taxa:** `scripts/01_core_taxa/01_select_core_taxa.R` uses abundance and metadata for 879 samples. It stops if its QA directory already exists.
2. **Host genetics:** set `CRC_GENETICS_ROOT` to the Linux data directory containing `inputs/`, `genotypes/`, `mgwas/` and `software/`. Run `scripts/02_host_genetics/01_mgwas.py`, then `scripts/02_host_genetics/02_clump_and_cohort_effects.py`. Clumping requires the validated mGWAS results and their `VALIDATED.json` checkpoint. Use `scripts/02_host_genetics/03_heritability.py` for heritability. Pass the organized genetics data root as the first argument to `scripts/02_host_genetics/04_disease_stratified.R`.
3. **Configuration and metabolomics:** `scripts/03_configuration_metabolomics/01_configuration_metabolomics.R` uses the 57-taxon selection, all-190 trait annotations, abundance, metadata and metabolite tables.
4. **Mediation:** `scripts/04_mediation_mr/01_mediation_bootstrap.R` loads `scripts/04_mediation_mr/mediation_calc.R` and uses six Unix workers. Supply the matched dosage, phenotype, PC and covariate tables.
5. **MR:** `scripts/04_mediation_mr/02_mendelian_randomization.R` requires normalized outcome statistics, lead/proxy exposures, build mappings and LD tables.
6. **Pfam:** run `scripts/05_functional_analysis/01_screen_pfam.R`, `scripts/05_functional_analysis/02_pfam_enrichment.py`, `scripts/05_functional_analysis/03_map_pfam_genes.py`, `scripts/05_functional_analysis/04_aggregate_pfam_tpm.py`, then `scripts/05_functional_analysis/05_community_associations.R`. Enrichment requires `--manifest`, `--pfam-hits`, `--feature-key` and `--outdir`; aggregation requires `--mapping`, `--tpm` and `--outdir`. Place enrichment tables under `Figure5/Tables/Fig5b`, mapping summaries under `Linked_Fig5CD_Analysis/Tables`, and aggregation outputs under `Linked_Fig5CD_Analysis/Tables_v2`.
7. **Prediction:** run `scripts/06_prediction/01_repeated_cv.R`, then `scripts/06_prediction/02_oof_shap.R`. For the combined model, run `scripts/06_prediction/03_compare_taxa_pfam.R`, then `scripts/06_prediction/04_paired_bootstrap.R`. Prediction uses 1,441 samples; community association models use 1,108 complete cases.
8. **Genotype strata:** `scripts/06_prediction/05_genotype_stratified.R` uses existing OOF predictions and hard-genotype calls. It evaluates rs13134510 and scans the 67 lead SNPs.

## Dependencies

- R 4.4.2; packages listed in `R_packages.txt`. Methods report glmnet 4.1-10 and MMUPHin 2.0.0.
- Python 3 and NumPy (`../requirements.txt`).
- PLINK 2 and GCTA 1.95.1 for host genetics; binaries are expected under `CRC_GENETICS_ROOT/software/`.
- CD-HIT 4.8.1 for protein clustering; Prodigal and Pfam annotations are prepared inputs.

Host-genetic Python scripts require Linux (`/proc/meminfo`, `fcntl` and process signals). Their resource gates can require more than 96 GiB of available RAM. Review these settings before running on a different server.

The package list is not a version-locked environment.

## Calculation notes

The current scripts retain two issues identified during code review: the MR input-directory variable `i` is reused as a loop index, and `ap_fast` does not group tied scores. The latter slightly changes AP bootstrap intervals and P values; ROC AUC is unaffected. Repeated-CV t intervals describe split variability, not population sampling uncertainty.
