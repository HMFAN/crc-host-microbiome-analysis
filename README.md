# CRC host–microbiome analysis

Analysis scripts for host genetics, gut microbial features and colorectal cancer.

| Folder | Contents |
|---|---|
| [01_core_taxa](scripts/01_core_taxa) | Core taxa selection |
| [02_host_genetics](scripts/02_host_genetics) | mGWAS, LD clumping, heritability and stratified associations |
| [03_configuration_metabolomics](scripts/03_configuration_metabolomics) | Microbial configuration and metabolomic associations |
| [04_mediation_mr](scripts/04_mediation_mr) | Mediation and Mendelian randomization |
| [05_functional_analysis](scripts/05_functional_analysis) | Pfam enrichment and community functional profiles |
| [06_prediction](scripts/06_prediction) | Cross-validation, SHAP, paired bootstrap and genotype strata |

The repository contains 20 analysis and helper scripts. Study data and plotting code are not included.

## Running the analyses

Use R 4.4.2 and Python 3. Host-genetic analyses also require PLINK 2 and GCTA on Linux.

Prepare the [inputs](docs/INPUTS.md), configure data paths, and run scripts from the repository root. For example:

```sh
Rscript scripts/06_prediction/01_repeated_cv.R
Rscript scripts/06_prediction/02_oof_shap.R
```

See [running order, dependencies and calculation notes](docs/RUNNING.md). Some inputs are intermediate results prepared outside this repository.
