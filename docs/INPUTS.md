# Inputs

Study data are not included. Input and output paths inside the scripts retain the original analysis layout; configure them for your data before running.

| Input family | Required content |
|---|---|
| Prediction metadata | 1,441 unique rows with IID, BioSample, Cohort, condition; 715 CRC and 726 CTR across 11 cohorts |
| Abundance matrix | Feature-by-sample adjusted abundance; sample columns match IID; finite, nonnegative values; separate unadjusted tables where explicitly read |
| Genetic trait definitions | All-190 trait/feature/taxonomy map; frozen retained 57 trait IDs and retained 71 SNP–trait associations |
| Genetics | 879-row metadata/phenotypes/PCs, 1,568,972-SNP prepared dosage/hard-call data and GRM; effect-allele dosage and hard calls are distinct inputs |
| Metabolomics | 172 matched Japanese participants, 421 standardized metabolites, covariates and biochemical class membership |
| Mediation | 879 matched samples, PCs, microbial phenotypes and effect-allele dosages; complete covariates identify the 605-sample adjusted model |
| MR | Current lead exposures, proxy exposures, allele/build records, normalized summary statistics for GCST012879, GCST90129505 and FinnGen DF13, and LD information |
| Comparative genomics | HQ20 genome/protein manifest, Pfam hits, normalized all-domain counts and protein FASTA inputs for clustering |
| Community Pfam | Gene annotation catalogue, gene TPM matrix, mapped Pfam/classes and model sample metadata; 1,108 complete observations for community models |
| Prediction Pfam | Pfam TPM and mapping summary for the full 1,441-sample prediction universe; do not replace with the 1,108-row model subset |

Preserve column names, factor meanings, effect-allele orientation, sample order, sample universes and taxon identifiers. The scripts contain many explicit checks for these properties. Do not weaken a sample-count assertion just to make a different dataset run.

