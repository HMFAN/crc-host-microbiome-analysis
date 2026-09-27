#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- "."
infile <- file.path(root, "02_Analysis/Bacteroides_HQ_wholegenome_Pfam/Server_Results/tables/Pfam_all_genomes_long.tsv")
outdir <- file.path(root, "03_Results_9/Figure5/Unbiased_Pfam_Screen")
tabdir <- file.path(outdir, "Tables")
qadir <- file.path(outdir, "QA")
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(qadir, recursive = TRUE, showWarnings = FALSE)

bc_id <- "ERR1293731.10"
dt <- fread(infile)
stopifnot(uniqueN(dt$genome_id) == 20L, uniqueN(dt$feature) == 3668L)
stopifnot(nrow(dt) == 20L * 3668L)

dt[, rank_min := frank(-proteins_per_1000, ties.method = "min"), by = feature]
dt[, rank_average := frank(-proteins_per_1000, ties.method = "average"), by = feature]

bc <- dt[genome_id == bc_id, .(
  feature, accession, description,
  Bc_proteins = proteins_with_domain,
  Bc_density_per_1000 = proteins_per_1000,
  Bc_rank_min = rank_min,
  Bc_rank_average = rank_average
)]

others <- dt[genome_id != bc_id, .(
  other_mean_density = mean(proteins_per_1000),
  other_median_density = median(proteins_per_1000),
  other_sd_density = sd(proteins_per_1000),
  other_max_density = max(proteins_per_1000),
  other_nonzero_n = sum(proteins_with_domain > 0),
  other_total_proteins = sum(proteins_with_domain)
), by = feature]

all_prev <- dt[, .(
  genome_prevalence_of_20 = sum(proteins_with_domain > 0),
  total_proteins_all20 = sum(proteins_with_domain)
), by = feature]

res <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), list(bc, others, all_prev))
res[, Bc_to_other_mean_ratio := fifelse(other_mean_density == 0,
                                         fifelse(Bc_density_per_1000 > 0, Inf, NA_real_),
                                         Bc_density_per_1000 / other_mean_density)]
res[, Bc_minus_other_mean := Bc_density_per_1000 - other_mean_density]
res[, log2_ratio_pc01 := log2((Bc_density_per_1000 + 0.1) / (other_mean_density + 0.1))]
res[, feature_scope := fifelse(genome_prevalence_of_20 >= 10,
                               "broadly distributed (>=10/20)",
                               "restricted (<10/20)")]
res[, primary_core := Bc_proteins >= 3 & Bc_rank_min == 1 & Bc_to_other_mean_ratio >= 2]
res[, extended := Bc_proteins >= 3 & Bc_rank_min <= 3 & Bc_to_other_mean_ratio >= 1.5]
res[, evidence_tier := fcase(
  primary_core & genome_prevalence_of_20 >= 10, "Core-A: rank 1, >=2x, broadly distributed",
  primary_core & genome_prevalence_of_20 < 10, "Core-B: rank 1, >=2x, restricted",
  extended, "Extended: rank <=3, >=1.5x",
  default = "Not shortlisted"
)]

setorder(res, -primary_core, -extended, Bc_rank_min, -Bc_to_other_mean_ratio, -Bc_proteins, feature)
fwrite(res, file.path(tabdir, "Bc_all3668_Pfam_unbiased_metrics.tsv"), sep = "\t", na = "NA")

core <- res[primary_core == TRUE]
setorder(core, feature_scope, -Bc_to_other_mean_ratio, -Bc_proteins, feature)
fwrite(core, file.path(tabdir, "Bc_core_unbiased_Pfam.tsv"), sep = "\t", na = "NA")

ext <- res[extended == TRUE]
setorder(ext, Bc_rank_min, -Bc_to_other_mean_ratio, -Bc_proteins, feature)
fwrite(ext, file.path(tabdir, "Bc_extended_unbiased_Pfam.tsv"), sep = "\t", na = "NA")

qa <- data.table(
  check = c("genomes", "Pfam_features", "matrix_rows", "Bc_rows", "core_features", "extended_features"),
  value = c(uniqueN(dt$genome_id), uniqueN(dt$feature), nrow(dt), nrow(bc), nrow(core), nrow(ext)),
  expected = c(20, 3668, 73360, 3668, NA, NA),
  pass = c(uniqueN(dt$genome_id) == 20, uniqueN(dt$feature) == 3668,
           nrow(dt) == 73360, nrow(bc) == 3668, TRUE, TRUE)
)
fwrite(qa, file.path(qadir, "unbiased_Pfam_screen_QA.tsv"), sep = "\t", na = "NA")

cat(sprintf("All Pfam: %d\nCore: %d (broad %d; restricted %d)\nExtended: %d\n",
            nrow(res), nrow(core),
            core[genome_prevalence_of_20 >= 10, .N],
            core[genome_prevalence_of_20 < 10, .N],
            nrow(ext)))
