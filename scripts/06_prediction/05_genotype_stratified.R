#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(digest)
})

set.seed(20260924)

root <- "03_Results_9/Figure6/PanelC_Fig6a_FrozenOOF_Genotype"
td <- file.path(root, "Tables")
qd <- file.path(root, "QA")
dir.create(td, recursive = TRUE, showWarnings = FALSE)
dir.create(qd, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  fig6a_oof = "03_Results_9/RepeatedCV_mGWAS57/Tables/04_sample_OOF_predictions.tsv",
  hard_genotype = "03_Results_9/Exploratory_Genotype_Subgroup_CRC/Inputs/lead67_hard_GT.raw"
)
stopifnot(all(file.exists(inputs)))

p <- fread(inputs["fig6a_oof"])
g <- fread(inputs["hard_genotype"], na.strings = c("", "NA"))
prob_cols <- paste0("probability_repeat", 1:10)
gt_cols <- setdiff(names(g), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
stopifnot(nrow(p) == 1441L, uniqueN(p$IID) == 1441L,
          all(prob_cols %in% names(p)), length(gt_cols) == 67L,
          uniqueN(g$IID) == nrow(g), sum(p$IID %in% g$IID) == 879L)

d <- merge(p, g[, c("IID", gt_cols), with = FALSE], by = "IID", all = FALSE, sort = FALSE)
stopifnot(nrow(d) == 879L)
d[, population := fifelse(Cohort == "JAP", "JAP",
                           fifelse(Cohort %chin% c("CHN1", "CHN2", "CHN3"),
                                   "Pooled CHN", NA_character_))]

auc_safe <- function(y, pr) {
  if (length(y) < 2L || uniqueN(y) < 2L) return(NA_real_)
  as.numeric(roc(y, pr, quiet = TRUE, direction = "<")$auc)
}

target_snp <- "4:147366360:C:T"
target_col <- grep(paste0("^", gsub(":", "\\:", target_snp), "_[ACGT]+$"), gt_cols, value = TRUE)
stopifnot(length(target_col) == 1L, sub(".*_", "", target_col) == "T")

target <- d[!is.na(population) & !is.na(get(target_col))]
target[, dosage := get(target_col)]
target[, genotype := fifelse(dosage == 0, "Non-carrier", "T-carrier")]
target[, genotype := factor(genotype, levels = c("Non-carrier", "T-carrier"))]
target[, population := factor(population, levels = c("Pooled CHN", "JAP"))]
stopifnot(nrow(target[population == "Pooled CHN"]) == 174L,
          nrow(target[population == "JAP"]) == 83L)

target_qc <- target[, .(
  n = .N, CRC = sum(CRC == 1), control = sum(CRC == 0),
  CRC_prevalence = mean(CRC)
), by = .(population, genotype)]

target_auc <- target[, .(
  ensemble_AUC = auc_safe(CRC, mean_OOF_probability),
  n = .N, CRC = sum(CRC == 1), control = sum(CRC == 0)
), by = .(population, genotype)]

target_long <- melt(
  target,
  id.vars = c("IID", "CRC", "population", "genotype"),
  measure.vars = prob_cols,
  variable.name = "repeat_label", value.name = "OOF_probability"
)
target_long[, repetition := as.integer(sub("probability_repeat", "", repeat_label))]
target_repeat_auc <- target_long[, .(
  AUC = auc_safe(CRC, OOF_probability)
), by = .(population, genotype, repetition)]
target_repeat_wide <- dcast(target_repeat_auc, population + repetition ~ genotype,
                            value.var = "AUC")
target_repeat_wide[, delta_AUC := `T-carrier` - `Non-carrier`]

target_repeat_summary <- target_repeat_wide[, .(
  mean_delta_AUC = mean(delta_AUC), SD_delta_AUC = sd(delta_AUC),
  median_delta_AUC = median(delta_AUC), min_delta_AUC = min(delta_AUC),
  max_delta_AUC = max(delta_AUC), positive_repeats = sum(delta_AUC > 0),
  total_repeats = .N
), by = population]

roc_coords <- rbindlist(lapply(seq_len(nrow(target_auc)), function(i) {
  pop <- target_auc$population[i]
  gen <- target_auc$genotype[i]
  z <- target[population == pop & genotype == gen]
  rr <- roc(z$CRC, z$mean_OOF_probability, quiet = TRUE, direction = "<")
  data.table(population = pop, genotype = gen,
             FPR = 1 - rr$specificities, TPR = rr$sensitivities)
}))

B <- 2000L
target_boot <- rbindlist(lapply(levels(target$population), function(pop) {
  z <- target[population == pop]
  strata <- split(seq_len(nrow(z)), interaction(z$genotype, z$CRC, drop = TRUE))
  vals <- numeric(B)
  for (bb in seq_len(B)) {
    ii <- unlist(lapply(strata, function(j) sample(j, length(j), replace = TRUE)),
                 use.names = FALSE)
    zz <- z[ii]
    vals[bb] <- auc_safe(zz[genotype == "T-carrier"]$CRC,
                         zz[genotype == "T-carrier"]$mean_OOF_probability) -
      auc_safe(zz[genotype == "Non-carrier"]$CRC,
               zz[genotype == "Non-carrier"]$mean_OOF_probability)
  }
  data.table(population = pop, bootstrap = seq_len(B), delta_AUC = vals)
}))

target_delta <- rbindlist(lapply(levels(target$population), function(pop) {
  a <- target_auc[population == pop]
  observed <- a[genotype == "T-carrier"]$ensemble_AUC -
    a[genotype == "Non-carrier"]$ensemble_AUC
  bb <- target_boot[population == pop]$delta_AUC
  data.table(population = pop, delta_AUC = observed,
             CI_low = unname(quantile(bb, 0.025)),
             CI_high = unname(quantile(bb, 0.975)),
             bootstrap_replicates = B)
}))

scan_one <- function(gt_col) {
  snp <- sub("_[ACGT]+$", "", gt_col)
  effect_allele <- sub(".*_", "", gt_col)
  z0 <- d[!is.na(population) & !is.na(get(gt_col)),
          .(IID, CRC, population, dosage = get(gt_col), mean_OOF_probability,
            probability_repeat1, probability_repeat2, probability_repeat3,
            probability_repeat4, probability_repeat5, probability_repeat6,
            probability_repeat7, probability_repeat8, probability_repeat9,
            probability_repeat10)]
  if (!nrow(z0)) return(NULL)
  z0[, genotype := fifelse(dosage == 0, "Non-carrier", "T-carrier")]

  out <- list()
  for (pop in c("Pooled CHN", "JAP")) {
    z <- z0[population == pop]
    qc <- z[, .(n = .N, CRC = sum(CRC == 1), control = sum(CRC == 0)), by = genotype]
    getn <- function(grp, col) {
      v <- qc[genotype == grp, get(col)]
      if (length(v)) v else 0L
    }
    anc <- auc_safe(z[genotype == "Non-carrier"]$CRC,
                    z[genotype == "Non-carrier"]$mean_OOF_probability)
    acar <- auc_safe(z[genotype == "T-carrier"]$CRC,
                     z[genotype == "T-carrier"]$mean_OOF_probability)

    repeat_delta <- numeric(10)
    for (rr in 1:10) {
      pc <- paste0("probability_repeat", rr)
      repeat_delta[rr] <- auc_safe(z[genotype == "T-carrier"]$CRC,
                                   z[genotype == "T-carrier"][[pc]]) -
        auc_safe(z[genotype == "Non-carrier"]$CRC,
                 z[genotype == "Non-carrier"][[pc]])
    }
    min_cell <- min(getn("Non-carrier", "CRC"), getn("Non-carrier", "control"),
                    getn("T-carrier", "CRC"), getn("T-carrier", "control"))
    out[[pop]] <- data.table(
      SNP = snp, effect_allele = effect_allele, population = pop,
      n_noncarrier = getn("Non-carrier", "n"),
      CRC_noncarrier = getn("Non-carrier", "CRC"),
      control_noncarrier = getn("Non-carrier", "control"),
      n_carrier = getn("T-carrier", "n"),
      CRC_carrier = getn("T-carrier", "CRC"),
      control_carrier = getn("T-carrier", "control"),
      min_genotype_status_cell = min_cell,
      AUC_noncarrier = anc, AUC_carrier = acar,
      delta_AUC = acar - anc,
      mean_repeat_delta_AUC = mean(repeat_delta, na.rm = TRUE),
      SD_repeat_delta_AUC = sd(repeat_delta, na.rm = TRUE),
      positive_repeats = sum(repeat_delta > 0, na.rm = TRUE),
      valid_repeats = sum(is.finite(repeat_delta))
    )
  }
  rbindlist(out, fill = TRUE)
}

scan67 <- rbindlist(lapply(gt_cols, scan_one), fill = TRUE)
scan67[, stable_positive := valid_repeats == 10L & positive_repeats >= 8L]
scan67[, sufficient_for_screen := min_genotype_status_cell >= 10L]
scan67[, similar_high_discrimination_pattern :=
         sufficient_for_screen & AUC_carrier >= 0.80 & delta_AUC >= 0.10 & stable_positive]
setorder(scan67, population, -delta_AUC)

pooled_rankable <- scan67[population == "Pooled CHN" & sufficient_for_screen &
                            is.finite(delta_AUC)]
pooled_rankable[, delta_rank := frank(-delta_AUC, ties.method = "min")]
target_rank <- pooled_rankable[SNP == target_snp]
top_candidates <- pooled_rankable[order(delta_rank)][1:min(.N, 20)]

fwrite(target_qc, file.path(td, "01_target_sample_QC.tsv"), sep = "\t")
fwrite(target_auc, file.path(td, "02_target_ensemble_AUC.tsv"), sep = "\t")
fwrite(target_repeat_auc, file.path(td, "03_target_AUC_by_repeat.tsv"), sep = "\t")
fwrite(target_repeat_summary, file.path(td, "04_target_repeat_delta_summary.tsv"), sep = "\t")
fwrite(roc_coords, file.path(td, "05_target_ROC_coordinates.tsv"), sep = "\t")
fwrite(target_delta, file.path(td, "06_target_delta_AUC_bootstrap_summary.tsv"), sep = "\t")
fwrite(target_boot, file.path(td, "07_target_delta_AUC_bootstrap_replicates.tsv"), sep = "\t")
fwrite(scan67, file.path(td, "08_all67_frozen_Fig6a_genotype_scan.tsv"), sep = "\t")
fwrite(top_candidates, file.path(td, "09_top20_pooled_CHN_candidates.tsv"), sep = "\t")
fwrite(target_rank, file.path(td, "10_target_SNP_rank.tsv"), sep = "\t")

fwrite(data.table(
  item = c("prediction_source", "sample_universe", "repeats", "model_refit",
           "target_SNP", "target_effect_allele",
           "scan_candidate_rule"),
  value = c(
    inputs["fig6a_oof"],
    "Fig6a 1,441-sample OOF predictions intersected with 879 hard-genotyped samples",
    "10 repeated five-fold CV OOF predictions from Fig6a",
    "none",
    target_snp, "T",
    "pooled CHN: minimum genotype-by-status cell n>=10, carrier AUC>=0.80, delta AUC>=0.10, positive delta in >=8/10 repeats"
  )
), file.path(qd, "run_summary.tsv"), sep = "\t")

fwrite(data.table(role = names(inputs), path = unname(inputs),
                  sha256 = vapply(inputs, digest, "", file = TRUE, algo = "sha256")),
       file.path(qd, "input_manifest.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(qd, "sessionInfo.txt"))

cat("\nTarget sample QC\n")
print(target_qc)
cat("\nTarget ensemble AUC\n")
print(target_auc)
cat("\nTarget delta AUC bootstrap\n")
print(target_delta)
cat("\nTarget rank among pooled CHN SNPs with min genotype-status cell n>=10\n")
print(target_rank)
cat("\nOther similar pooled CHN candidates\n")
print(scan67[population == "Pooled CHN" & similar_high_discrimination_pattern & SNP != target_snp][order(-delta_AUC)])
