suppressPackageStartupMessages({
  library(data.table)
  library(lmerTest)
  library(lme4)
})

root <- "03_Results_9/Figure5/Linked_Fig5CD_Analysis"
tab <- file.path(root, "Tables_v2")
dir.create(tab, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  metadata = "03_Results_9/Figure5/Tables/Inputs/Model_samples.tsv",
  bc = "02_Analysis/BioSampleDedup_Reanalysis_20260906_v1/02_Core_rSGB_v2/01_all_rSGB_relative_abundance.tsv",
  pfam = file.path(tab, "Pfam38_abundance_TPM.tsv"),
  class38 = file.path(tab, "FunctionalClass38_abundance_TPM.tsv"),
  class11 = file.path(tab, "FunctionalClass_FDR11_abundance_TPM.tsv")
)
stopifnot(all(file.exists(paths)))

bc_feature <- paste0(
  "k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|",
  "f__Bacteroidaceae|g__Bacteroides|s__ERR1293731.10:Bacteroides_caccae"
)

m <- fread(paths[["metadata"]])
stopifnot(nrow(m) == 1108L, uniqueN(m$IID) == 1108L, uniqueN(m$BioSample) == 1108L)
m[, `:=`(
  disease_status = factor(condition, levels = c("CTR", "CRC")),
  sex = factor(Sex), Cohort = factor(Cohort),
  Age_z = as.numeric(scale(Age)), BMI_z = as.numeric(scale(BMI))
)]

a <- fread(paths[["bc"]])
br <- a[feature == bc_feature]
stopifnot(nrow(br) == 1L)
bc <- data.table(IID = names(br)[-1], Bc_relative_abundance = as.numeric(br[1, -1]))
m <- merge(m, bc, by = "IID", all.x = TRUE)
stopifnot(!anyNA(m$Bc_relative_abundance))
m[, Bc_z := as.numeric(scale(log10(Bc_relative_abundance + 1e-8)))]

read_module <- function(path, set_name) {
  z <- fread(path)
  setnames(z, 1L, "feature")
  long <- melt(z, id.vars = "feature", variable.name = "IID", value.name = "TPM")
  long[, IID := sub("\\.quant$", "", IID)]
  long[, set := set_name]
  long
}

long <- rbindlist(list(
  read_module(paths[["pfam"]], "Individual Pfam (Fig5a 38)"),
  read_module(paths[["class38"]], "Functional class (Fig5a 38)"),
  read_module(paths[["class11"]], "Functional class (Fig5b FDR11)")
), use.names = TRUE)
long <- merge(long, m, by = "IID")
long[, feature_log10 := log10(TPM + 1)]
long[, feature_z := as.numeric(scale(feature_log10)), by = .(set, feature)]

coef_term <- function(fit, term) {
  cf <- summary(fit)$coefficients
  if (!term %in% rownames(cf)) return(data.table(beta=NA_real_,SE=NA_real_,df=NA_real_,P=NA_real_))
  x <- cf[term,]
  data.table(beta=unname(x["Estimate"]), SE=unname(x["Std. Error"]),
             df=unname(x["df"]), P=unname(x["Pr(>|t|)"]))
}

res <- list(); i <- 0L
for (ss in unique(long$set)) for (ff in unique(long[set == ss]$feature)) {
  d <- long[set == ss & feature == ff]
  f0 <- lmer(feature_z ~ disease_status + Age_z + sex + BMI_z + (1|Cohort), data=d, REML=FALSE)
  f1 <- lmer(feature_z ~ disease_status + Bc_z + Age_z + sex + BMI_z + (1|Cohort), data=d, REML=FALSE)
  fb <- lmer(feature_z ~ Bc_z + disease_status + Age_z + sex + BMI_z + (1|Cohort), data=d, REML=FALSE)
  for (spec in c("CRC unadjusted for Bc","CRC adjusted for Bc","Association with Bc")) {
    rr <- switch(spec,
      "CRC unadjusted for Bc" = coef_term(f0,"disease_statusCRC"),
      "CRC adjusted for Bc" = coef_term(f1,"disease_statusCRC"),
      "Association with Bc" = coef_term(fb,"Bc_z"))
    i <- i + 1L
    rr[, `:=`(set=ss, feature=ff, model=spec, n=nrow(d),
              prevalence=sum(d$TPM>0)/nrow(d), median_TPM=median(d$TPM),
              singular=if (spec=="CRC unadjusted for Bc") isSingular(f0) else if (spec=="CRC adjusted for Bc") isSingular(f1) else isSingular(fb))]
    res[[i]] <- rr
  }
}
res <- rbindlist(res)
res[, `:=`(CI_low=beta-qt(.975,df)*SE, CI_high=beta+qt(.975,df)*SE)]
res[, FDR := p.adjust(P,"BH"), by=.(set,model)]
setcolorder(res, c("set","feature","model","beta","SE","df","P","FDR","CI_low","CI_high","n","prevalence","median_TPM","singular"))
fwrite(res, file.path(tab,"CRC_and_Bc_adjusted_models.tsv"), sep="\t")

wide <- dcast(res[model %chin% c("CRC unadjusted for Bc","CRC adjusted for Bc")],
              set+feature ~ model, value.var=c("beta","P","FDR"))
wide[, beta_retained_fraction := `beta_CRC adjusted for Bc` / `beta_CRC unadjusted for Bc`]
fwrite(wide, file.path(tab,"CRC_effect_before_after_Bc_adjustment.tsv"), sep="\t")

qa <- data.table(
  metric=c("model_samples","unique_IID","input_sets","input_features","missing_TPM","missing_Bc"),
  value=c(nrow(m),uniqueN(m$IID),uniqueN(long$set),uniqueN(long[,paste(set,feature)]),sum(is.na(long$TPM)),sum(is.na(long$Bc_z)))
)
fwrite(qa,file.path(tab,"model_QA.tsv"),sep="\t")
print(res[set != "Individual Pfam (Fig5a 38)"][order(set,model,FDR)])
