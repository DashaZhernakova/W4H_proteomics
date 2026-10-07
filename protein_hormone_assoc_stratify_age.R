my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/utility_functions.R")

library(ggplot2)
library(ggExtra)

set.seed(123)

out_basedir <- "results12/intensity_all_prots_220526/"

# READ PROTEINS
d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.phase_avg.txt", as.is = T, check.names = F, sep = "\t")

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$phase <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    dplyr::select(SampleID, ID, phase, everything())
}

d_wide$phase <- relevel(factor(d_wide$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
#d_wide$TP <- as.numeric(d_wide$phase)

all_phases <- c("F", "O", "EL", "LL")
all_prots <- colnames(d_wide)[! colnames(d_wide) %in% c("SampleID", "ID", "TP","phase")]


# READ COVARIATES
covariates <- read.delim("results12/covariates_olink_batch12.phase_avg.txt", sep = "\t", check.names = F, as.is = T)

if (! "ID" %in% colnames(covariates)){
  covariates$ID <- gsub("_.*", "", covariates$SampleID)
  covariates$phase <- gsub(".*_", "", covariates$SampleID)
  covariates <- covariates %>%
    dplyr::select(SampleID, ID, phase, everything())
}
covariates$phase <- relevel(factor(covariates$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
#covariates$TP <- as.numeric(covariates$phase)

covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})


covariate_names <- c("Age","BMI","batch", "from", "storage_months")
covariates$from <- relevel(as.factor(covariates$from), ref = "X")

covariates_under30 <- covariates[covariates$Age <= 30, ]
covariates_above30 <- covariates[covariates$Age > 30, ]

# READ PHENOTYPES
pheno <- read.delim("../../phenotypes/batch12/cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.phase_avg.txt", as.is = T, check.names = F, sep = "\t")            

all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
all_phenos <- c("INS", "HOMA_B", "HOMA_IR", "GL", "AST", "ALT", "TRI", "COL", "HDL", "LDL")
all_phenos_combined <- c(all_hormones, all_phenos)
covariate_names_pheno <- c("from", "Age", "BMI")

if (! "ID" %in% colnames(pheno)){
  pheno$ID <- gsub("_.*", "", pheno$SampleID)
  pheno$phase <- gsub(".*_", "", pheno$SampleID)
  pheno <- pheno %>%
    dplyr::select(SampleID, ID, phase, everything())
}

pheno$phase <- relevel(factor(pheno$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
#pheno$TP <- as.numeric(pheno$phase)


# RUN PROTEIN - HORMONE/PHENO ASSOCIATIONS IN 2 AGE SUBSETS

gam_res <- data.frame(matrix(nrow = length(all_prots) * length(all_hormones), ncol = 12))
colnames(gam_res) <- c("prot", "pheno", "pval_under30", "estimate_under30", "SE_under30", "n_under30", "n_samples_under30",
                       "pval_above30", "estimate_above30", "SE_above30", "n_above30", "n_ids_above30")

lmm_res <- data.frame(matrix(nrow = length(all_prots) * length(all_hormones), ncol = 16))
colnames(lmm_res) <- c("prot", "pheno", "estimate_under30", "pval_under30", "SE_under30", "tval_under30", "n_under30", "n_samples_under30", "singular_fit_under30",
                       "estimate_above30", "pval_above30", "SE_above30", "tval_above30", "n_above30", "n_ids_above30", "singular_fit_above30")

cnt <- 1
for (ph in all_hormones) {
  cat(ph, "\n")
  i = 1
  for (prot in all_prots){
    #res_gam_u30 <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates_under30, scale = T, adjust_timepoint = 'spline', anova_pval = F)
    #res_gam_a30 <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates_above30, scale = T, adjust_timepoint = 'spline', anova_pval = F)
    #gam_res[cnt,] <- c(prot, ph, unlist(res_gam_u30), unlist(res_gam_a30))
    
    res_lmm_u30 <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates_under30, scale = T, adjust_timepoint = 'cubic', report_singular = T)
    res_lmm_a30 <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates_above30, scale = T, adjust_timepoint = 'cubic', report_singular = T)
    
    lmm_res[cnt,] <- c(prot, ph, unlist(res_lmm_u30), unlist(res_lmm_a30))
    
    cnt <- cnt + 1
    i <- i + 1
  }
}

gam_res <- gam_res %>%
  na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval_under30) %>%
  mutate(BH_pval_under30 = p.adjust(pval_under30, method = 'BH')) %>%
  mutate(BH_pval_above30 = p.adjust(pval_above30, method = 'BH')) %>%
  relocate(BH_pval_under30, .after = pval_under30) %>%
  relocate(BH_pval_above30, .after = pval_above30)

lmm_res <- lmm_res %>%
  na.omit(lmm_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval_under30) %>%
  mutate(BH_pval_under30 = p.adjust(pval_under30, method = 'BH')) %>%
  mutate(BH_pval_above30 = p.adjust(pval_above30, method = 'BH')) %>%
  relocate(BH_pval_under30, .after = pval_under30) %>%
  relocate(BH_pval_above30, .after = pval_above30)

write.table(gam_res, file = paste0(out_basedir, "prot_vs_hormones_stratified_by_age.txt"), quote = F, sep = "\t", row.names = FALSE)

write.table(lmm_res, file = paste0(out_basedir, "prot_vs_hormones_stratified_by_age.LMM.txt"), quote = F, sep = "\t", row.names = FALSE)

