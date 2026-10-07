my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("//mnt/sannaLAB-Temp/dasha/olink/batch2/")
source("scripts/utility_functions.R")

library(ggplot2)
library(rmcorr)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)
library(pheatmap)
library(corrplot)
library(patchwork)

  library(future) 
  library(future.apply)

set.seed(123)

out_basedir <- "results/intensity_shared_subsampled_110526/"

d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.phase_avg.txt", as.is = T, check.names = F, sep = "\t")

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$phase <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    dplyr::select(SampleID, ID, phase, everything())
}

d_wide$phase <- relevel(factor(d_wide$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

covariates <- read.delim("data/covariates_olink_batch12.phase_avg.txt", sep = "\t", check.names = F, as.is = T)

if (! "ID" %in% colnames(covariates)){
  covariates$ID <- gsub("_.*", "", covariates$SampleID)
  covariates$phase <- gsub(".*_", "", covariates$SampleID)
  covariates <- covariates %>%
    dplyr::select(SampleID, ID, phase, everything())
}
covariates$phase <- relevel(factor(covariates$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})


# Make a dataframe with proteins adjusted for all covariates per visit
covariate_names <- c("Age","BMI","batch", "from", "storage_months")
covariates$from <- relevel(as.factor(covariates$from), ref = "X")

pheno <- read.delim("data/cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.phase_avg.txt", as.is = T, check.names = F, sep = "\t")            

all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
all_phenos <- c("INS", "HOMA_B", "HOMA_IR", "GL", "AST", "ALT", "TRI", "COL", "HDL", "LDL")
covariate_names_pheno <- c("from", "Age", "BMI")

if (! "ID" %in% colnames(pheno)){
  pheno$ID <- gsub("_.*", "", pheno$SampleID)
  pheno$phase <- gsub(".*_", "", pheno$SampleID)
  pheno <- pheno %>%
    dplyr::select(SampleID, ID, phase, everything())
}

pheno$phase <- relevel(factor(pheno$phase, levels = c("F", "O", "EL", "LL")), ref = "F")


shared_prots <- colnames(d_wide)[colSums(is.na(d_wide)) < 50]
d_wide_shared <- d_wide[ ,shared_prots]
all_prots <- colnames(d_wide_shared)[! colnames(d_wide_shared) %in% c("ID", "TP","phase","SampleID")]

d_wide_b2 <- d_wide[d_wide$SampleID %in% covariates[covariates$batch == 'batch2', "SampleID"],]

# remove overlapping prots from the b2 dataset
d_wide_b2 <- d_wide_b2[, !colnames(d_wide_b2) %in% all_prots]

dim(d_wide)
dim(d_wide_b2)
dim(d_wide_shared)

d_wide_full <- d_wide
d_wide <- d_wide_shared

all_phases <- c("F", "O", "EL", "LL")

################################################################################
# Functions
################################################################################

subsample_d_wide <- function(d_wide, d_wide_b2){
  nsamples_b2 <- nrow(d_wide_b2)
  nindiv_b2 <- length(unique(d_wide_b2$ID))
  indiv_subset <- sample(unique(d_wide$ID), size = nindiv_b2, replace = F)
  d_wide_subset_tmp <- d_wide[d_wide$ID %in% indiv_subset,]
  sample_subset <- sample(d_wide_subset_tmp$SampleID, size = nsamples_b2, replace = F)
  
  return(d_wide[d_wide$SampleID %in% sample_subset,])                    
}

run_prot_pheno_gam <- function(d_wide, pheno, covariates, round_num){
  all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
  all_prots <- colnames(d_wide)[! colnames(d_wide) %in% c("SampleID", "ID", "phase", "TP")]
  
  gam_res <- data.frame(matrix(nrow = length(all_prots) * length(all_hormones), ncol = 8))
  colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")
  
  cnt <- 1
  for (ph in all_hormones) {
    cat(ph, "\n")
    
    pb <- txtProgressBar(min = 1, max = length(all_prots), style = 3)
    i = 1
    for (prot in all_prots){
      res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F, longitudinal = T)
      
      #res_lmm <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic', longitudinal = T)
      
      gam_res[cnt,] <- c(prot, ph, unlist(res_gam), res_lmm$pval)
      cnt <- cnt + 1
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  gam_res <- gam_res %>%
    na.omit(gam_res) %>%
    mutate(across(-c(pheno, prot), as.numeric)) %>% 
    arrange(pval) %>%
    mutate(BH_pval = p.adjust(pval, method = 'BH'))
  
    
  cat("Round", round_num, ":", nrow(gam_res[gam_res$BH_pval < 0.05,]), "\n")
  
  write.table(gam_res, file = paste0(out_basedir, "prot_vs_hormones.round", round_num, ".txt"), quote = F, sep = "\t", row.names = FALSE)

}

run_prot_pheno_gam_parallel <- function(d_wide, pheno, covariates, round_num){
  all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
  all_prots <- colnames(d_wide)[! colnames(d_wide) %in% c("SampleID", "ID", "phase", "TP")]
  
  # 1. Create a grid of all combinations to process
  tasks <- expand.grid(prot = all_prots, ph = all_hormones, stringsAsFactors = FALSE)
  
  # 2. Setup parallel backend (use most available cores)
  plan(multisession, workers = parallel::detectCores() - 1)
  
  # 3. Use future_lapply to run calculations in parallel
  results_list <- future_lapply(1:nrow(tasks), function(idx) {
    prot <- tasks$prot[idx]
    ph <- tasks$ph[idx]
    
    # Run your GAM function
    res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, 
                                        scale = T, adjust_timepoint = 'spline', 
                                        anova_pval = F, longitudinal = T)
    
    # Return as a named vector or list
    return(c(prot = prot, pheno = ph, unlist(res_gam))) 
  }, future.seed = TRUE) # Important for reproducibility/random processes in GAMs
  
  # 4. Combine results
  gam_res <- as.data.frame(do.call(rbind, results_list))
  colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples")
  
  # 5. Clean up
  gam_res <- gam_res %>%
    mutate(across(c(pval, estimate, SE, n, n_samples), as.numeric)) %>% 
    na.omit() %>%
    arrange(pval) %>%
    mutate(BH_pval = p.adjust(pval, method = 'BH'))
  
  # Save...
  write.table(gam_res, file = paste0(out_basedir, "prot_vs_hormones.round", round_num, ".txt"), 
              quote = F, sep = "\t", row.names = FALSE)
}

################################################################################
# Run associations
################################################################################

round_num = 1
d_wide_subs <- subsample_d_wide(d_wide, d_wide_b2)
dim(d_wide_subs)

run_prot_pheno_gam_parallel(d_wide_subs, pheno, covariates, round_num)
