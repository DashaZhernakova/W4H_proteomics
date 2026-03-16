args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
library(dplyr)

source("/mnt/sannaLAB-Temp/dasha/olink/batch2/W4H_proteomics/utility_functions_generic_association_analysis.R")

fname1 <- args[1]
fname2 <- args[2]
covar_fname <- args[3]
prot_list_fname <- args[4]
out_fname <- args[5]

run_lmm <- T
run_gam <- F

d1 <- read_file_add_phase(fname1)
d2 <- read_file_add_phase(fname2)
covariates  <- read_file_add_phase(covar_fname)

feature1_list <- as.character(read.delim(prot_list_fname, as.is = T, check.names = F, sep = "\t", header = F))

all_features2 <- colnames(d2)[! colnames(d2) %in% c("ID", "SampleID", "TP", "phase")]

cat("Running association analysis for ", length(feature1_list), " features\n")
cat("File 1: ", fname1, "\n")
cat("File 2: ", fname2, "\n")
cat ("Covariates: ", colnames(covariates)[2:ncol(covariates)], "\n")

if (run_gam){
  gam_res <- data.frame(matrix(nrow = length(feature1_list) * length(all_features2), ncol = 7))
  colnames(gam_res) <- c("feature1", "feature2", "pval", "edf_round", "fval", "n", " n_samples")

  cnt <- 1
  for (prot in feature1_list) {
    cat(prot, "\n")
    for (ph in all_features2){
      res_gam <- gam_prot_pheno_adj_covar(d1, d2, prot, ph, covariates, scale = T, adjust_timepoint = 'spline')
      gam_res[cnt,] <- c(prot, ph, unlist(res_gam))
      cnt <- cnt + 1
    }
  }
  
  write.table(gam_res, file = paste0(out_fname, ".GAM_results.txt"), quote = F, sep = "\t", row.names = FALSE)

}
if (run_lmm) {
  lmm_res <- data.frame(matrix(nrow = length(prot_list) * (ncol(d2) -4), ncol = 8))
  colnames(lmm_res) <- c("feature1", "feature2", "estimate", "pval", "se", "tval", "N", "N_unique")

  cnt <- 1
  for (prot in feature1_list) {
    cat(prot, "\n")
    for (ph in all_features2){
      res_lmm <- lmm_pheno_prot_adj_covar(d1, d2, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
      lmm_res[cnt,] <- c(prot, ph, unlist(res_lmm))
      cnt <- cnt + 1
    }
  }
  write.table(lmm_res, file = paste0(out_fname, ".LMM_results.txt"), quote = F, sep = "\t", row.names = FALSE)
}

cat("Finished!\n")
