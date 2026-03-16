library(dplyr)
library(lme4)
library(lmerTest)
library(mgcv)

read_file_add_phase <- function(path, add_phase = T){
  d <- read.delim(path, sep = "\t", check.names = F, as.is = T)
  
  if (! "ID" %in% colnames(d)){
    d$ID <- gsub("_.*", "", d$SampleID)
    if (add_phase){
      d$phase <- gsub(".*_", "", d$SampleID)
      d <- d %>%
        dplyr::select(SampleID, ID, phase, everything())
    } else {
      d$TP <- gsub(".*_", "", d$SampleID)
      d <- d %>%
        dplyr::select(SampleID, ID, TP, everything())
    }
  }
  if (add_phase) d$phase <- relevel(factor(d$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
  
  d[] <- lapply(d, function(col) {
    if (length(unique(col)) < 3) {
      return(factor(col))
    } else {
      return(col)
    }
  })
  
  return(d)
}


gam_prot_pheno_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, rm_outliers = F, adjust_timepoint = 'spline', adjust_pheno = 'linear', anova_pval = F, predict = F, add_age_interaction = F, longitudinal = T){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    pheno$TP <- as.numeric(pheno$phase)
    pheno$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase <- NULL
  }
  
  if(! ("TP" %in% colnames(covariates) || "phase" %in% colnames(covariates)) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID", "TP"))
    d_subs$SampleID.y <- NULL
    covariates$TP = NULL
  }
  
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$ID <- as.factor(d_subs$ID)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
  # generate the GAM formula
  if (longitudinal){
    if (adjust_timepoint == 'spline'){
      fo_gam <- paste("prot ~ s(pheno) + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
      fo_gam_null <- paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
    } else if (adjust_timepoint == 'linear') {
      fo_gam <- paste("prot ~ s(pheno) + TP + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
      fo_gam_null <- paste("prot ~ TP + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
    } else if (adjust_timepoint == 'none') {
      fo_gam <- paste("prot ~ s(pheno) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
      fo_gam_null <- paste("prot ~  s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
    } else {
      stop ("Wrong adjust_timepoint argument. Should be one of spline, linear or none.")
    }
  } else {
    fo_gam <- paste("prot ~ s(pheno) + ", paste(covariate_names, collapse = "+"))
    fo_gam_null <- paste("prot ~ ", paste(covariate_names, collapse = "+"))
  }
  
  
  if (add_age_interaction) {
    if (! "Age" %in% colnames(d_subs)) {cat ("No Age covariate provided for the interaction!\n")}
    fo_gam <- paste0(fo_gam, " + pheno * Age")
    d_subs$Age <- scale(d_subs$Age)
  }
  
  # Linear relation between protein and phenotype
  if (adjust_pheno != 'spline'){
    fo_gam <- gsub("s\\(pheno\\)", "pheno", fo_gam)
    
    model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
    
    est <- summary(model)$p.table["pheno","Estimate"]
    se <- summary(model)$p.table["pheno","Std. Error"]
    pval <- summary(model)$p.table["pheno","Pr(>|t|)"]
    
    if (add_age_interaction) {
      interaction_pval <- summary(model)$p.table["pheno:Age","Pr(>|t|)"]
      return(list(pval = pval,  est = est, se = se, n = nrow(d_subs), n_samples = length(unique(d_subs$ID)), age_inter_pval = interaction_pval))
    } 
    return(list(pval = pval,  est = est, se = se, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
  }
  
  # NON-linear relation between protein and phenotype
  model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
  
  edf <- round(summary(model)$s.table["s(pheno)","edf"])
  fval <- summary(model)$s.table["s(pheno)","F"]
  
  if (anova_pval){
    model0 <- gam(as.formula(fo_gam_null), data = d_subs, method = 'REML')
    an <- anova.gam(model, model0)
    pval <- an$`Pr(>F)`[2]
  } else {
    pval <- summary(model)$s.table["s(pheno)","p-value"]
  }
  
  return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}


lmm_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic", longitudinal = T){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    pheno$TP <- as.numeric(pheno$phase)
    pheno$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase <- NULL
  }
  
  if(! ("TP" %in% colnames(covariates) || "phase" %in% colnames(covariates)) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID", "TP"))
    d_subs$SampleID.y <- NULL
    covariates$TP = NULL
  }
  
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
  if (longitudinal){
    if (adjust_timepoint == 'cubic'){
      fo_lmm <- as.formula(paste("prot ~ poly(TP, 3) + pheno +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
      fo_lmm_base <- as.formula(paste("prot ~ poly(TP, 3) + ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
    } else if (adjust_timepoint == 'linear') {
      fo_lmm <- as.formula(paste("prot ~ TP + pheno +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
      fo_lmm_base <- as.formula(paste("prot ~ TP + ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
    } else if (adjust_timepoint == 'none') {
      fo_lmm <- as.formula(paste("prot ~ pheno +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
      fo_lmm_base <- as.formula(paste("prot ~ ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
    } else {
      stop ("Wrong adjust_timepoint argument. Should be one of cubic, linear or none.")
    }
    model <- lmer(fo_lmm, data = d_subs)
    
  } else {
    fo_lmm <- as.formula(paste("prot ~ pheno +", paste(covariate_names, collapse = "+")))
    fo_lmm_base <- as.formula(paste("prot ~ ", paste(covariate_names, collapse = "+")))
    model <- lm(fo_lmm, data = d_subs)
  }
  est <- summary(model)$coefficients["pheno", "Estimate"]
  se <- summary(model)$coefficients["pheno","Std. Error"]
  tval <- summary(model)$coefficients["pheno","t value"]
  pval <- summary(model)$coefficients["pheno","Pr(>|t|)"]
  
  return(list(estimate = est, pval = pval, se = se, tval = tval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

