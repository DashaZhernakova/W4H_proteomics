
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
#my_colors <- c("#eda048", "#4d9bb0", "#075c62", "#86e6ca")
#my_colors <- c("#f5b026", "#e26413", "#9a4020", "#63552d")
#my_colors <- c("#097054", "#FFDE00", "#6599FF", "#FF9900")

#library(tidyverse)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(tidyverse)
library(mgcv)
library(ggtext)
library(RColorBrewer)
library(limma)


gam_prot_tp_adj_covar <- function(d_wide, prot, covariates, scale = F, rm_outliers = F, predict = T, anova_pval = F, n_points = 20){
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase = NULL
  }
  
  d_subs <- inner_join(d_wide[,c(prot, "SampleID", "ID", "TP")], covariates, by = c("SampleID", "ID", "TP"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$ID <- as.factor(d_subs$ID)
  d_subs <- na.omit(d_subs)
  
  if (rm_outliers) d_subs <- remove_outliers_zscore(d_subs, "prot")
  
  if (scale) d_subs$prot <- scale(d_subs$prot)
  
  covariate_names = colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
  fo_gam <- as.formula(paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+")))
  fo_gam_null <- as.formula(paste("prot ~ s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+")))
  
  model <- gam(fo_gam, data = d_subs,  method = 'REML')
  if (anova_pval){
    model0 <- gam(fo_gam_null, data = d_subs, method = 'REML')
    an <- anova.gam(model, model0)
    pval <- an$`Pr(>F)`[2]
  } else {
    pval <- summary(model)$s.table["s(TP)","p-value"]
  }
  edf <- summary(model)$s.table["s(TP)","edf"]
  fval <- summary(model)$s.table["s(TP)","F"]
  
  if (predict){
    covar_means <- as.data.frame(lapply(covariates[,covariate_names], function(x) {
      if(is.numeric(x)) {
        mean(x, na.rm = TRUE)
      } else {
        levels(x)[1]  # Use first factor level
      }
    }))
    
    new_data <- expand.grid(
      TP = seq(1, 4, length.out = n_points),
      ID = unique(d_subs$ID),
      predicted = NA
    ) %>%
      bind_cols(
        covar_means[1,] 
      )
    
    predictions <- predict.gam(model, newdata = new_data,  exclude = "s(ID)", se.fit = T)
    new_data$predicted <- predictions$fit
    new_data$SE <- predictions$se.fit
    new_data$lower <- new_data$predicted - 1.96 * new_data$SE
    new_data$upper <- new_data$predicted + 1.96 * new_data$SE
    
    new_data2 <- unique(new_data[,c("TP", "predicted", "lower", "upper")])
    
    return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID)), predicted = new_data2$predicted, lower = new_data2$lower, upper = new_data2$upper))
  } 
  return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

gam_prot_pheno_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, rm_outliers = F, adjust_timepoint = 'spline', adjust_pheno = 'linear', anova_pval = F, predict = F, add_age_interaction = F, longitudinal = T){
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
  d_subs$ID <- as.factor(d_subs$ID)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
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
  #print(fo_gam)
  #print(str(d_subs))
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
  
  if (predict){
    covar_means <- as.data.frame(lapply(covariates[,covariate_names], function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else as.factor(2)))
    
    new_data <- expand.grid(
      TP = seq(1, 4),
      ID = unique(d_subs$ID),
      pheno = seq(min(d_subs$pheno), max(d_subs$pheno), length.out = 50),
      predicted = NA
    ) %>%
      bind_cols(
        covar_means %>%
          slice(1)   # to use the first row of covar_means
      )
    predictions <- predict.gam(model, newdata = new_data,  exclude = "s(ID)", se.fit = T)
    new_data$predicted <- predictions$fit
    new_data$SE <- predictions$se.fit
    new_data$lower <- new_data$predicted - 1.96 * new_data$SE
    new_data$upper <- new_data$predicted + 1.96 * new_data$SE
    
    new_data2 <- unique(new_data[,c("TP", "pheno","predicted", "lower", "upper")])
    new_data2$TP <- as.factor(new_data2$TP)
    
    return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID)), new_data = new_data2))
  } 
  
  return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

# test for association all hormones together
gam_prot_all_pheno_together_adj_covar <- function(d_wide, pheno, prot, covariates, scale = F, adjust_timepoint = 'spline'){
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    pheno$TP <- as.numeric(pheno$phase)
    pheno$phase <- NULL
  }
  
  colnames(pheno) <- gsub("17BES", "X17BES",colnames(pheno))
  if(! "TP" %in% colnames(covariates) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno, by = c("SampleID", "TP", "ID")),
                         covariates, by = c("SampleID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno, by = c("SampleID", "TP", "ID")),
                         covariates, by = c("SampleID","ID", "TP"))
    covariates$TP = NULL
  }
  colnames(d_subs)[1:4] <- c("SampleID", "ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$ID <- as.factor(d_subs$ID)
  
  #if (scale) {
  #  d_subs$prot <- scale(d_subs$prot)
  #  d_subs$pheno <- scale(d_subs$pheno)
  #}

  pheno_names <- colnames(pheno)[!colnames(pheno ) %in% c("ID", "SampleID", "TP")]
  if (adjust_timepoint == 'spline'){
    fo_gam <- paste("prot ~  s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
    fo_gam_null <- paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
  } else if (adjust_timepoint == 'linear') {
    fo_gam <- paste("prot ~  TP + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
    fo_gam_null <- paste("prot ~ TP + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
  } else if (adjust_timepoint == 'none') {
    fo_gam <- paste("prot ~  s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
    fo_gam_null <- paste("prot ~  s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of spline, linear or none.")
  }
  
  # Linear relation between protein and phenotype
  
  model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
  
  ests <- summary(model)$p.table[pheno_names,"Estimate"]
  ses <- summary(model)$p.table[pheno_names,"Std. Error"]
  pvals <- summary(model)$p.table[pheno_names,"Pr(>|t|)"]
  
  return(list(pvals = pvals,  ests = ests, ses = ses, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}


lmm_prot_tp_poly3_adj_covar <- function(d_wide, prot, covariates, scale = F){
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase = NULL
  }
  d_subs <- inner_join(d_wide[,c(prot, "SampleID","ID", "TP")], covariates, by = c("SampleID", "ID", "TP"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  covariate_names = colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
  fo_lmm <- as.formula(paste("prot ~ poly(TP,3) +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  model <- lmer(fo_lmm, data = d_subs)
  fo_lmm_base <- as.formula(paste("prot ~ ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  model0 <- lmer(fo_lmm_base, data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return(pval)
}

lmm_prot_tp_factor_adj_covar <- function(d_wide, prot, covariates, scale = F){
  d_subs <- inner_join(d_wide[,c(prot, "ID", "TP")], covariates, by = c("ID"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  fo_lmm <- as.formula(paste("prot ~ TP +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  model <- lmer(fo_lmm, data = d_subs)
  fo_lmm_base <- as.formula(paste("prot ~ ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  model0 <- lmer(fo_lmm_base, data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return(pval)
}

lmm_pheno_prot_no_adj_covar <- function(d_wide, pheno, prot, ph,  scale = F, adjust_timepoint = "cubic"){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  if (adjust_timepoint == 'cubic'){
    fo_lmm <- as.formula("prot ~ poly(TP, 3) + pheno + (1|ID)")
    fo_lmm_base <- as.formula("prot ~ poly(TP, 3) + (1|ID)")
  } else if (adjust_timepoint == 'linear') {
    fo_lmm <- as.formula("prot ~ TP + pheno + (1|ID)")
    fo_lmm_base <- as.formula("prot ~ TP +  (1|ID)")
  } else if (adjust_timepoint == 'none') {
    fo_lmm <- as.formula("prot ~ pheno + (1|ID)")
    fo_lmm_base <- as.formula("prot ~  (1|ID)")
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of cubic, linear or none.")
  }
  
  model <- lmer(fo_lmm, data = d_subs)
  model0 <- lmer(fo_lmm_base, data = d_subs)
  
  est <- summary(model)$coefficients["pheno", "Estimate"]
  se <- summary(model)$coefficients["pheno",2]
  tval <- summary(model)$coefficients["pheno",3]
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return(list(estimate = est, pval = pval, se = se, tval = tval))
}

lmm_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic", longitudinal = T){
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

lmm_prot_tp_interaction_pheno_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic"){
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("ID"))
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  if (adjust_timepoint == 'cubic'){
    fo_lmm <- as.formula(paste("prot ~ poly(TP, 3) + pheno +  poly(TP, 3) * pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ poly(TP, 3) + pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else if (adjust_timepoint == 'linear') {
    fo_lmm <- as.formula(paste("prot ~ TP + pheno + TP * pheno + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ TP + pheno + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of cubic or linear.")
  }
  
  model <- lmer(fo_lmm, data = d_subs)
  model0 <- lmer(fo_lmm_base, data = d_subs)
  
  #est <- summary(model)$coefficients["pheno", "Estimate"]
  #se <- summary(model)$coefficients["pheno",2]
  #tval <- summary(model)$coefficients["pheno",3]
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return( pval)
}

lm_per_tp_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F){
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
  
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("SampleID"))
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }

  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP","phase")]
  res_table <- data.frame()
  for (tp in unique(d_subs$TP)){

    fo_lm <- as.formula(paste("prot ~ pheno +", paste(covariate_names, collapse = "+")))
    model <- lm(fo_lm, data = d_subs[d_subs$TP == tp,])
    coefs <- summary(model)$coefficients
    res_table <- rbind(res_table, c(ph, prot, tp, coefs['pheno', 1], coefs['pheno', 4]))
  }
  colnames(res_table) <- c("pheno", "prot", "TP","estimate", "pval")
  if (phases){
    res_table <- res_table %>%
      mutate(
        phase = case_when(
            TP == "1" ~ "F",
            TP == "2" ~ "O", 
            TP == "3" ~ "EL",
            TP == "4" ~ "LL",
            .default = "other"
        ),
        .before = TP
      ) 
    res_table$TP <- NULL
  }
  return(res_table)
}

lmm_pheno_prot_inter <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  model <- lmer(pheno ~ prot + TP + prot*TP + (1|ID), data = d_subs)
  #est <- summary(model)$coefficients["prot", "Estimate"]
  model0 <- lmer(pheno ~ prot + TP + (1|ID), data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  return(pval)
}

get_ICC <- function(d_wide, prot){
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
  }
  
   d_subs <- d_wide[,c(prot, "ID", "TP")]
   colnames(d_subs)[1] <- "prot"
   
   d_subs$TP <- as.numeric(d_subs$TP)
   d_subs <- na.omit(d_subs)
   d_subs$prot <- scale(d_subs$prot)
   
   m <- lmer(prot ~ 1 + TP + (1|ID), data = d_subs)
   
   vc <- as.data.frame(VarCorr(m))
   var_ID <- vc$vcov[vc$grp == "ID"]  # Variance due to random effect (ID)
   var_residual <- vc$vcov[vc$grp == "Residual"]  # Residual variance
   total_var <- var_ID + var_residual  # Total variance (excluding fixed effects)
   prop_ID <- var_ID / total_var  # Proportion of variance explained by ID
   
   R2m <- performance::r2(m)$R2_marginal
   
   return (list(ICC = prop_ID, var_tp = R2m))
}

fit_lmm_poly3_adj_covar <- function(d_wide, prot, n = 10, covariates, scale = F, poly_raw = F){
  d_subs <- inner_join(d_wide[,c(prot, "SampleID", "ID", "TP")], covariates, by = c("ID"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$ID <- as.factor(d_subs$ID)
  if (scale) d_subs[,"prot"] <- scale(d_subs[,"prot"])
  
  # make all columns with less than 3 unique values as factors  
  d_subs[] <- lapply(d_subs, function(col) {
    if (length(unique(col)) < 3) {
      return(factor(col))
    } else {
      return(col)
    }
  })
  
  if (poly_raw){
    fo_lmm <- as.formula(paste("prot ~ poly(TP, 3, raw = TRUE) +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else {
    fo_lmm <- as.formula(paste("prot ~ poly(TP,3) +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  }
  lmm_fit <- lmer(fo_lmm, data = d_subs)
  
  covar_means <- as.data.frame(lapply(covariates[,-1], function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else as.factor(2)))
  new_data <- cbind(data.frame(TP = seq(1,4, length.out = n),  predicted = NA), covar_means)
  new_data$predicted <- predict(lmm_fit, newdata = new_data, re.form = NA)
  
  coef <- summary(lmm_fit)$coefficients[grepl("poly\\(TP", row.names(summary(lmm_fit)$coefficients)),1]
  return(list("predicted" = new_data$predicted, "coefficients" = coef))
}

run_limma<-function(joined_data, tp1, tp2) {
  df <-joined_data[joined_data$phase %in% c(tp1, tp2),]
  df$ID <- as.factor(df$ID)
  df$SampleID <- NULL
  df$phase <- factor(df$phase, levels = c(tp1,tp2))
  
  # design a model 
  formula <- reformulate(termlabels = c("0 + as.factor(phase)", covariate_names), 
                         response = NULL)
  design<-model.matrix(formula, data = df)
  colnames(design)[c(1,2)] <- c("phase1", "phase2")
  
  # specify the pairing
  corfit <- duplicateCorrelation(t(df[,all_prots]), design, block = df$ID)
  
  # make contrast - what to compare
  contrast<- makeContrasts(Diff = phase2 - phase1, levels=design)
  
  # apply linear model to each protein
  # Robust regression provides an alternative to least squares regression that works with less restrictive assumptions. Specifically, it provides much better regression coefficient estimates when outliers are present in the data
  fit<-lmFit(t(df[,all_prots]), design=design,  method="robust", correlation =
               corfit$consensus )
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=length(all_prots), adjust.method="BH", confint=TRUE)
  #DE_results$Bonferroni_signif <- ifelse(DE_results$P.Value < 0.05 / nrow(DE_results), T, F)
  return(DE_results)
}

run_wilcox <- function(joined_data_adj_covar, tp1, tp2) {
  joined_data_adj_covar$SampleID <- NULL
  wilcox_pvals <- data.frame(matrix(ncol = 3))
  colnames(wilcox_pvals) <- c("TP1_TP2", "prot", "wilcox_pval")
  cnt <- 1
  for (prot in all_prots){
    df <-joined_data_adj_covar[joined_data_adj_covar$phase %in% c(tp1, tp2), c("ID", "phase", prot)]
    df_wide <- na.omit(my_pivot_wider(df, row_names = "ID", names_from = "phase", values_from = prot))
    pval <- wilcox.test(df_wide[,1], df_wide[,2], paired = T)$p.value
    wilcox_pvals[cnt,] <- c(paste0(tp1, "_", tp2), prot, pval)
    cnt <- cnt + 1
  }
  wilcox_pvals$wilcox_pval <- as.numeric(wilcox_pvals$wilcox_pval)
  wilcox_pvals$BH_qval <- p.adjust(wilcox_pvals$wilcox_pval, method = 'BH')
  #wilcox_pvals$Bonferroni_signif <- ifelse(wilcox_pvals$wilcox_pval < 0.05 / nrow(wilcox_pvals), T, F)
  return(wilcox_pvals)
}

compare_variances_levene <- function(data) {
  data_long <- data %>%
    pivot_longer(-c(SampleID, ID, phase), 
                 names_to = "protein", values_to = "abundance")
  
  results <- data_long %>%
    group_by(protein) %>%
    summarise(
      levene_p = leveneTest(abundance ~ phase)$`Pr(>F)`[1],
      # Variance statistics
      var_by_phase = list({
        group_by(cur_data(), phase) %>%
          summarise(variance = var(abundance), .groups = 'drop')
      }),
      .groups = 'drop'
    ) %>%
    filter(!is.na(levene_p)) %>%
    mutate(
      adj_p = p.adjust(levene_p, "BH"),
      significant = adj_p < 0.05,
      max_var = map_dbl(var_by_phase, ~max(.x$variance)),
      min_var = map_dbl(var_by_phase, ~min(.x$variance)),
      fold_change = max_var / min_var
    ) %>%
    select(-var_by_phase) %>%
    arrange(levene_p)
  
  return(results)
}

get_mean_per_tp <- function(d_wide, prot){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  mean_prot_by_TP <- d_subs %>%
    group_by(TP) %>%
    summarize(mean_prot = mean(prot, na.rm = TRUE))
  return(mean_prot_by_TP)
}

correlate_per_id <- function(d_wide, pheno, prot, ph){
  d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  
  cor_res <- d_subs %>% 
    group_by(ID) %>%
    summarise(correl = cor(prot, pheno, method = 'spearman'))
  
  return(cor_res)
}

my_pivot_wider <- function(d, row_names, names_from, values_from){
  d2 <- d[,c(row_names, names_from, values_from)] %>%
    pivot_wider(names_from = {{names_from}}, values_from = {{values_from}})
  d2 <- as.data.frame(d2)
  row.names(d2) <- d2[,row_names]
  d2[,row_names ] <- NULL
  return(d2)
}

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


regress_covariates_lmm_phase <- function(data, covar_data, covars_longitudinal = T, keep_scale = F){
  phases = F
  if(! "TP" %in% colnames(data) & "phase" %in% colnames(data)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    data$TP <- as.numeric(data$phase)
    data$phase <- NULL
  }
  
  if (!"SampleID" %in% colnames(covar_data) & covars_longitudinal) {
    covar_data <- cbind(paste0(covar_data$ID, "_",covar_data$TP), covar_data)
    colnames(covar_data)[1] <- "SampleID"
  }
  
  d_adj <- data[,c("SampleID", "ID", "TP")]
  
  data[,"TP"] <- NULL
  covar_data[,"TP"] <- NULL
  
  covar_names = colnames(covar_data)[! colnames(covar_data) %in% c("SampleID", "ID", "TP", "phase")]
  cnt <- 1
  for (ph in colnames(data)[3: (ncol(data))]){
    if (covars_longitudinal){
      covar_data$ID <- NULL
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    } else {
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "ID"))
    }
    colnames(subs)[3] <- 'pheno'
    
    if (length(unique(subs$ID)) == length(subs$ID)) { # if no repeated measurements
      fo_lm <- as.formula(paste("pheno ~ ", paste(covar_names, collapse = "+")))
      lm_fit <- lm(fo_lm, data = subs)
      if (!keep_scale){
        subs[,ph] <- residuals(lm_fit)
      } else { # keep the original scale and global mean
        intercept <- coef(lm_fit)[1]
        subs[,ph] <- subs$pheno - (predict(lm_fit) - intercept)
      }

    } else {
      fo_lmm <- as.formula(paste("pheno ~ ", paste(covar_names, collapse = "+"), "+ (1|ID)"))
      lmm_fit <- lmer(fo_lmm, data = subs)
      if (!keep_scale){
        subs[,ph] <- subs$pheno - lme4:::predict.merMod(lmm_fit, re.form = NA)
      } else { # keep the original scale and global mean
        intercept <- fixef(lmm_fit)[1]
        predicted <- lme4:::predict.merMod(lmm_fit, re.form = NA)
        subs[,ph] <- subs$pheno - (predicted - intercept)
      }
    }
    d_adj <- left_join(d_adj, subs[, c("SampleID", ph)], by = "SampleID")
  }
  
  if(phases){
    cat ("renaming TP to phase\n")
    d_adj <- rename_TP_to_phase(d_adj)
  }
  
  return(d_adj)
}

rename_TP_to_phase <- function(d) {
  if (! "TP" %in% colnames(d)) {
    cat("error during converting visit to phase: no TP column!\n")
    return (d)
  }
  d %>%
    mutate(
      phase = factor(
        case_when(
          TP == 1 ~ "F",
          TP == 2 ~ "O", 
          TP == 3 ~ "EL",
          TP == 4 ~ "LL",
          .default = "other"
        ),
        levels = c("F", "O", "EL", "LL", "other")  # Specify factor levels
      ),
      .before = TP
    ) %>%
    dplyr::select(-TP)
}


plot_together <- function(d_wide = NULL, pheno = NULL, prot, ph, annot = "", method = "gam", scale = T, trajectories = NULL){
  if(! is.null(trajectories)){
    prot_name <- sym(prot)
    ph_name <- sym(ph)
    traj_t <- as.data.frame(t(trajectories)) %>%
      rownames_to_column(var = "TP")
    traj_t$TP <- as.numeric(traj_t$TP)
    g <- ggplot(traj_t) +
      geom_line(aes(x = TP, y = !!prot_name), color = my_colors[2]) +
      geom_line(aes(x = TP, y = !!ph_name), color = my_colors[3]) +
      labs(
        x = "Timepoint",
        y = "",
        title = paste0(
          "<span style='color:", my_colors[2], "'>", prot_name, "</span> vs ",
          "<span style='color:", my_colors[3], "'>", ph_name, "</span>"
        )
      ) +
      theme_minimal() + theme(plot.title = element_markdown())
    return(g)
  }
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    pheno$TP <- as.numeric(pheno$phase)
    d_wide$phase <- NULL
    pheno$phase <- NULL
  }
  
  if (is.character(d_wide$TP)) d_wide$TP <- as.numeric(d_wide$TP)
  d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], pheno[,c("ID", "TP", ph)], by = c("ID", "TP"))
  colnames(d_subs) <- c("ID", "TP", "prot", "pheno")
  
  d_subs$TP <- as.numeric(d_subs$TP)
  if (scale){
    d_subs$prot <- scale( d_subs$prot)
    d_subs$pheno <- scale( d_subs$pheno)
  }
  plot_title <- paste0(ph, " - ", prot)
  if (annot != "") plot_title <- paste0(plot_title, ", ", annot)
  if (method == "smooth"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  }  else if (method == "poly3"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  } else if(method == 'boxplot'){
    d_subs_long <- d_subs[,-1] %>% pivot_longer(names_to = 'type', cols = c('prot', 'pheno'), values_to = 'value')
    colnames(d_subs_long)[3] <- "value"
    d_subs_long$TP <- as.factor(d_subs_long$TP)
    g <- ggplot(d_subs_long, aes(x = TP,  fill = type, y = value)) +
      geom_boxplot(position = position_dodge(width = 0.85), width = 0.8) + 
      labs(x = "Timepoint ", y = "", title = plot_title) +
      theme_minimal() +
      stat_summary(
        fun = median,
        geom = 'line',
        aes(group = type, color = type),
        position = position_dodge(width = 0.85),
        linewidth = 1
      ) +
      scale_fill_manual(values = c(my_colors[c(2,3)])) + 
      scale_color_manual(values = c('darkgoldenrod4', 'deepskyblue4'))
  } else if (method == 'gam'){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(method = 'gam', formula=y ~ s(x, k = 4) , color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(method = 'gam', formula=y ~ s(x, k = 4), color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  }
  
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
  
}

plot_traj_many_prots2 <- function(prot_trajs = NULL, prots, colored = T, signif = NULL, title = "", phases = F){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  tmp <- prot_trajs[prots,]
  colnames(tmp) <- colnames(prot_trajs)
  d_subs <- tmp %>%
    rownames_to_column(var = 'prot') %>%
    pivot_longer(cols = -prot, names_to = 'TP')
  
  
  d_subs$TP <- as.numeric(d_subs$TP)
  
  if (colored){
    g <- ggplot(d_subs, aes(x = TP, color = prot, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F) +
      theme_minimal()
  } else {
    g <- ggplot(d_subs, aes(x = TP, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, alpha = 0.5) +
      theme_minimal() 
    if (!is.null(signif)){
      g <- g + geom_line(data = d_subs[d_subs$prot %in% signif,],aes(x = TP, y = value, group = prot), stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F,  color = 'red')
    }
    
    if (title != ""){
      g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
    }
  }
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
}

plot_traj_prots_and_pheno <- function(d_wide, pheno, prots, ph, title = "", method = 'gam', prot_trajs = NULL, ph_trajs = NULL, phases = F){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  if (!is.null(prot_trajs)) {
    res_trajs <- rbind(ph_trajs, prot_trajs) %>%
      rownames_to_column(var = 'feature')
    
    res_trajs[1, "feature"] <- "pheno"
  } else {
    res_trajs <- data.frame(matrix(nrow = length(prots) + 1, ncol = 101))
    ph_fit <- fit_lmm_poly3_adj_covar(pheno, ph, n = 100, covariates, scale = T)
    res_trajs[1,] <- c("pheno", ph_fit$predicted)
    
    cnt <- 2
    for (prot in prots){
      if (method == 'lmm') {
        fit <- fit_lmm_poly3_adj_covar(d_wide, prot, n = 100, covariates, scale = T)
      } else if (method == 'gam') {
        fit <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = T)
      } else {
        stop("Error! Wrong method, should be gam or lmm.")
      }
      res_trajs[cnt,] <- c(prot, fit$predicted)
      cnt <- cnt + 1
    }
    colnames(res_trajs) <- c("feature", seq(1,4, length.out = 100))
    
  }
  d_subs <- res_trajs %>%
    pivot_longer(cols = -feature,names_to = 'TP')
  
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$value <- as.numeric(d_subs$value)
  d_subs$feature_type <- ifelse(d_subs$feature == 'pheno', "phenotype" ,"proteins")  
  g <- ggplot(d_subs, aes(x = TP, color = feature_type, y = value, group = feature)) +
    geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, linewidth = 0.5, alpha = 0.4) +
    geom_line(data = d_subs[d_subs$feature_type == 'phenotype',], color = my_colors[2], stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, linewidth = 1) +
    theme_minimal() +
    scale_color_manual(values = my_colors[c(2,3)])
  
  if (title != ""){
    g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
  }
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
}

plot_medians_prots_and_pheno <- function(d_wide, pheno, prots, ph, title = "", scale = T,phases = F){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  d_subs <- inner_join(pheno[,c("ID", "TP", ph)], d_wide[,c("ID", "TP", prots)], by = c("ID", "TP"))
  colnames(d_subs)[1:3] <- c("ID", "TP", "pheno")
  
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs$ID <- NULL
  
  d_subs <- na.omit(d_subs)
  if (scale){
    d_subs[,-1] <- scale(d_subs[,-1])
  }
  
  medians <- aggregate(. ~ TP, data=d_subs, FUN=median) %>%
    pivot_longer(-TP, names_to = 'feature')
  
  
  medians$TP <- as.numeric(medians$TP)
  medians$value <- as.numeric(medians$value)
  medians$feature_type <- ifelse(medians$feature == 'pheno', "phenotype" ,"proteins")  
  g <- ggplot(medians, aes(x = TP, color = feature_type, y = value, group = feature)) + 
    geom_point() + 
    geom_line(linewidth = 0.5, alpha = 0.4) + 
    geom_line(data = medians[medians$feature_type == 'phenotype',], color = my_colors[2], linewidth = 1) +
    theme_minimal() +
    scale_color_manual(values = my_colors[c(2,3)])
  
  if (title != ""){
    g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
  }
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
}

make_radian_plot <- function(d_wide, prots, value = 'mean'){
  library(ggradar)
  if (value == 'mean'){
    d_long <- d_wide %>%
      pivot_longer(names_to = 'Assay', cols = -c('ID', 'TP', 'SampleID'), values_to = 'NPX')
    d_subs <- d_long[d_long$Assay %in% prots,]
    mean_prot_by_TP <- d_subs %>%
      group_by(TP, Assay) %>%
      summarize(value = 1 + mean(NPX, na.rm = TRUE))
  }
  
  data_wide <- mean_prot_by_TP %>%
    pivot_wider(names_from = Assay, values_from = value)
  
  data_wide$TP <- as.factor(data_wide$TP)
  
  g <- ggradar(
    data_wide,
    group.colours = my_colors,
    legend.title = "Timepoints",
    axis.label.size = 2,
    grid.label.size = 4,
    group.line.width = 1,
    group.point.size = 0,
    background.circle.colour = "white",
    gridline.mid.colour = "gray"
  )
  g
}

plot_clusters <- function(cl, method = "", num_k = "", colored = F, signif = NULL, save_pdf = T, add_cluster_name = F, out_path = NA, prot_trajs = NULL){
  plot_list = list()
  for (cluster in unique(cl)){
    title = ifelse(add_cluster_name, paste0(cluster, ", N = ", length(cl[cl == cluster])), "")
    plot_list[[cluster]] <- plot_traj_many_prots2(prot_trajs, names(cl[cl == cluster]), colored = colored, signif, title)
  }
  
  pdf_path = ifelse(is.na(out_path), 
                    paste0("../plots/clustering_signif_v2/", method, "_k", num_k, ".pdf"),
                    out_path)
  
  
  ncols = 4
  nrows = ceiling(length(unique(cl))/4)
  
  if(save_pdf) pdf(pdf_path, width = 4*ncols, height = 4*nrows)
  grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)
  
  
  #if (length(unique(cl)) < 10){
  #  grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
  #} else if(length(unique(cl)) < 10) {
  #  grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)  
  #} else {
  #  grid.arrange(grobs = plot_list, ncol = 5, nrow = 5)  
  #}
  if(save_pdf) dev.off()
  
}

plot_association_heatmap <- function(assoc_df, prot_subs, rows = 'pheno', cols = 'prot', vals = 'estimate', signif_vals = 'BH_pval', transpose = F, cutrows = NA, cutcols = NA, cluster_cols = T, col_order = NULL){
  assoc_df_wide <- my_pivot_wider(assoc_df[assoc_df$prot %in% prot_subs,], rows, cols, vals)
  signif_labels <- my_pivot_wider(assoc_df[assoc_df$prot %in% prot_subs,], rows, cols, signif_vals)
  signif_labels <- ifelse(signif_labels < 0.05, "*", "")
  
  if (transpose){
    assoc_df_wide <- as.data.frame(t(assoc_df_wide))
    signif_labels <- as.data.frame(t(signif_labels))
    fontsize_row = 8
    fontsize_col = 10
  }
  if (!is.null(col_order)) {
    assoc_df_wide <- assoc_df_wide[, col_order, drop = FALSE]
    signif_labels <- signif_labels[, col_order, drop = FALSE]
  }
  
  max_val <- max(abs(min(assoc_df_wide)), max(assoc_df_wide))
  breaksList = seq(-max_val, max_val, by = 0.01)
  colorList <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
  h <- pheatmap(assoc_df_wide, display_numbers = signif_labels, fontsize_number = 10, 
                fontsize_col = fontsize_col, fontsize_row = fontsize_row,
                color = colorList, breaks = breaksList, cutree_rows = cutrows, 
                cutree_cols = cutcols, cluster_cols = cluster_cols)
  
  h
}

plot_association_volcano <- function(assoc_df){
  if ("PROG" %in% assoc_df$pheno) desired_order <- c("PROG", "17BES", "LH", "FSH","PRL")
  if (!"PROG" %in% assoc_df$pheno) desired_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
  ggplot(assoc_df, aes(x = estimate, y = -log10(BH_pval))) +
    geom_point(aes(color = BH_pval < 0.05), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("grey", "red"), 
                       labels = c("FALSE" = "Not significant", "TRUE" = "FDR < 0.05")) +
    ggrepel::geom_text_repel(
      data = subset(assoc_df, BH_pval < 0.05),
      aes(label = prot), 
      size = 2,
      max.overlaps = 20
    ) +
    facet_wrap(~ factor(pheno, levels = desired_order)) +
    theme_minimal() +
    labs(color = "Significance")
}

# scatter colored by visit to see the relationship at each visit
scatter_col_tp <- function(d_wide, pheno, prot, ph, scale = F, add_points = F){
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    pheno$TP <- as.numeric(pheno$phase)
    d_wide$phase <- NULL
    pheno$phase <- NULL
  }
  
  if ("SampleID" %in% colnames(d_wide)){
    d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
    colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  } else {
    d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], pheno[,c("ID", "TP" ,ph)], by = c("ID", "TP"))
    colnames(d_subs) <- c("ID", "TP", "prot", "pheno")
  }
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if(scale){
    d_subs$pheno <- scale(d_subs$pheno)
    d_subs$prot <- scale(d_subs$prot)
  }
  g <- ggplot(d_subs, aes(x = prot, y = pheno, colour = TP)) + 
    geom_smooth(method = 'lm', alpha = 0.2) + 
    stat_smooth(method = 'lm', se = F) +
    theme_minimal() +
    labs(x = prot, y = ph, 
         title = paste0(ph, " - ", prot))  +
    scale_color_manual(values = my_colors)
  
  if (add_points) g <- g + geom_point()
  if(phases) {
    g <-g + scale_color_manual(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL"),
      values = my_colors
    )
  }
  g
}

plot_boxplot_with_traj <- function(d, d_adj, prot, covariates, add_pval = T) {
  res_gam <- gam_prot_tp_adj_covar(d, prot, covariates, scale = F, predict = T)
  
  if("phase" %in% colnames(d_adj)) {
    d_adj$TP <- as.numeric(d_adj$phase)
  }
  
  traj <- data.frame(TP = seq(1,4, length.out = length(res_gam$predicted)), 
                     pheno = res_gam$predicted, 
                     lower = res_gam$lower, 
                     upper = res_gam$upper)
  
  traj <- full_join(d_adj[,c("TP", "ID", prot)], traj, by = "TP")
  colnames(traj)[3] <- "values"
  
  # Scale the adjusted values using raw data statistics
  d_raw_subset <- inner_join(d[,c(prot, "SampleID")], covariates, by = "SampleID")
  d_raw_subset <- na.omit(d_raw_subset) 
  raw_mean <- mean(d_raw_subset[[prot]], na.rm = T)
  raw_sd   <- sd(d_raw_subset[[prot]], na.rm = T)
  
  traj$values <- (traj$values - raw_mean) / raw_sd
  traj$pheno <- (traj$pheno - raw_mean) / raw_sd
  traj$lower <- (traj$lower - raw_mean) / raw_sd
  traj$upper <- (traj$upper - raw_mean) / raw_sd
  
  
  p <- ggplot(traj) + 
    geom_boxplot(aes(x = TP, y = values, group = TP), width = 0.1, color = my_colors[3], outliers = F) +
    geom_line(aes(x = TP, y = pheno), color = my_colors[4]) + 
    geom_ribbon(aes(x = TP, ymin = lower, ymax = upper), alpha = 0.2, fill = my_colors[4]) +
    geom_jitter(aes(x = TP, y = values), alpha = 0.2, width = 0.1, color = my_colors[3]) +
    labs(x = "Phase ", y = prot) +
    theme_minimal() + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
 if (add_pval){
   p <- p + ggtitle(paste0("GAM P = ", formatC(res_gam$pval, digits = 2)))
 }
  p 
}




plot_boxplot_with_traj <- function(d_adj, prot) {
  d_adj_gam <- na.omit(d_adj[,c("SampleID", "ID", "phase", prot)])
  colnames(d_adj_gam)[4] <- "prot"
  if("phase" %in% colnames(d_adj_gam)) {
    d_adj_gam$TP <- as.numeric(d_adj_gam$phase)
  }
  d_adj_gam$ID <- as.factor(d_adj_gam$ID)

  model <- gam(prot ~ s(TP, k = 4) + s(ID, bs = 're'), 
               data = d_adj_gam, method = 'REML')
  
  new_data <- data.frame(
    TP = seq(1, 4, length.out = 20),
    ID = unique(d_adj_gam$ID)[1] 
  )
  
  predictions <- predict(model, newdata = new_data, 
                         exclude = "s(ID)", se.fit = TRUE)
  
  traj <- data.frame(
    TP = new_data$TP,
    pheno = predictions$fit,
    lower = predictions$fit - 1.96 * predictions$se.fit,
    upper = predictions$fit + 1.96 * predictions$se.fit
  )
  
  p <- ggplot(d_adj_gam) +
    geom_boxplot(aes(x = as.numeric(TP), y = prot, group = TP), 
                 width = 0.1, color = my_colors[3], outliers = FALSE) +
    geom_line(data = traj, aes(x = TP, y = pheno), color = my_colors[4]) +
    geom_ribbon(data = traj, aes(x = TP, ymin = lower, ymax = upper), 
                alpha = 0.2, fill = my_colors[4]) +
    geom_jitter(aes(x = as.numeric(TP), y = prot), 
                alpha = 0.2, width = 0.1, color = my_colors[3]) +
    labs(x = "Phase", y = paste0(prot, "\n(covariate-adjusted)")) +
    theme_minimal() +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  
  return(p)
}
