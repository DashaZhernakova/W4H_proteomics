d_wide_100s <- read.delim("data/olink_clean_CVD+INF_rm_below_lod_more_100_samples_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
d_wide_100s$TP <- as.numeric(d_wide_100s$TP)

d_wide_30i <- read.delim("data/olink_clean_CVD+INF_rm_below_lod_more_30_indiv_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
d_wide_30i$TP <- as.numeric(d_wide_30i$TP)

d_wide_0.5perc <- read.delim("data/olink_clean_CVD+INF_rm_below_lod_prot_with_more_half_samples_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
d_wide_0.5perc$TP <- as.numeric(d_wide_0.5perc$TP)

d_wide <- read.delim("data/olink_clean_CVD+INF_rm_below_lod_keep_NA.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
d_wide_0.5perc$TP <- as.numeric(d_wide_0.5perc$TP)


################################################################################
# GAM protein vs phenotype 
################################################################################
all_prot <- colnames(d_wide_30i)[4:ncol(d_wide_30i)]

gam_res <- data.frame(matrix(nrow = length(all_prot) * (ncol(pheno) -4), ncol = 9))
colnames(gam_res) <- c("prot", "pheno", "pval", "edf_round", "fval", "n", " n_samples", "pval_100s", "pval_0.5perc")


cnt <- 1
for (ph in colnames(pheno)[4:ncol(pheno)]) {
  cat(ph, "\n")
  for (prot in all_prot){
    res_gam <- gam_prot_pheno_adj_covar(d_wide_30i, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F)
    pr1 <- NA
    pr2 <- NA
    if (prot %in% colnames(d_wide_100s)[4:ncol(d_wide_100s)]) pr1 <- res_gam[['pval']]
    if (prot %in% colnames(d_wide_0.5perc)[4:ncol(d_wide_0.5perc)]) pr2 <- res_gam[['pval']]
    gam_res[cnt,] <- c(prot, ph, unlist(res_gam), pr1, pr2)
    cnt <- cnt + 1
  }
}

gam_res2 <- gam_res %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
gam_res2$BH_pval <- p.adjust(gam_res2$pval)
gam_res2$BH_pval_100s <- p.adjust(gam_res2$pval_100s)
gam_res2$BH_pval_0.5perc <- p.adjust(gam_res2$pval_0.5perc)

gam_res <- gam_res[order(gam_res$pval_spline),]

old_res <- read.delim("results/pheno_batch2_prot_rm_outliers_4sd/prot_vs_pheno_linear_gam_adj_covar.txt", as.is = T, check.names = F, sep = "\t")
old_res <- old_res[,c(1,2,3,4,6,7, ncol(old_res) - 2)]

merged <- full_join(gam_res2, old_res, by = c("prot", "pheno"))
merged <- merged[merged$pheno != 'Glucose',]
nrow(merged[merged$BH_pval_spline < 0.05 & (is.na(merged$BH_pval) | merged$BH_pval > 0.05),])
nrow(merged[merged$BH_pval_spline > 0.05 & (!is.na(merged$BH_pval) & merged$BH_pval < 0.05), ])


################################################################################
# Protein vs TP GAM and LMM
################################################################################

gam_res_prot_tp <- data.frame(matrix(nrow = length(all_prot), ncol = 8))
colnames(gam_res_prot_tp) <- c("prot", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "pval_100s", "pval_0.5perc")

cnt <- 1
for (prot in all_prot){
  res_gam <- gam_prot_tp_adj_covar(d_wide_30i, prot, covariates, scale = T, predict = F)
  pr1 <- NA
  pr2 <- NA
  if (prot %in% colnames(d_wide_100s)[4:ncol(d_wide_100s)]) pr1 <- res_gam[['pval']]
  if (prot %in% colnames(d_wide_0.5perc)[4:ncol(d_wide_0.5perc)]) pr2 <- res_gam[['pval']]
  
  gam_res_prot_tp[cnt,] <- c(prot, unlist(res_gam), pr1, pr2)
  cnt <- cnt + 1
}

gam_res_prot_tp <- gam_res_prot_tp %>%
  mutate(across(-c( prot), as.numeric)) 

gam_res_prot_tp$BH_pval <- p.adjust(gam_res_prot_tp$gam_pval, method = 'BH')

old_res_tp <- read.delim("results/pheno_batch2_prot_rm_outliers_4sd/prot_vs_tp_gam_adj_age_bmi_preg_storage.txt", as.is = T, check.names = F, sep = "\t")
merged <- full_join(gam_res_prot_tp, old_res_tp, by = "prot")

nrow(merged[merged$gam_BH_pval < 0.05,])
nrow(merged[!is.na(merged$BH_pval) & merged$BH_pval < 0.05,])
nrow(merged[merged$gam_BH_pval < 0.05 & (is.na(merged$BH_pval) | merged$BH_pval > 0.05),])
nrow(merged[merged$gam_BH_pval > 0.05 & (!is.na(merged$BH_pval) & merged$BH_pval < 0.05), ])

d_wide_old <- read.delim("data/olink_clean_CVD+INF_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))

ggplot(d_wide_old, aes(x = TP, y = IL32, color = ID, group = ID)) + geom_point() + geom_line() + theme(legend.position = None)

ggplot(d_wide_old, aes(x = TP, y = IL32, color = ID, group = ID)) + 
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) + 
  theme_minimal() + theme(legend.position = "none")

ggplot(d_wide_30i, aes(x = TP, y = IL32, color = ID, group = ID)) + 
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) + 
  theme_minimal() + theme(legend.position = "none")
