my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/utility_functions.R")

library(ggplot2)
library(rmcorr)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)
library(pheatmap)
library(corrplot)
library(patchwork)

set.seed(123)

out_basedir <- "results12/intensity_shared_prots_visit_061125/"

d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$TP <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    select(SampleID, ID, TP, everything())
}

d_wide$TP <- as.numeric(d_wide$TP)

covariates <- read.delim("results12/covariates_per_id_olink_batch12.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))

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

shared_prots <- colnames(d_wide)[colSums(is.na(d_wide)) < 50]
d_wide_shared <- d_wide[ ,shared_prots]
all_prots <- colnames(d_wide_shared)[! colnames(d_wide_shared) %in% c("ID", "TP","TP","SampleID")]

d_wide_b2 <- d_wide[d_wide$ID %in% covariates[covariates$batch == 'batch2', "ID"],]

# remove overlapping prots from the b2 dataset
d_wide_b2 <- d_wide_b2[, !colnames(d_wide_b2) %in% all_prots]

dim(d_wide)
dim(d_wide_b2)
dim(d_wide_shared)

d_wide_full <- d_wide
d_wide <- d_wide_shared

all_TPs <- 1:4

joined_data <- full_join(covariates, d_wide, by = c("ID"), relationship = 'one-to-many')
d_wide_adj_covar <- regress_covariates_lmm(d_wide, covariates, covars_longitudinal = F)
joined_data_adj_covar <- full_join(covariates, d_wide_adj_covar, by = c("ID"), relationship = 'one-to-many')

write.table(d_wide_adj_covar, file = paste0(out_basedir, "olink_clean_adj_covariates.txt"), quote = F, sep = "\t", row.names = FALSE)


################################################################################
# ICC for each protein
################################################################################

icc <- data.frame(matrix(nrow = (ncol(d_wide) -4), ncol = 3))
colnames(icc) <- c("prot", "ICC", "R2_TP")
cnt <- 1
for (prot in all_prots){
  res <- get_ICC(d_wide, prot)
  icc[cnt,] <- c(prot, unlist(res))
  cnt <- cnt + 1
}

icc <- na.omit(icc) %>%
  mutate(across(-c( prot), as.numeric)) 
icc <- icc[order(icc$ICC, decreasing = F),]

cat("ICC ranges from", min(icc$ICC), "to", max(icc$ICC), "with a median of", median(icc$ICC), "\n")
ggplot(icc, aes(x=ICC, y =R2_TP)) + geom_point() + theme_minimal() + xlab("ICC (variance explained by ID)") + ylab("Marginal R2 (var explained by TP)")

write.table(icc, file = paste0(out_basedir, "ICC_per_protein.txt"), quote = F, sep = "\t", row.names = FALSE)

pdf(paste0(out_basedir, "plots/ICC_per_protein.pdf"))
p1 <- ggplot(icc, aes(x = ICC)) + 
  geom_density() + 
  theme_minimal() 

p2 <- ggplot(icc, aes(y = ICC)) + 
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x=element_blank())
p1 + p2

ggplot(icc, aes(x=ICC, y =R2_TP)) + geom_point() + theme_minimal() + xlab("ICC (variance explained by ID)") + ylab("Marginal R2 (var explained by TP)")

p3 <- ggplot(icc, aes(y = ICC)) + 
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x=element_blank())

dev.off()


################################################################################
# Differentially expressed proteins between visits
################################################################################
library("limma")

run_limma<-function(joined_data, tp1, tp2) {
  df <-joined_data[joined_data$TP %in% c(tp1, tp2),]
  df$ID <- as.factor(df$ID)
  df$SampleID <- NULL
  df$TP <- factor(df$TP, levels = c(tp1,tp2))
  
  # design a model 
  formula <- reformulate(termlabels = c("0 + as.factor(TP)", covariate_names), 
                         response = NULL)
  design<-model.matrix(formula, data = df)
  colnames(design)[c(1,2)] <- c("TP1", "TP2")
  
  # specify the pairing
  corfit <- duplicateCorrelation(t(df[,all_prots]), design, block = df$ID)
  
  # make contrast - what to compare
  contrast<- makeContrasts(Diff = TP2 - TP1, levels=design)
  
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
  
  return(DE_results)
}

run_wilcox <- function(joined_data_adj_covar, tp1, tp2) {
  joined_data_adj_covar$SampleID <- NULL
  wilcox_pvals <- data.frame(matrix(ncol = 3))
  colnames(wilcox_pvals) <- c("TP1_TP2", "prot", "wilcox_pval")
  cnt <- 1
  for (prot in all_prots){
    df <-joined_data_adj_covar[joined_data_adj_covar$TP %in% c(tp1, tp2), c("ID", "TP", prot)]
    df_wide <- na.omit(my_pivot_wider(df, row_names = "ID", names_from = "TP", values_from = prot))
    pval <- wilcox.test(df_wide[,1], df_wide[,2], paired = T)$p.value
    wilcox_pvals[cnt,] <- c(paste0(tp1, "_", tp2), prot, pval)
    cnt <- cnt + 1
  }
  wilcox_pvals$wilcox_pval <- as.numeric(wilcox_pvals$wilcox_pval)
  wilcox_pvals$BH_qval <- p.adjust(wilcox_pvals$wilcox_pval, method = 'BH')
  return(wilcox_pvals)
}

limma_res_all <- data.frame()
wilcox_res_all <- data.frame()
TP_comb <- t(combn(all_TPs, 2))
for (i in 1:nrow(TP_comb)){
  tp1 = TP_comb[i,1]
  tp2 = TP_comb[i,2]
  tp1 <- as.character(tp1)
  tp2 <- as.character(tp2)
  
  # limma
  limma_res <- run_limma(joined_data, tp1, tp2) %>%
    rownames_to_column(var = 'prot')
  if (nrow(limma_res[limma_res$adj.P.Val < 0.05,]) > 0) {
    limma_res_all <- rbind(limma_res_all, cbind(paste0(tp1, "_", tp2), limma_res))
  }
  
  # wilcoxon
  wilcox_res <- run_wilcox(joined_data_adj_covar, tp1, tp2)
  wilcox_res_all <- rbind(wilcox_res_all, wilcox_res)
}


colnames(limma_res_all)[1] <- "TP1_TP2"
limma_res_all$sign <- ifelse(limma_res_all$adj.P.Val < 0.05, T, F)

#all(limma_res_all[limma_res_all$adj.P.Val < 0.05, "prot"] %in% signif$prot)
write.table(limma_res_all, file = paste0(out_basedir, "limma_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(wilcox_res_all, file = paste0(out_basedir, "wilcox_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)

# Plot significant DEPs

# stacked barplot
results <- limma_res_all %>%
  mutate(
    Regulation = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < 0 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

summary_data <- results %>%
  filter(Regulation != "Not significant") %>%  # Exclude non-significant proteins
  group_by(TP1_TP2, Regulation) %>%
  summarize(Count = n(), .groups = "drop") %>%
  mutate(TP1_TP2 = factor(TP1_TP2, 
                                levels = c("1_2", "1_3", "1_4", "2_3", "2_4", "3_4")))

pdf(paste0(out_basedir, "limma_DEPs_barplot.pdf"))
ggplot(summary_data, aes(x = TP1_TP2, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Visit Comparison",
    y = "Number of DEP",
    fill = "Effect direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors[c(3,2)])
dev.off()

limma_res_all$signif_direction <- ifelse(limma_res_all$sign, ifelse(limma_res_all$logFC < 0, -1, 1), 0)
limma_heatmap <- my_pivot_wider(limma_res_all[limma_res_all$signif_direction != 0,], row_names = "prot", names_from = "TP1_TP2", values_from = "signif_direction")
limma_heatmap[is.na(limma_heatmap)] <- 0
pheatmap(limma_heatmap, color = c(my_colors[3], "white", my_colors[2]), legend = F, cluster_cols = F, filename = paste0(out_basedir, "limma_DEPs_heatmap.pdf"))
dev.off()

################################################################################
# Protein vs TP GAM and LMM
################################################################################

gam_res_prot_tp <- data.frame(matrix(nrow = length(all_prots), ncol = 7))
colnames(gam_res_prot_tp) <- c("prot", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")

cnt <- 1
for (prot in all_prots){
  res_gam <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = F)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(d_wide, prot, covariates)
  gam_res_prot_tp[cnt,] <- c(prot, unlist(res_gam), res_lmm)
  cnt <- cnt + 1
}

gam_res_prot_tp <- na.omit(gam_res_prot_tp) %>%
  mutate(across(-c( prot), as.numeric)) 

gam_res_prot_tp$gam_BH_pval <- p.adjust(gam_res_prot_tp$gam_pval, method = 'BH')
gam_res_prot_tp$lmm_BH_pval <- p.adjust(gam_res_prot_tp$lmm_pval, method = 'BH')

gam_res_prot_tp <- gam_res_prot_tp[order(gam_res_prot_tp$gam_pval),]
gam_res_prot_tp$gam_bonf_sign <- ifelse(gam_res_prot_tp$gam_pval < 0.05/60,T,F)
gam_res_prot_tp$gam_edf_round <- round(gam_res_prot_tp$gam_edf)
gam_res_prot_tp$gam_BH_sign <- ifelse(gam_res_prot_tp$gam_BH_pval < 0.05,T,F)

write.table(gam_res_prot_tp, file = paste0(out_basedir, "prot_vs_tp_gam.txt"), quote = F, sep = "\t", row.names = FALSE)
signif <- gam_res_prot_tp[gam_res_prot_tp$gam_BH_pval < 0.05,]

cat ("Number of proteins that change significantly with time:", nrow(signif), "\n")
cat("Of them, the number of proteins with a non-linear change: ", nrow(signif[signif$gam_edf_round > 1,]), "\n")
cat("Of them, the number of proteins showing a significant association with time also in LMMs:", nrow(signif[signif$lmm_BH_pval < 0.05,]))

