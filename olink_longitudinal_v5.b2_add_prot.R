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
library(RColorBrewer)

set.seed(123)

out_basedir <- "results12/intensity_batch2_prots_261125/"

d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.phase_avg.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$phase <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    select(SampleID, ID, phase, everything())
}

d_wide$phase <- relevel(factor(d_wide$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

covariates <- read.delim("results12/covariates_olink_batch12.phase_avg.txt", sep = "\t", check.names = F, as.is = T, colClasses = c(ID = "character"))

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
covariate_names <- c("Age","BMI", "from", "storage_months")


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
d_wide <- d_wide_b2

all_prots <- colnames(d_wide)[! colnames(d_wide) %in% c("ID", "TP","phase","SampleID")]

all_phases <- c("F", "O", "EL", "LL")

covariates$batch <- NULL

joined_data <- full_join(covariates, d_wide, by = c("SampleID", "ID", "phase"))
d_wide_adj_covar <- regress_covariates_lmm(d_wide, covariates, covars_longitudinal = T)
joined_data_adj_covar <- full_join(covariates, d_wide_adj_covar, by = c("SampleID", "ID", "phase"))

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
ggplot(icc, aes(x=ICC, y =R2_TP)) + geom_point() + theme_minimal() + xlab("ICC (variance explained by ID)") + ylab("Marginal R2 (var explained by phase)")

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

ggplot(icc, aes(x=ICC, y =R2_TP)) + geom_point() + theme_minimal() + xlab("ICC (variance explained by ID)") + ylab("Marginal R2 (var explained by phase)")

p3 <- ggplot(icc, aes(y = ICC)) + 
  geom_boxplot() + 
  theme_minimal() +
  theme(axis.text.x=element_blank())

dev.off()


################################################################################
# Differentially expressed proteins between visits
################################################################################

limma_res_all <- data.frame()
#wilcox_res_all <- data.frame()
phase_comb <- t(combn(all_phases, 2))
for (i in 1:nrow(phase_comb)){
  print(i)
  tp1 = phase_comb[i,1]
  tp2 = phase_comb[i,2]
  #tp1 <- as.character(tp1)
  #tp2 <- as.character(tp2)
  
  # limma
  limma_res <- run_limma(joined_data, tp1, tp2) %>%
    rownames_to_column(var = 'prot')
  if (nrow(limma_res[limma_res$adj.P.Val < 0.05,]) > 0) {
    limma_res_all <- rbind(limma_res_all, cbind(paste0(tp1, "_", tp2), limma_res))
  }
  
  # wilcoxon
  #wilcox_res <- run_wilcox(joined_data_adj_covar, tp1, tp2)
  #wilcox_res_all <- rbind(wilcox_res_all, wilcox_res)
}


colnames(limma_res_all)[1] <- "phase1_phase2"
limma_res_all$sign <- ifelse(limma_res_all$adj.P.Val < 0.05, T, F)

write.table(limma_res_all, file = paste0(out_basedir, "limma_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)
#write.table(wilcox_res_all, file = paste0(out_basedir, "wilcox_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)

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
  group_by(phase1_phase2, Regulation) %>%
  summarize(Count = n(), .groups = "drop") %>%
  mutate(phase1_phase2 = factor(phase1_phase2, 
                                levels = c("F_O", "F_EL", "F_LL", "O_EL", "O_LL", "EL_LL")))

pdf(paste0(out_basedir, "plots/limma_DEPs_barplot.pdf"))
ggplot(summary_data, aes(x = phase1_phase2, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Phase Comparison",
    y = "Number of DEP",
    fill = "Effect direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors[c(3,2)])
dev.off()

limma_res_all$signif_direction <- ifelse(limma_res_all$sign, ifelse(limma_res_all$logFC < 0, -1, 1), 0)
limma_heatmap <- my_pivot_wider(limma_res_all[limma_res_all$signif_direction != 0,], row_names = "prot", names_from = "phase1_phase2", values_from = "signif_direction")
limma_heatmap[is.na(limma_heatmap)] <- 0
pheatmap(limma_heatmap, color = c(my_colors[3], "white", my_colors[2]), legend = F, cluster_cols = F, filename = paste0(out_basedir, "plots/limma_DEPs_heatmap.pdf"))
dev.off()

################################################################################
# Protein vs TP GAM and LMM
################################################################################

gam_res_prot_tp <- data.frame(matrix(nrow = length(all_prots), ncol = 7))
colnames(gam_res_prot_tp) <- c("prot", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")
pb <- txtProgressBar(min = 0, max = length(all_prots), style = 3)

cnt <- 1
for (prot in all_prots){
  res_gam <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = F)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(d_wide, prot, covariates)
  gam_res_prot_tp[cnt,] <- c(prot, unlist(res_gam), res_lmm)
  cnt <- cnt + 1
  setTxtProgressBar(pb, cnt)
}
close(pb)

gam_res_prot_tp <- na.omit(gam_res_prot_tp) %>%
  mutate(across(-c( prot), as.numeric)) 

gam_res_prot_tp$gam_BH_pval <- p.adjust(gam_res_prot_tp$gam_pval, method = 'BH')
gam_res_prot_tp$lmm_BH_pval <- p.adjust(gam_res_prot_tp$lmm_pval, method = 'BH')

gam_res_prot_tp <- gam_res_prot_tp[order(gam_res_prot_tp$gam_pval),]
gam_res_prot_tp$gam_bonf_sign <- ifelse(gam_res_prot_tp$gam_pval < 0.05/60,T,F)
gam_res_prot_tp$gam_edf_round <- round(gam_res_prot_tp$gam_edf)
gam_res_prot_tp$gam_BH_sign <- ifelse(gam_res_prot_tp$gam_BH_pval < 0.05,T,F)

write.table(gam_res_prot_tp, file = paste0(out_basedir, "prot_vs_tp_gam.txt"), quote = F, sep = "\t", row.names = FALSE)
gam_res_prot_tp <- read.delim(paste0(out_basedir, "prot_vs_tp_gam.txt"), sep = "\t", as.is = T, check.names = F)
signif <- gam_res_prot_tp[gam_res_prot_tp$gam_BH_pval < 0.05,]

cat ("Number of proteins that change significantly with time:", nrow(signif), "\n")
cat("Of them, the number of proteins with a non-linear change: ", nrow(signif[signif$gam_edf_round > 1,]), "\n")
cat("Of them, the number of proteins showing a significant association with time also in LMMs:", nrow(signif[signif$lmm_BH_pval < 0.05,]))


################################################################################
#  read phenotype 
################################################################################

pheno <- read.delim("../../phenotypes/batch12/cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.phase_avg.txt", as.is = T, check.names = F, sep = "\t")            

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

pheno_adjusted <- regress_covariates_lmm_phase(pheno, subset(covariates, select = -c(storage_months)), covars_longitudinal = T)
pheno_adjusted$TP <- as.numeric(pheno_adjusted$phase)
## Hormone and lipid trajectories
all_phenos_combined <- c(all_hormones, all_phenos)


################################################################################
# GAM protein vs phenotype 
################################################################################

gam_res <- data.frame(matrix(nrow = length(all_prots) * length(all_hormones), ncol = 8))
colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

cnt <- 1
for (ph in all_hormones) {
  cat(ph, "\n")
  longitudinal <- ifelse(ph %in% c("TSH", "FT4", "TST"), F, T) 
  
  pb <- txtProgressBar(min = 1, max = length(all_prots), style = 3)
  i = 1
  for (prot in all_prots){
    res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F, longitudinal = longitudinal)
    
    res_lmm <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic', longitudinal = longitudinal)
    
    gam_res[cnt,] <- c(prot, ph, unlist(res_gam), res_lmm$pval)
    cnt <- cnt + 1
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

gam_res <- na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 
gam_res$BH_pval <- p.adjust(gam_res$pval)
gam_res$BH_lmm_pval <- p.adjust(gam_res$lmm_pval)
gam_res <- gam_res[order(gam_res$pval),]

cat("Number of BH significant associations:\n")
cat(nrow(gam_res[gam_res$BH_pval < 0.05,]), "\n")

write.table(gam_res, file = paste0(out_basedir, "prot_vs_pheno_spline_gam_adj_covar.shared_prots.all_res.BH.txt"), quote = F, sep = "\t", row.names = FALSE)

gam_res <- read.delim(paste0(out_basedir,"prot_vs_allpheno_spline_gam.batch2_prots.txt"), as.is = T, check.names = F, sep = "\t")
colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

gam_res <- na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) 

gam_res_hormones <- gam_res[gam_res$pheno %in% all_hormones,]
gam_res_phenos <- gam_res[gam_res$pheno %in% all_phenos,]

gam_res_hormones$BH_pval <- p.adjust(gam_res_hormones$pval)
gam_res_phenos$BH_pval <- p.adjust(gam_res_phenos$pval)

gam_res_hormones$BH_lmm_pval <- p.adjust(gam_res_hormones$lmm_pval)
gam_res_phenos$BH_lmm_pval <- p.adjust(gam_res_phenos$lmm_pval)

write.table(gam_res_hormones, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(gam_res_phenos, file = paste0(out_basedir, "prot_vs_phenotypes.gam.spline.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)

nrow(gam_res_hormones[gam_res_hormones$BH_pval < 0.05,])


### Heatmaps

# signif in at least 1 hormone
prot_subs <- gam_res_hormones[gam_res_hormones$BH_pval < 0.05, "prot"]
h <- plot_association_heatmap(gam_res_hormones, prot_subs, transpose = T, cluster_cols = F)
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.BH0.05.transposed.pdf"), width = 4, height = 10, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()


# volcano
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.volcano.pdf"), width = 8, height = 6, useDingbats = F)
plot_association_volcano(gam_res_hormones)
dev.off()

### Phenotypes

# signif for at least 1 pheno
prot_subs <- gam_res_phenos[gam_res_phenos$BH_pval < 0.05, "prot"]
col_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
gam_res_phenos$pheno <- factor(gam_res_phenos$pheno, levels = col_order)
h <- plot_association_heatmap(gam_res_phenos, prot_subs, transpose = T, 
                              cluster_cols = F, col_order = col_order)
pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.BH0.05.transposed.pdf"), width = 6, height = 10, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.volcano.pdf"), width = 8, height = 8, useDingbats = F)
plot_association_volcano(gam_res_phenos)
dev.off()

