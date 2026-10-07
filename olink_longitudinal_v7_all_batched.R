my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/all_batches/")
source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/utility_functions.R")

library(ggplot2)
library(ggExtra)
library(rmcorr)
library(dplyr)
library(lme4)
library(grid)
library(gridExtra)
library(pheatmap)
library(corrplot)
library(patchwork)
library(ggrepel)

set.seed(123)

theme_set(theme_bw(base_size = 12))

out_basedir <- "results/primary_analysis_250926/"

d_wide <- read.delim("data/olink_all_batches.lod150_wide_rm_outliers_4sd.adj_technical.phase_avg.txt.gz", as.is = T, check.names = F, sep = "\t")
#d_wide <- read.delim("data/olink_all_batches.lod150_wide_rm_outliers_4sd.phase_avg.txt.gz", as.is = T, check.names = F, sep = "\t")

d_wide <- d_wide %>% 
  mutate( ID = sub("_.*", "", SampleID), 
          phase = relevel( factor(sub(".*_", "", SampleID), levels = c("F", "O", "EL", "LL")), ref = "F" )
  ) %>% 
  select(SampleID, ID, phase, everything())

covariates <- read.delim("results/covariates_olink_all_batches.phase_avg.txt", sep = "\t", check.names = F, as.is = T)

covariates <- covariates %>% 
  mutate( ID = sub("_.*", "", SampleID), 
          phase = relevel( factor(sub(".*_", "", SampleID), levels = c("F", "O", "EL", "LL")), ref = "F" )
  ) %>% 
  select(SampleID, ID, phase, everything())

covariates[] <- lapply(covariates, function(col) { if (length(unique(col)) < 4) factor(col) else col })

covariates$from <- relevel(as.factor(covariates$from), ref = "Trieste")
covariates$batch <- relevel(as.factor(covariates$batch), ref = "batch2")
covariates$season <- relevel(as.factor(covariates$season), ref = "Winter")

#covariate_names <- c("Age","BMI","batch", "from", "storage_months", "median_NPX", "season")
covariate_names <- c("Age","BMI","from")

all_phases <- c("F", "O", "EL", "LL")
all_prots <- colnames(d_wide)[! colnames(d_wide) %in% c("SampleID", "ID", "TP","phase")]

shared_prots <- colnames(d_wide)[colSums(is.na(d_wide)) < 20]
shared_prots <- shared_prots[! shared_prots %in% c("SampleID", "ID", "TP","phase")]
explore_ht_prots <- all_prots[! all_prots %in% shared_prots]

# Make a dataframe with proteins adjusted for all covariates per visit
joined_data <- full_join(covariates[,c("SampleID", "ID", "phase", covariate_names)], d_wide, by = c("SampleID", "ID", "phase"))
d_wide_adj_covar <- read.delim("data/olink_all_batches.lod150_wide_rm_outliers_4sd.adj_all_covars.phase_avg.txt", as.is = T, check.names = F, sep = "\t")
d_wide_adj_covar <- d_wide_adj_covar %>% 
  mutate( ID = sub("_.*", "", SampleID), 
          phase = relevel( factor(sub(".*_", "", SampleID), levels = c("F", "O", "EL", "LL")), ref = "F" )
  ) %>% 
  select(SampleID, ID, phase, everything())


# make a joined data without any covariate adjustment
d_wide_noadj <- read.delim("data/olink_all_batches.lod150_wide_rm_outliers_4sd.phase_avg.txt.gz", as.is = T, check.names = F, sep = "\t")

d_wide_noadj <- d_wide_noadj %>% 
  mutate( ID = sub("_.*", "", SampleID), 
          phase = relevel( factor(sub(".*_", "", SampleID), levels = c("F", "O", "EL", "LL")), ref = "F" )
          ) %>% 
  select(SampleID, ID, phase, everything())

joined_data_noadj <- full_join(covariates, d_wide_noadj, by = c("SampleID", "ID", "phase"))




#joined_data_adj_covar <- full_join(covariates, d_wide_adj_covar, by = c("SampleID") )

#joined_data_b23 <- joined_data[joined_data$batch %in% c("batch2", "batch3"), c("SampleID", "ID", "phase", covariate_names, explore_ht_prots)]
#joined_data_b23$batch <- droplevels(joined_data_b23$batch)
################################################################################
# variance differences for each protein
################################################################################

var_res <- compare_variances_levene(d_wide)

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

icc <- icc %>%
  mutate(
    outlier = abs(ICC - mean(ICC, na.rm = TRUE)) > 4 * sd(ICC, na.rm = TRUE) |
      abs(R2_TP - mean(R2_TP, na.rm = TRUE)) > 4 * sd(R2_TP, na.rm = TRUE)
  )

p <- ggplot(icc, aes(x = ICC, y = R2_TP)) +
  geom_point() +
  geom_text_repel(
    data = subset(icc, outlier),
    aes(label = prot),
    size = 3
  ) +
  theme_minimal() +
  xlab("ICC (variance explained by sample ID)") +
  ylab("Marginal R2 (variance explained by phase)")

# Add marginal boxplots
ggMarginal(p, type = "boxplot", margins = "both", 
           size = 5, fill = 'white', alpha = 0.7)

dev.off()


################################################################################
# Differentially expressed proteins between visits
################################################################################

# calculate mean abundance per phase for all proteins
mean_abund_per_phase_adj_covar_long <- d_wide_adj_covar %>%
  pivot_longer(cols = all_of(all_prots), 
               names_to = "Protein", 
               values_to = "Level") %>%
  group_by(Protein, phase) %>%
  summarise(mean_prot = mean(Level, na.rm = TRUE), .groups = "drop")

mean_abund_per_phase_adj_covar <- mean_abund_per_phase_adj_covar_long %>%
  pivot_wider(names_from = phase, values_from = mean_prot)
write.table(mean_abund_per_phase_adj_covar, file = paste0(out_basedir, "mean_abundance_per_phase_adj_covar.txt"), quote = F, sep = "\t", row.names = FALSE)

# Run Limma DAP analysis; handle proteins measured not in all batches: 
#    do not use batch as a covariate if it's measured in only one batch
#    for the few proteins measured in only 1 batch run the limma model on all proteins measured in that batch to enable eBayes fit work correctly, but remove those from the results

protein_batch_info <- lapply(
    all_prots, function(prot) { 
    observed <- !is.na(joined_data_noadj[[prot]]) 
    observed_batches <- unique(joined_data_noadj$batch[observed]) 
    data.frame( protein = prot, n_samples = sum(observed), n_batches = length(observed_batches), measured_batches = paste(observed_batches, collapse = ",") ) 
    }
  ) %>% 
  bind_rows()

write.table(protein_batch_info, file = paste0(out_basedir, "protein_info_num_samples_batches.txt"), quote = F, sep = "\t", row.names = FALSE)
# covariates_no_batch <- setdiff(covariate_names, "batch")
# 
# table(protein_batch_info$measured_batches)
# 
# protein_groups <- split(protein_batch_info$protein, protein_batch_info$measured_batches)
# base_covariates <- setdiff(covariate_names, "batch")
# 
# phase_comb <- t(combn(all_phases, 2))
# 
# all_results_full_adjustment <- list()
# result_index <- 1
# 
# for (i in seq_len(nrow(phase_comb))) {
# 
#   tp1 <- phase_comb[i, 1]
#   tp2 <- phase_comb[i, 2]
# 
#   comparison_results <- list()
#   group_index <- 1
# 
#   for (pattern in names(protein_groups)) {
# 
#     target_proteins <- protein_groups[[pattern]]
#     batches <- strsplit(pattern, ",", fixed = TRUE)[[1]]
# 
#     # Proteins measured in every batch used by this model
#     fit_proteins <- protein_batch_info$protein[
#       vapply(
#         protein_batch_info$measured_batches,
#         function(x) {
#           measured <- strsplit(x, ",", fixed = TRUE)[[1]]
#           all(batches %in% measured)
#         },
#         logical(1)
#       )
#     ]
# 
#     analysis_data <- joined_data_noadj %>%
#       filter(batch %in% batches) %>%
#       mutate(batch = droplevels(factor(batch)))
# 
#     model_covariates <- if (length(batches) > 1) {
#       c(base_covariates, "batch")
#     } else {
#       base_covariates
#     }
# 
#     limma_res <- run_limma(
#       joined_data = analysis_data,
#       tp1 = tp1,
#       tp2 = tp2,
#       all_prots = fit_proteins,
#       covariate_names = model_covariates
#     ) %>%
#       rownames_to_column(var = "prot") %>%
#       filter(prot %in% target_proteins) %>%
#       mutate(measured_batches = pattern)
# 
#     comparison_results[[group_index]] <- limma_res
#     group_index <- group_index + 1
#   }
# 
#   comparison_results <- bind_rows(comparison_results) %>%
#     mutate(
#       tp1 = tp1,
#       tp2 = tp2,
#       phase1_phase2 = paste(tp1, tp2, sep = "_"),
#       # Recalculate BH across all proteins for this comparison
#       adj.P.Val = p.adjust(P.Value, method = "BH")
#     )
# 
#   all_results_full_adjustment[[result_index]] <- comparison_results
#   result_index <- result_index + 1
# }
# 
# limma_res_all_full_adjustment <- bind_rows(all_results_full_adjustment)
# limma_res_all_full_adjustment$adj.P.Val <- p.adjust(limma_res_all_full_adjustment$P.Value, method = 'BH')
# limma_res_all_full_adjustment$sign <- ifelse(limma_res_all_full_adjustment$adj.P.Val < 0.05,T,F)
# 

# if we don't adjust for batch the code is simpler:

phase_comb <- t(combn(all_phases, 2))

limma_res_all <- data.frame()

for (i in seq_len(nrow(phase_comb))) {

  tp1 <- phase_comb[i, 1]
  tp2 <- phase_comb[i, 2]

  limma_res <- run_limma(joined_data, tp1, tp2, all_prots, c("from", "Age", "BMI")) %>%
    rownames_to_column(var = 'prot')
  
  limma_res_all <- rbind(limma_res_all, cbind(paste0(tp1, "_", tp2), limma_res))
}
colnames(limma_res_all)[1] <- "phase1_phase2"
limma_res_all$adj.P.Val <- p.adjust(limma_res_all$P.Value, method = 'BH')
limma_res_all$sign <- ifelse(limma_res_all$adj.P.Val < 0.05,T,F)
table(limma_res_all$sign)


# add means per phase to limma res
limma_res_with_means <- limma_res_all %>%
  separate(phase1_phase2, into = c("phase1", "phase2"), sep = "_", remove = FALSE) %>%
  left_join(mean_abund_per_phase_adj_covar_long, by = c("prot" = "Protein", "phase1" = "phase")) %>%
  rename(phase1_mean = mean_prot) %>%
  left_join(mean_abund_per_phase_adj_covar_long, by = c("prot" = "Protein", "phase2" = "phase")) %>%
  rename(phase2_mean = mean_prot) %>%
  select(-phase1, -phase2)

#all(limma_res_all[limma_res_all$adj.P.Val < 0.05, "prot"] %in% signif$prot)
write.table(limma_res_with_means, file = paste0(out_basedir, "limma_DEPs_withmeans_v2.txt"), quote = F, sep = "\t", row.names = FALSE)
#write.table(wilcox_res_all, file = paste0(out_basedir, "wilcox_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)

# cmp with previous results
limma_res_b12 <- read.delim("../batch12/results12/intensity_all_prots_220526/limma_DEPs_withmeans_fixed_paired_samples.txt", as.is = T, sep = "\t", check.names = F)
cmp <- full_join(limma_res_b12, limma_res_with_means, by = c("prot", "phase1_phase2"), suffix = c(".b12", ".b123"))
table(cmp[,c("sign.b12", "sign.b123")])

unique(cmp[,c("prot", )])
#         sign.b123
# sign.b12 FALSE  TRUE
# FALSE 14062   320
# TRUE     43   287
sign_b12 <- unique(cmp[cmp$sign.b12 == T, "prot"])
sign_b123 <- unique(cmp[cmp$sign.b123 == T, "prot"])

length(setdiff(sign_b12,sign_b123))
# [1] 16
length(setdiff(sign_b123,sign_b12))
# [1] 153

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


barplot <- ggplot(summary_data, aes(x = phase1_phase2, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Phase Comparison",
    y = "Number of DEP",
    fill = "Effect direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors[c(3,2)])

pdf(paste0(out_basedir, "plots/limma_DEPs_barplot.pdf"))
print(barplot)
dev.off()

ggplot(limma_res_all, aes(x = logFC, y = -log10(P.Value), color = sign)) + geom_point(alpha = 0.5)

# CSF3

csf3_plot <- ggplot(d_wide_adj_covar, aes(x = phase, y = CSF3, group = phase)) + 
  geom_boxplot(width = 0.3, color = my_colors[3], outliers = F) + 
  geom_jitter(alpha = 0.2, width = 0.1, color = my_colors[3]) + 
  theme_minimal() +
  ylab("CSF3 adjusted levels")

pdf(paste0(out_basedir, "plots/limma_DEPs_barplot_CSF3.pdf"), width = 3, height = 8)
print((barplot + theme(legend.position="bottom"))/csf3_plot)
dev.off()

# PROK1
prok1_plot <- ggplot(d_wide_adj_covar, aes(x = phase, y = PROK1, group = phase)) + 
  geom_boxplot(width = 0.3, color = my_colors[3], outliers = F) + 
  geom_jitter(alpha = 0.2, width = 0.1, color = my_colors[3]) + 
  theme_minimal() +
  ylab("PROK1 adjusted levels")
pdf(paste0(out_basedir, "plots/limma_DEPs_barplot_PROK1.pdf"), width = 3, height = 8)
print((barplot + theme(legend.position="bottom"))/prok1_plot)
dev.off()

limma_res_all$signif_direction <- ifelse(limma_res_all$sign, ifelse(limma_res_all$logFC < 0, -1, 1), 0)

# colored heatmap of FCs
prot_subs <- unique(limma_res_all[limma_res_all$sign == T, "prot"])
limma_heatmap2 <- my_pivot_wider(limma_res_all[limma_res_all$prot %in% prot_subs,], "prot", "phase1_phase2", "logFC")
limma_heatmap2[is.na(limma_heatmap2)] <- 0
signif_labels <- my_pivot_wider(limma_res_all[limma_res_all$prot %in% prot_subs,], "prot", "phase1_phase2", "adj.P.Val")
signif_labels <- ifelse(signif_labels < 0.05, "*", "")
signif_labels[is.na(signif_labels)] <- ""

too_strong_DAPs = c('PROK1', "CXCL13", "NOS3", "FGF21")
max_val <- max(abs(min(limma_heatmap2[!row.names(limma_heatmap2) %in% too_strong_DAPs,])), max(limma_heatmap2[!row.names(limma_heatmap2) %in% too_strong_DAPs,]))
limma_heatmap2[limma_heatmap2 > max_val] <- max_val
breaksList = seq(-max_val, max_val, by = 0.01)
if(!0 %in% breaksList) breaksList <- sort(c(breaksList, 0))

full_palette <- rev(brewer.pal(n = 11, name = "RdYlBu"))
full_palette[ceiling(length(full_palette)/2)] <- "#FFFFFF"
colorList <- colorRampPalette(full_palette)(length(breaksList))

num_k = 5
h2 <- pheatmap(limma_heatmap2, display_numbers = signif_labels, fontsize_number = 8, fontsize_row = 4,
               cluster_cols = F, cluster_rows = T,
               color = colorList, breaks = breaksList, cellwidth = 12, cellheight = 4,
               cutree_rows = num_k)

clusters <- as.data.frame(cutree(h2$tree_row,num_k)) %>%
  rownames_to_column("prot")
colnames(clusters)[2] <- "cluster"
row.names(clusters) <- clusters$prot
clusters_ordered <- clusters[h2$tree_row$labels[h2$tree_row$order], ]
# Positions where the cluster changes
gaps_row <- which(
  clusters_ordered$cluster[-1] != clusters_ordered$cluster[-nrow(clusters_ordered)]
)

# heatmap of mean values
mat_mean_abund_per_phase_adj_covar <- mean_abund_per_phase_adj_covar %>%
  filter(Protein %in% prot_subs) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

max_val <- max(abs(min(mat_mean_abund_per_phase_adj_covar[!row.names(mat_mean_abund_per_phase_adj_covar) %in% too_strong_DAPs,])), max(mat_mean_abund_per_phase_adj_covar[!row.names(mat_mean_abund_per_phase_adj_covar) %in% too_strong_DAPs,]))
mat_mean_abund_per_phase_adj_covar[mat_mean_abund_per_phase_adj_covar > max_val] <- max_val
breaksList = seq(-max_val, max_val, by = 0.01)
if(!0 %in% breaksList) breaksList <- sort(c(breaksList, 0))

full_palette <- rev(brewer.pal(n = 11, name = "RdYlBu"))
full_palette[ceiling(length(full_palette)/2)] <- "#FFFFFF"
colorList <- colorRampPalette(full_palette)(length(breaksList))

mat_mean_abund_per_phase_adj_covar <- mat_mean_abund_per_phase_adj_covar[,c("F", "O", "EL", "LL")]
h3 <- pheatmap(mat_mean_abund_per_phase_adj_covar[h2$tree_row$labels[h2$tree_row$order],], fontsize_row = 4, fontsize_number = 8, cluster_cols = F, cluster_rows = F,
               color = colorList, breaks = breaksList,cellwidth = 12, cellheight = 4, gaps_row = gaps_row)



pdf(paste0(out_basedir, "plots/limma_DEPs_heatmap.pdf"), width = 8, height = 35)
combined <- arrangeGrob(h2$gtable, h3$gtable, ncol = 2, nrow = 1)
grid::grid.draw(combined)
dev.off()


################################################################################
# Proteins that change more than by 2 sds
################################################################################
prots_strong_changes <- d_wide %>%
  select(-SampleID) %>%
  pivot_longer(-c(ID, phase), names_to = "Protein", values_to = "Val") %>%
  group_by(Protein, phase) %>%
  summarise(Mean = mean(Val, na.rm = TRUE), SD = sd(Val, na.rm = TRUE), .groups = "drop") %>%
  group_by(Protein) %>%
  filter(phase != "EL", abs(Mean - Mean[phase == "F"]) > SD[phase == "F"])
prots_strong_changes

################################################################################
# Protein vs TP GAM and LMM
################################################################################

gam_res_list <- vector("list", length(all_prots))
pb <- txtProgressBar(min = 0, max = length(all_prots), style = 3)

cnt <- 1
for (prot in all_prots){
  res_gam <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = F)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(d_wide, prot, covariates, scale = T, report_singular = T)
  gam_res_list[[i]] <- data.frame(prot = prot, gam_pval = res_gam$pval, gam_edf = res_gam$edf,gam_fval = res_gam$fval,n = res_gam$n, n_samples = res_gam$n_samples, lmm_pval = res_lmm$pval, lmm_singular = res_lmm$singular)
  cnt <- cnt + 1
  setTxtProgressBar(pb, cnt)
}
close(pb)
write.table(gam_res_prot_tp, file = paste0(out_basedir, "/prot_vs_tp_gam.txt"), quote = F, sep = "\t", row.names = F)

gam_res_prot_tp <- bind_rows(gam_res_list) %>%
  mutate(BH_pval = p.adjust(gam_pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = gam_pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval) %>%
  arrange(gam_pval) %>%
  mutate(gam_edf_round = round(gam_edf))

signif = gam_res_prot_tp[gam_res_prot_tp$BH_pval < 0.05,]
cat ("Number of proteins that change significantly with time:", nrow(signif), "\n")
cat("Of them, the number of proteins with a non-linear change: ", nrow(signif[signif$gam_edf_round > 1,]), "\n")
cat("Of them, the number of proteins showing a significant association with time also in LMMs:", nrow(signif[signif$lmm_BH_pval < 0.05,]))

write.table(gam_res_prot_tp, file = paste0(out_basedir, "/prot_vs_tp_gam.txt"), quote = F, sep = "\t", row.names = F)

####################NULL################################################################################
# Cluster longitudinal trajectories: GAM
################################################################################

n_points = 100

# Take only significant non-linear trajectories:
all_sign_prots <- signif$prot

prot_trajs <- data.frame(matrix(nrow = length(all_sign_prots) , ncol = n_points))
row.names(prot_trajs) <- all_sign_prots
colnames(prot_trajs) <- seq(1,4, length.out = n_points)

# Precompute all GAM fits once
pb <- txtProgressBar(min = 0, max = length(all_sign_prots), style = 3)
cnt <- 1
for (prot in all_sign_prots) {
  gam_fit <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = T, n_points = n_points)
  prot_trajs[prot,] <- gam_fit$predicted
  cnt <- cnt + 1
  setTxtProgressBar(pb, cnt)
}
close(pb)

write.table(prot_trajs, file = paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_", length(all_sign_prots), "_prots.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)



source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/clustering of trajectories.R")

################################################################################
# GAM  phenotype vs TP
################################################################################
#my_colors <- c("#eddb6d", "#ed9f47", "#6C9EBF", "#3B6E8F")
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

pheno_adjusted <- regress_covariates_lmm_phase(pheno, subset(covariates, select = -c(storage_months, batch)), covars_longitudinal = T)
pheno_adjusted$TP <- as.numeric(pheno_adjusted$phase)
## Hormone and lipid trajectories
all_phenos_combined <- c(all_hormones, all_phenos)

gam_res_pheno_tp <- data.frame(matrix(nrow = length(all_phenos_combined), ncol = 7))
colnames(gam_res_pheno_tp) <- c("pheno", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")

pheno_trajs <- data.frame(matrix(nrow = length(all_phenos_combined) , ncol = n_points))
row.names(pheno_trajs) <- all_phenos_combined
colnames(pheno_trajs) <- seq(1,4, length.out = n_points)

plot_list <- list()
cnt <- 1
for (ph in c(all_hormones, all_phenos)){
  print(ph)
  res_gam <- gam_prot_tp_adj_covar(pheno, ph, subset(covariates, select = -c(storage_months,batch)), scale = T, predict = T, n_points = n_points)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(pheno, ph, subset(covariates, select = -c(storage_months,batch)))
  gam_res_pheno_tp[cnt,] <- c(ph, unlist(res_gam[1:5]), res_lmm)
  cnt <- cnt + 1
  
  pheno_trajs[ph,] <- res_gam$predicted
  traj <- data.frame(TP = seq(1,4, length.out = length(res_gam$predicted)), pheno = res_gam$predicted, lower = res_gam$lower, upper = res_gam$upper)
  traj <- full_join(pheno_adjusted[,c("TP", "ID", ph)], traj, by = "TP")
  colnames(traj)[3] <- "values"
  traj$values <- scale(traj$values)
  
  ph_name = ph
  if (ph == "PROG") ph_name = "P4"
  if (ph == "17BES") ph_name = "E2"
  
  p_str = paste0("P = ", formatC(res_gam$pval, digits = 3))
  if(res_gam$pval == 0) {
    title_str <- bquote(.(ph_name)~";"~P < 2.2 %*% 10^-16)
  } else if(res_gam$pval < 0.001) {
    sci_val <- formatC(res_gam$pval, format = "e", digits = 2)
    coeff <- as.numeric(strsplit(sci_val, "e")[[1]][1])
    exp_val <- as.integer(strsplit(sci_val, "e")[[1]][2])
    title_str <- bquote(.(ph_name)~";"~P == .(coeff) %*% 10^.(exp_val))
  } else {
    title_str <- paste0(ph_name, "; P = ", formatC(res_gam$pval, digits = 3))
  }
  
  plot_list[[ph]] <- ggplot(traj) + 
    geom_line(aes(x = TP, y = pheno), color = "black", alpha = 0.7) + 
    geom_ribbon(aes(x = TP, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") +
    geom_boxplot(aes(x = TP, y = values, group = TP), width = 0.3, color = my_colors[3], outliers = F, size = 0.7) +
    geom_jitter(aes(x = TP, y = values), alpha = 0.1, width = 0.1, color = my_colors[3]) +
    labs(x = "Phase ", y = ph_name) +
    ggtitle(title_str) + 
    theme_minimal() + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  
  
}
hormone_order = c("17BES", "LH", "FSH", "PROG", "PRL")
pdf(paste0(out_basedir,"/plots/hormones_gam_withpoints3.pdf"), height = 7, width = 8)
wrap_plots(plot_list[hormone_order], ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = "bold"))
dev.off()

pdf(paste0(out_basedir,"/plots/phenotypes_gam_withpoints3.pdf"), height = 11, width = 7)
grid.arrange(grobs = plot_list[names(plot_list) %in% all_phenos], ncol = 3, nrow = 4)  
dev.off()

gam_res_pheno_tp <- na.omit(gam_res_pheno_tp) %>%
  mutate(across(-c( pheno), as.numeric)) 

gam_res_pheno_tp$gam_BH_pval <- p.adjust(gam_res_pheno_tp$gam_pval, method = 'BH')
gam_res_pheno_tp$lmm_BH_pval <- p.adjust(gam_res_pheno_tp$lmm_pval, method = 'BH')

gam_res_pheno_tp <- gam_res_pheno_tp[order(gam_res_pheno_tp$gam_pval),]
gam_res_pheno_tp$gam_edf_round <- round(gam_res_pheno_tp$gam_edf)
gam_res_pheno_tp$gam_BH_sign <- ifelse(gam_res_pheno_tp$gam_BH_pval < 0.05,T,F)

write.table(gam_res_pheno_tp, file = paste0(out_basedir, "pheno_vs_tp_gam.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(pheno_trajs, file = paste0(out_basedir, "pheno_traj_gam.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)


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

gam_res <- gam_res %>%
  na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval) 

cat("Number of BH significant associations:\n")
cat(nrow(gam_res[gam_res$BH_pval < 0.05,]), "\n")

write.table(gam_res, file = paste0(out_basedir, "prot_vs_allpheno_spline_gam.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)

b1_res <- read.delim("results12/intensity_shared_prots_261125/prot_vs_allpheno_spline_gam.shared_prots.txt", as.is = T, check.names = F, sep = "\t")
colnames(b1_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")
b2_res <- read.delim("results12/intensity_batch2_prots_261125/prot_vs_allpheno_spline_gam.batch2_prots.txt", as.is = T, check.names = F, sep = "\t")

gam_res <- rbind(b1_res, b2_res)

gam_res <- gam_res[!gam_res$pheno %in% c("TSH", "FT4", "TST"),]

gam_res <- gam_res %>%
  na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval) 

gam_res_hormones <- gam_res[gam_res$pheno %in% all_hormones,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)

gam_res_phenos <- gam_res[gam_res$pheno %in% all_phenos,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)

nrow(gam_res_hormones[gam_res_hormones$BH_pval < 0.05,])
nrow(gam_res_phenos[gam_res_phenos$BH_pval < 0.05,])

gam_res_hormones$pheno <- gsub("PROG", "P4", gam_res_hormones$pheno)
gam_res_hormones$pheno <- gsub("17BES", "E2", gam_res_hormones$pheno)

write.table(gam_res_hormones, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(gam_res_phenos, file = paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.txt"), quote = F, sep = "\t", row.names = FALSE)


### Heatmaps
gam_res_hormones <- read.delim(paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.txt"), sep = "\t", as.is = T, check.names = F)
gam_res_phenos <- read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.txt"), sep = "\t", as.is = T, check.names = F)

# signif in at least 1 hormone
prot_subs <- gam_res_hormones[gam_res_hormones$BH_pval < 0.05, "prot"]
h <- plot_association_heatmap(gam_res_hormones, prot_subs, transpose = T, cluster_cols = F, fontsize = 6)
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.BH0.05.transposed_2.pdf"), width = 3, height = 15, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

# signif in at least 2 hormones
prot_subs2 <- with(gam_res_hormones[gam_res_hormones$BH_pval < 0.05, ], 
                   names(which(table(prot) >= 2)))
h <- plot_association_heatmap(gam_res_hormones, prot_subs2, transpose = T, cluster_cols = F)
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.BH0.05.transposed.min2horm_2.pdf"), width = 4, height = 5, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

# volcano
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.volcano_no_PRL.pdf"), width = 7, height = 10, useDingbats = F)
gam_res_hormones_no_PRL <- gam_res_hormones[!(gam_res_hormones$prot == "PRL" & gam_res_hormones$pheno == 'PRL'), ]
plot_association_volcano(gam_res_hormones_no_PRL)
dev.off()

### Phenotypes

# signif for at least 1 pheno
prot_subs <- gam_res_phenos[gam_res_phenos$BH_pval < 0.05, "prot"]
col_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
gam_res_phenos$pheno <- factor(gam_res_phenos$pheno, levels = col_order)
h <- plot_association_heatmap(gam_res_phenos, prot_subs, transpose = T, 
                              cluster_cols = F, col_order = col_order, fontsize = 5)
pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.BH0.05.transposed_2.pdf"), width = 6, height = 14, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

# signif for at least 2 pheno
prot_subs2 <- with(gam_res_phenos[gam_res_phenos$BH_pval < 0.05, ], 
                   names(which(table(prot) >= 2)))
col_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
gam_res_phenos$pheno <- factor(gam_res_phenos$pheno, levels = col_order)
h <- plot_association_heatmap(gam_res_phenos, prot_subs, transpose = T, 
                              cluster_cols = F, col_order = col_order, fontsize = 6)
pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.BH0.05.transposed.min2pheno_2.pdf"), width = 6, height = 12, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.volcano.pdf"), width = 8, height = 8, useDingbats = F)
plot_association_volcano(gam_res_phenos)
dev.off()



hormone_assoc_count <- as.data.frame(table(gam_res_hormones[gam_res_hormones$BH_pval < 0.05, "pheno"]))
pdf(paste0(out_basedir, "plots/hormone_assoc_barplot.pdf"), width = 4, height = 6, useDingbats = F)

ggplot(hormone_assoc_count, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Hormone", y = "Number of associated proteins") +
  theme_minimal()
dev.off()


# Upset
library(UpSetR)
upset_matrix <- gam_res_hormones %>%
  filter(BH_pval < 0.05) %>%
  dplyr::select(prot, pheno) %>%
  distinct() %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    id_cols = prot,
    names_from = pheno,
    values_from = value,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("prot")

colnames(upset_matrix) <- rename_p4_e2(colnames(upset_matrix))
pdf(paste0(out_basedir, "plots/hormone_assoc_upset.pdf"), width = 6, height = 4, useDingbats = F,onefile = FALSE)
# Create plot
upset(upset_matrix,
      nsets = ncol(upset_matrix),
      nintersects = 20,
      order.by = "freq",
      mainbar.y.label = "Number of Proteins",
      sets.x.label = "Proteins per Hormone")


dev.off()

###
#### no phase correction
###
gam_res <- unique(read.delim(paste0(out_basedir,"prot_vs_allpheno_noTP_gam.shared_prots.txt"), as.is = T, check.names = F, sep = "\t"))
colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

gam_res <- gam_res[!gam_res$pheno %in% c("TSH", "FT4", "TST"),]

gam_res <- gam_res %>%
  na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval)

gam_res_hormones <- gam_res[gam_res$pheno %in% all_hormones,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)


gam_res_phenos <- gam_res[gam_res$pheno %in% all_phenos,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)

nrow(gam_res_hormones[gam_res_hormones$BH_pval < 0.05,])

write.table(gam_res_hormones, file = paste0(out_basedir, "prot_vs_hormones.gam.noTP.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(gam_res_phenos, file = paste0(out_basedir, "prot_vs_phenotypes.gam.noTP.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)


### Heatmaps
prot_subs <- gam_res_hormones[gam_res_hormones$BH_pval < 0.1, "prot"]
h <- plot_association_heatmap(gam_res_hormones, prot_subs, transpose = T, cutrows = 3)
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.noTP.BH0.1.transposed.pdf"), width = 5, height = 15)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

prot_subs <- gam_res_phenos[gam_res_phenos$BH_pval < 0.1, "prot"]
h <- plot_association_heatmap(gam_res_phenos, prot_subs, transpose = T)
pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.noTP.BH0.1.transposed.pdf"), width = 6, height = 10)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()


### Association with all hormones in the model

gam_prot_all_pheno_together_adj_covar(d_wide, pheno[c("SampleID", "ID","phase", all_hormones)], "PROK1", covariates)

################################################################################
# GAM hormones vs phenotypes all with all
################################################################################

gam_res_hormones_pheno <- data.frame(matrix(nrow = length(all_phenos_combined) * length(all_phenos_combined), ncol = 8))
colnames(gam_res_hormones_pheno) <- c("hormone", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

tested_comb <- c()
cnt <- 1
for (hormone in all_phenos_combined) {
  cat(hormone, "\n")
  pb <- txtProgressBar(min = 1, max = length(all_phenos), style = 3)
  i = 1
  for (ph in all_phenos_combined){
    if ( (paste(hormone, ph) %in% tested_comb) ||  (paste(ph, hormone) %in% tested_comb)) next
    if (hormone == ph) next
    res_gam <- gam_prot_pheno_adj_covar(pheno, pheno, hormone, ph, subset(covariates, select = -c(storage_months,batch)), scale = T, adjust_timepoint = 'spline')
    res_lmm <- lmm_pheno_prot_adj_covar(pheno, pheno, hormone, ph, subset(covariates, select = -c(storage_months,batch)), scale = T, adjust_timepoint = 'cubic')
    
    gam_res_hormones_pheno[cnt,] <- c(hormone, ph, unlist(res_gam), res_lmm$pval)
    cnt <- cnt + 1
    i <- i + 1
    
    tested_comb <- c(tested_comb, paste(hormone, ph))
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

gam_res_hormones_pheno <- gam_res_hormones_pheno %>%
  na.omit(gam_res_hormones_pheno) %>%
  mutate(across(-c(pheno, hormone), as.numeric)) %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval) %>%
  arrange(pval)

cat("Number of BH significant associations:\n")
cat(nrow(gam_res_hormones_pheno[gam_res_hormones_pheno$BH_pval < 0.05,]), "\n")

write.table(gam_res_hormones_pheno, file = paste0(out_basedir, "all_pheno_vs_all_pheno_spline_gam.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)




################################################################################
# GAM hormones vs proteins adjust for PRS
################################################################################
prs <- read.delim("data/merged_protein_PRS.tsv", sep = "\t", as.is = T, check.names = F)
prs[grepl("^[0-9]",prs$IID), "IID"] <- paste0("X", prs[grepl("^[0-9]",prs$IID), "IID"])

prs <- prs[, colSums(is.na(prs)) <= 100]

all_prots_prs <- intersect(all_prots, colnames(prs))
tasks <- expand.grid(prot = all_prots_prs, ph = all_phenos, stringsAsFactors = FALSE)

# 2. Setup parallel backend (use most available cores)
plan(multisession, workers = parallel::detectCores() - 1)

# 3. run calculations in parallel
results_list <- future_lapply(1:nrow(tasks), function(idx) {
  prot <- tasks$prot[idx]
  ph <- tasks$ph[idx]
  cov_prs <- left_join(covariates, prs[,c("IID", prot)], by = c("ID" = "IID"))
  colnames(cov_prs)[ncol(cov_prs)] <- "PRS"
  res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, cov_prs, 
                                      scale = T, adjust_timepoint = 'spline', 
                                      anova_pval = F, longitudinal = T)
  
  return(c(prot = prot, pheno = ph, unlist(res_gam))) 
}, future.seed = TRUE) # Important for reproducibility/random processes in GAMs

# 4. Combine results
gam_res_prs <- as.data.frame(do.call(rbind, results_list))
colnames(gam_res_prs) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples")

# 5. Clean up
gam_res_prs <- gam_res_prs %>%
  mutate(across(c(pval, estimate, SE, n, n_samples), as.numeric)) %>% 
  na.omit() %>%
  arrange(pval) %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH'))

write.table(gam_res_prs, file = paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.PRS.txt"), 
            quote = F, sep = "\t", row.names = FALSE)



gam_res_prs = read.delim(paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.PRS.txt"), sep = "\t", as.is = T, check.names = F)
gam_res_prs_ph = read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.PRS.txt"), sep = "\t", as.is = T, check.names = F)

gam_res_hormones <- read.delim(paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.txt"), sep = "\t", as.is = T, check.names = F)
gam_res_phenos <- read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.txt"), sep = "\t", as.is = T, check.names = F)

cmp_hormones <- inner_join(gam_res_hormones[,c("prot", "pheno", "pval", "estimate", "n_samples")], gam_res_prs[,c("prot", "pheno", "pval", "estimate", "n_samples")], by = c("pheno", "prot"))
cmp_pheno <- inner_join(gam_res_phenos[,c("prot", "pheno", "pval", "estimate", "n_samples")], gam_res_prs_ph[,c("prot", "pheno", "pval", "estimate", "n_samples")], by = c("pheno", "prot"))

formatC(cor(cmp_hormones$estimate.x, cmp_hormones$estimate.y, use = 'complete.obs'), digits = 2)
formatC(cor(cmp_pheno$estimate.x, cmp_pheno$estimate.y, use = 'complete.obs'), digits = 2)


p1 <- ggplot(cmp_hormones, aes(estimate.x, estimate.y)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Hormone - protein associations\ncomparison of estimates") + 
  theme_minimal()

p2 <- ggplot(cmp_hormones, aes(-log10(pval.x), -log10(pval.y) )) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Hormone - protein associations\ncomparison of -log10(P)") +
  theme_minimal()

p3 <- ggplot(cmp_pheno, aes(estimate.x, estimate.y)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Phenotype - protein associations\ncomparison of estimates") + 
  theme_minimal()

p4 <- ggplot(cmp_pheno, aes(-log10(pval.x), -log10(pval.y) )) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Phenotype - protein associations\nComparison of -log10(P)") +
  theme_minimal()

pdf(paste0(out_basedir, "plots/PRS_comparison_associations.pdf"), height = 8, width = 8, useDingbats = F)
(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 14, face="bold"))
dev.off()

write.table(cmp_hormones, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.CMP.PRS.txt"), 
            quote = F, sep = "\t", row.names = FALSE)
write.table(cmp_pheno, file = paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.CMP.PRS.txt"), 
            quote = F, sep = "\t", row.names = FALSE)


################################################################################
# GAM hormones vs proteins non linear
################################################################################

gam_res_spline <- unique(read.delim(paste0(out_basedir,"prot_vs_allpheno_spline_gam.pheno_spline.shared_prots.txt"), as.is = T, check.names = F, sep = "\t"))
#colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

gam_res_spline <- gam_res_spline[!gam_res_spline$pheno %in% c("TSH", "FT4", "TST"),]

gam_res_spline <- gam_res_spline %>%
  na.omit(gam_res_spline) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval)

gam_res_spline_hormones <- gam_res_spline[gam_res_spline$pheno %in% all_hormones,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)


gam_res_spline_phenos <- gam_res_spline[gam_res_spline$pheno %in% all_phenos,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)

nrow(gam_res_hormones[gam_res_hormones$BH_pval < 0.05,])
nrow(gam_res_hormones[gam_res_hormones$edf_round > 1 & gam_res_hormones$BH_pval < 0.05,])

write.table(gam_res_spline_hormones, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.pheno_spline.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(gam_res_spline_phenos, file = paste0(out_basedir, "prot_vs_phenotypes.gam.spline.pheno_spline.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)


write.table(gam_res_spline_hormones[gam_res_spline_hormones$edf_round > 1,], file = paste0(out_basedir, "prot_vs_hormones.gam.spline.pheno_spline.shared_prots.only_nonlinear.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(gam_res_spline_phenos[gam_res_spline_phenos$edf_round > 1,], file = paste0(out_basedir, "prot_vs_phenotypes.gam.spline.pheno_spline.shared_prots.only_nonlinear.txt"), quote = F, sep = "\t", row.names = FALSE)


cmp_h <- left_join(gam_res_hormones, gam_res_spline_hormones, by = c("pheno", "prot"))
nrow(cmp_h[cmp_h$BH_pval.x > 0.05 & cmp_h$BH_pval.y < 0.05,])





################################################################################
# Network
################################################################################


gam_res_hormones <- read.delim(paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_phenos<- read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_hormones_pheno <- read.delim(paste0(out_basedir, "all_pheno_vs_all_pheno_spline_gam.shared_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_hormones_pheno <- gam_res_hormones_pheno[gam_res_hormones_pheno$prot != gam_res_hormones_pheno$pheno,]
#colnames(gam_res_hormones_pheno)[c(1,2)] <-c("prot", "pheno")

gam_res_hormones$prot <- gsub("PRL", "PRL_olink", gam_res_hormones$prot)
gam_res_phenos$prot <- gsub("PRL", "PRL_olink", gam_res_phenos$prot)  

causal_1 <- read.delim(paste0(out_basedir, "causality/hormone_vs_phenotypes.causality_wide.linear.txt"), as.is = T, sep = "\t", check.names  = FALSE)
causal_2 <- read.delim(paste0(out_basedir, "causality/prot_vs_hormones.causality_wide.linear.txt"), as.is = T, sep = "\t", check.names  = FALSE)
causal_3 <- read.delim(paste0(out_basedir, "causality/prot_vs_phenotypes.causality_wide.linear.txt"), as.is = T, sep = "\t", check.names  = FALSE)

# rm hormone - pheno pairs whose association is not nominally significant:
sign_pairs <- c(
  paste(gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$prot, gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$pheno, sep = "_"),
  paste(gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$pheno, gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$prot, sep = "_")
)
causal_1 <- causal_1[paste(causal_1$phenotype, causal_1$protein, sep = "_") %in% sign_pairs,]

causal <- rbind(causal_1, causal_2, causal_3)

causal_links <- causal[,1:3] %>%
  separate_rows(summary, sep = ";") %>%               # one row per timepoint
  mutate(summary = str_remove(summary, "^\\d+:")) %>% # remove "1:", "2:", etc.
  filter(summary != "nothing significant") %>%        # drop empty timepoints
  rowwise() %>%
  mutate(
    pairs = if (summary == "both directions significant") {
      list(c(paste(phenotype, protein, sep = " → "),
             paste(protein, phenotype, sep = " → ")))
    } else {
      list(summary)
    }
  ) %>%
  unnest(pairs) %>%
  separate(pairs, into = c("cause", "consequence"), sep = " → ") %>%
  select(cause, consequence) %>%
  unique()

write.table(causal_links, file = paste0(out_basedir, "causality/causal_links.v2.txt"), quote = F, sep = "\t", row.names = F)


#
network_data <- rbind(gam_res_hormones[gam_res_hormones$BH_pval < 0.05, c("prot", "pheno", "estimate")],
                      gam_res_phenos[gam_res_phenos$BH_pval < 0.05, c("prot", "pheno", "estimate")],
                      gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05 , c("prot", "pheno", "estimate")])

# add hormone-pheno links that do not have a significant association, but do have a significant causal link
additional_hormones_pheno <- rbind(
  inner_join(gam_res_hormones_pheno[c("prot", "pheno", "estimate")], causal_links, by = c("prot" = "cause", "pheno" = "consequence")),
  inner_join(gam_res_hormones_pheno[c("prot", "pheno", "estimate")], causal_links, by = c("prot" = "consequence", "pheno" = "cause")))

network_data <- rbind(network_data, additional_hormones_pheno)


causal_forward <- paste(causal_links$cause, causal_links$consequence)

network_data_swapped <- network_data %>%
  mutate(
    # Check existence of both directions
    is_forward = paste(prot, pheno) %in% causal_forward,
    is_reverse = paste(pheno, prot) %in% causal_forward,
    
    # Determine if this edge has a defined direction
    has_direction = is_forward | is_reverse,
    
    # Determine how many rows we need.
    # If it is bi-directional (both forward and reverse exist), we need 2 rows.
    # Otherwise, we keep 1 row.
    n_rows = ifelse(is_forward & is_reverse, 2, 1)
  ) %>%
  # Expand the dataframe: this duplicates the rows where n_rows is 2.
  # .id = "row_id" creates a column (1 or 2) to distinguish the copies.
  uncount(n_rows, .id = "row_id") %>%
  mutate(
    # Logic to determine if we swap the nodes:
    # Swap if:
    # 1. It is strictly a reverse relationship (Reverse is T, Forward is F)
    # 2. OR It is the second copy of a bi-directional relationship (row_id == 2)
    should_swap = (is_reverse & !is_forward) | (row_id == 2),
    
    prot_new = ifelse(should_swap, pheno, prot),
    pheno_new = ifelse(should_swap, prot, pheno)
  ) %>%
  # Select and rename final columns
  dplyr::select(prot = prot_new, pheno = pheno_new, estimate, has_direction) %>%
  # Ensure no exact duplicates 
  distinct()


nodes_data <- unique(rbind(data.frame(feature = gam_res_hormones[gam_res_hormones$BH_pval < 0.05,"prot"], type = "protein"),
                           data.frame(feature = gam_res_hormones[gam_res_hormones$BH_pval < 0.05,"pheno"], type = "hormone"),
                           data.frame(feature = gam_res_phenos[gam_res_phenos$BH_pval < 0.05,"prot"], type = "protein"),
                           data.frame(feature = gam_res_phenos[gam_res_phenos$BH_pval < 0.05,"pheno"], type = "phenotype")))

nodes_data$feature <- rename_p4_e2(nodes_data$feature)
network_data_swapped$prot <- rename_p4_e2(network_data_swapped$prot)
network_data_swapped$pheno <- rename_p4_e2(network_data_swapped$pheno)


write.table(network_data_swapped, file = paste0(out_basedir, "network/network.spline.edges.causality2.with_pheno-pheno.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(nodes_data, file = paste0(out_basedir, "network/network.spline.nodes.txt"), quote = F, sep = "\t", row.names = FALSE)

