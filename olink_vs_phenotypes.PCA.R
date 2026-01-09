
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
my_colors <- c("#F65C00", "#eddb6d", "#006DB3", "#00A072")
my_colors <- c("#b71f57", "#96d1aa", "#099197", "#112f2c")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12/")

library(ggplot2)
library(dplyr)
library(patchwork)
library(lubridate)
library(purrr)
library(tibble)

set.seed(123)
out_basedir <- "results12/"

d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
batch_info <- read.delim("data/batch_info.txt", as.is = T, check.names = F, sep = "\t")

d_wide$TP <- gsub(".*_","", d_wide$SampleID)
d_wide$ID <- gsub("_.*","", d_wide$SampleID)
d_wide <- d_wide %>% select(SampleID, ID, TP, everything())

batch2_shared_prots <- colnames(d_wide)[colSums(is.na(d_wide)) < 50]
d_wide_shared <- d_wide[ ,batch2_shared_prots]
d_wide_b2 <- d_wide[d_wide$SampleID %in% batch_info[batch_info$Batch == 'batch2', "SampleID"],]

dim(d_wide)
dim(d_wide_b2)
dim(d_wide_shared)


################################################################################
# PCA on proteins with missing data (nipals) on shared proteins
################################################################################

res_pca <- run_pca_nipals_per_tp(d_wide_shared, nPCs = 70)
num_pcs_80_b12 <- res_pca$num_pcs_80
pca_per_tp_b12 <- res_pca$pca_per_tp

max(as.numeric(num_pcs_80_b12))
# [1] 60
write.table(pca_per_tp_b12, file = paste0(out_basedir, "olink_batch12_shared_prot_rm_outliers_4sd.PCA.txt"), quote = F, sep = "\t", row.names = FALSE)

################################################################################
# PCA on proteins with missing data (nipals) on batch2 data
################################################################################
res_pca_b2 <- run_pca_nipals_per_tp(d_wide_b2, nPCs = 60)

num_pcs_80_b2 <- res_pca_b2$num_pcs_80
pca_per_tp_b2 <- res_pca_b2$pca_per_tp

max(as.numeric(num_pcs_80_b2))
# [1] 49

write.table(pca_per_tp_b2, file = paste0(out_basedir, "olink_batch12_only_batch2_rm_outliers_4sd.PCA.txt"), quote = F, sep = "\t", row.names = FALSE)


################################################################################
# Technical covariates
################################################################################

#
# plates, well position. 
#
# plate_pos <- read.delim("data/plate_info_batch2.txt", as.is = T, check.names = F, sep = "\t")
# plate_pos$PlateID <- as.factor(plate_pos$PlateID)
# plate_pos$Well_pos1 <- as.factor(plate_pos$Well_pos1)
# plate_pos$Well_pos2 <- as.factor(plate_pos$Well_pos2)
# plate_pos$WellID <- NULL
# plate_pos <- plate_pos %>%
#   separate_wider_delim(cols = 'SampleID', delim = "_", names = c("ID", "TP"))
# 
# tech_correl <-run_kruskal_test_each_TP(pca_per_tp, plate_pos)
# 
# bonf_cutoff <- 0.05/nrow(tech_correl)
# # [1] 0.0004166667
# tech_correl$Bonferroni_sign <- ifelse(tech_correl$KW_test_pval < bonf_cutoff, T, F)
# 
# plate_pos$TP <- as.numeric(plate_pos$TP)
# tech_merged <- full_join(plate_pos, pca_per_tp, by = c("ID", "TP"))
# 
# pdf(paste0(out_basedir, "plots/correlations_with_covariates/well_pos_vs_PCs.pdf"))
# ggplot(tech_merged, aes(x = Well_pos1, y = PC1, group = Well_pos1)) + geom_boxplot() + theme_bw() + facet_wrap(~TP)
# dev.off()
# write.table(tech_correl, file = paste0(out_basedir, "correlations_with_covariates/plate_vs_PCs_KW.txt"), quote = F, sep = "\t", row.names = FALSE)

#
# Season and storage time.
#


collect_date <- read.delim("../../phenotypes/batch12/date_of_collection_clean.txt", as.is = T, check.names = F, sep = "\t")
collect_date$date_collection <- as.Date(collect_date$date_collection, format = "%Y-%m-%d")
colnames(collect_date) <- gsub("SampleID", "ID", colnames(collect_date))
collect_date$SampleID <- paste0(collect_date$ID, "_", collect_date$visit_number)

batch_info <- read.delim("data/batch_info.txt", as.is = T, check.names = F, sep = "\t")
collect_date <- left_join(collect_date, batch_info, by = "SampleID")
colnames(collect_date) <- gsub("Batch", "batch", colnames(collect_date))

collect_date$shipment_date <- as.Date(ifelse(collect_date$batch == 'batch1', "24/09/2024", "01/09/2025"), format = "%d/%m/%Y")

collect_date$storage_months<- interval(collect_date$date_collection, collect_date$shipment_date) %/% months(1)
collect_date$storage_quarters <- round(collect_date$storage_months / 4)

colnames(collect_date) <- gsub("visit_number","TP",colnames(collect_date))
getSeason <- function(DATES) {
  WS <- as.Date("2025-12-21", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2025-3-20",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2025-6-21",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2025-9-23",  format = "%Y-%m-%d") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2025-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Autumn")))
}
collect_date$season <- getSeason(collect_date$date_collection)

collect_date$date_collection <- NULL
collect_date$shipment_date <- NULL

write.table(collect_date, file = paste0(out_basedir, "correlations_with_covariates/season_storage_time.txt"), quote = F, sep = "\t", row.names = FALSE)
collect_date$season <- factor(collect_date$season, levels = c("Winter", "Spring", "Summer", "Autumn"))
collect_date$season_num <- as.numeric(collect_date$season)

season_kw <-  run_kruskal_test_each_TP(pca_per_tp_b12, collect_date[,c("ID", "TP","season")])
storage_lm <-  run_lm_each_TP(pca_per_tp_b12, collect_date[,c("ID", "TP","storage_months")])

write.table(season_kw, file = paste0(out_basedir, "correlations_with_covariates/shared_prots_season_KW.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(storage_lm, file = paste0(out_basedir, "correlations_with_covariates/shared_prots_storage_lm_per_tp.txt"), quote = F, sep = "\t", row.names = FALSE)

tmp <- left_join(pca_per_tp_b12, collect_date[,c("ID", "TP","season", "storage_months")], by = c("ID", "TP"))
ggplot(tmp, aes(y = PC5, x = PC6, color = season)) + 
  geom_point() + 
  scale_color_manual(values = my_colors) +
  theme_minimal()

pls <- plot_PCA_boxplot(pca_per_tp_b12, collect_date[,c("ID", "storage_months")], 'PC6', 'storage_months')
pdf(paste0(out_basedir, "correlations_with_covariates/plots/shared_prots_PC_vs_storage_months.pdf"), width = 10, height = 5)
pls[[1]] + pls[[2]]
dev.off()

pls <- plot_PCA_boxplot(pca_per_tp_b12, collect_date[,c("ID", "season", "storage_quarters")], 'PC1', 'season')
pdf(paste0(out_basedir, "correlations_with_covariates/plots/shared_prots_PC_vs_season.pdf"), width = 10, height = 5)
pls[[1]] + pls[[2]]
dev.off()


# check if correcting for storage months also removes the season effect
d_wide_shared_adj <- regress_covariates_lmm(d_wide_shared, collect_date[,c("SampleID", "storage_months", "season_num")], covars_longitudinal = T)
pca_per_tp_adj <- run_pca_nipals_per_tp(d_wide_shared_adj)$pca_per_tp
season_kw_adj <-  run_kruskal_test_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","season")])
storage_lm_adj <-  run_lm_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","storage_months")])

pca_per_tp_annot <- left_join(pca_per_tp_b12, batch_info, by = "SampleID")
pca_per_tp_adj_annot <- left_join(pca_per_tp_adj, batch_info, by = "SampleID")
ggplot(pca_per_tp_annot, aes(PC1, PC2, color = geo)) + geom_point() + theme_minimal()
ggplot(pca_per_tp_adj_annot, aes(PC1, PC2, color = geo)) + geom_point() + theme_minimal()

#### Get protein vs storage association

prots <- colnames(d_wide_shared)[! colnames(d_wide_shared) %in% c("SampleID", "ID", "TP")]
tmp <- collect_date[,c("SampleID", "storage_months")]
tmp$SampleID <- gsub("^X", "", tmp$SampleID)
res_table <- data.frame(matrix(nrow = length(prots), ncol = 1))
colnames(res_table) <- "storage_pval"
row.names(res_table) <- prots
for (prot in prots){
  res_table[prot, "storage_pval"] <- simple_lmm_no_covariates(d_wide_shared, tmp, prot, "storage_months")
}

res_table$BH_pval <- p.adjust(res_table$storage_pval, method = 'BH')
res_table <- res_table[order(res_table$storage_pval),]
write.table(res_table, file = "results12/correlations_with_covariates/association_with_storage_time.txt", sep = "\t", quote = F)

################################################################################
# Main covariates
################################################################################


pheno_0 <- read.delim("../../phenotypes/batch12/questionnaire_201025_selected_visit_0.txt", sep = "\t", check.names = F, as.is = T)
pheno_0[pheno_0 == 'NA'] <- NA
pheno_0$Visit_number <- NULL
pheno_0$Code <- NULL
#pheno_0$Patient_id <- NULL
#colnames(pheno_0) <- gsub("Visit_number", "TP", colnames(pheno_0),fixed = T)
pheno_0$ID <- gsub("X", "", pheno_0$ID )
pheno_0 <- pheno_0[pheno_0$ID %in% d_wide_shared$ID,]

factor_counts <- lapply(pheno_0, function(column) {
  if (length(unique(column)) < 5) {
    return(table(column, useNA = "ifany"))
  } else {
    return(NULL)
  }
})

for (col_name in names(factor_counts)) {
  cat("Factor levels for", col_name, ":\n")
  print(factor_counts[[col_name]])
  cat("\n")
}

pheno_0$caffe_quantita_die <- NULL
pheno_0$defecazioni_media <- NULL
pheno_0$diarrea_episodi <- NULL
pheno_0$costipazione_episodi <- NULL
pheno_0$giorni_mestruazione <- NULL
pheno_0$giorni_tra_cicli <- NULL

pheno_0$Age_category <- NULL
pheno_0$BMI_ranges <- NULL

continuous_pheno <- c("Age", "BMI", "Waist_circumference", "Hip_circumference","WHR", "numero_gravidanze","altezza_cm", "peso_kg", "anni_stop_pillola")


# Kruskal Wallis test on categorical covariates
covar_kw_res <- run_kruskal_test_each_TP(pca_per_tp_b12, pheno_0[,! colnames(pheno_0) %in% continuous_pheno])
colnames(covar_kw_res)[4] <- "pval"

# Linear regression on Age and BMI
pheno_cont_0 <- pheno_0[,c("ID", continuous_pheno)]

covar_lm_res <- run_lm_each_TP(pca_per_tp, pheno_cont_0)

colnames(season_kw)[4] <- "pval"

combined_covar_res <- rbind(covar_kw_res, subset(covar_lm_res, select = -c(beta)),
                            season_kw, subset(storage_lm, select = -c(beta)) )




num_tests_per_tp <- nrow(combined_covar_res[combined_covar_res$TP == '1',])
combined_covar_res$bonf_sign <- ifelse(combined_covar_res$pval < 0.05/num_tests_per_tp, T, F)

combined_covar_res$logp <- -log10(combined_covar_res$pval)
combined_covar_res$pval <- as.numeric(combined_covar_res$pval)

combined_covar_res <- combined_covar_res[order(combined_covar_res$pval),]

write.table(combined_covar_res, file = paste0(out_basedir, 'correlations_with_covariates/PC_vs_covariates_visit_0.txt'), quote = F, sep = "\t", row.names = FALSE)

col_order <- unique(combined_covar_res$covariate)
row_order <- paste("PC", seq(1,10), sep = "")

color_palette <- colorRampPalette(c("white", "#E6F2FF", "#99CCFF", "#3399FF", "#0066CC", "#004C99"))(100)
breaks <- seq(0, max(combined_covar_res$logp, na.rm = TRUE), length.out = 100)
color_fn <- colorRamp2(breaks, color_palette)

plot_list <- list()
for (tp in 1:4){
  logp_mat <- my_pivot_wider(combined_covar_res[combined_covar_res$TP == tp,c("PC", "covariate", "logp")], row_names = 'PC', names_from = 'covariate', values_from = 'logp')
  pval_mat <- my_pivot_wider(combined_covar_res[combined_covar_res$TP == tp,c("PC", "covariate", "pval")], row_names = 'PC', names_from = 'covariate', values_from = 'pval')
  logp_mat <- logp_mat[row_order, col_order]
  pval_mat <- pval_mat[row_order, col_order]
  labels_matrix <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  labels_matrix[pval_mat < 0.05/num_tests_per_tp] <- "*"
  
  plot_list[[tp]] <- Heatmap(
    logp_mat,
    name = "-log10(p)",
    border = TRUE,
    rect_gp = gpar(col = "grey", lwd = 1),
    col = color_fn,
    cluster_rows = F,
    cluster_columns = F,
    row_names_gp = gpar(fontsize = 14),  # Row font size
    column_names_gp = gpar(fontsize = 14),  # Column font size
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(labels_matrix[i, j] != "") {
        grid.text(labels_matrix[i, j], x, y, gp = gpar(fontsize = 12))
      }
    },
    column_title = paste("Visit", tp),
    show_heatmap_legend = if(tp == 1) TRUE else FALSE  
  )
}
pdf(paste0(out_basedir, 'correlations_with_covariates/plots/corrplots_PC_vs_covariates.pdf'), width = 20, height  = 5)
draw(plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]], gap = unit(1, "cm"))
dev.off()


################################################################################
# Correct for all covariates and rerun the associations
################################################################################

clean_pheno <- read.delim("../../phenotypes/batch12/cleaned_questionnaire_251125.csv", as.is = T, sep = ",", check.names = F)
clean_pheno_subs <- unique(clean_pheno[,c("ID", "from","Age", "BMI")])
collect_date[grepl("^[0-9]",collect_date$SampleID), "SampleID"] <- paste0("X", collect_date[grepl("^[0-9]",collect_date$SampleID), "SampleID"])
collect_date[grepl("^[0-9]",collect_date$ID), "ID"] <- paste0("X", collect_date[grepl("^[0-9]",collect_date$ID), "ID"])

all_covar <- left_join(collect_date[,c("SampleID", "ID","batch","storage_months")], clean_pheno_subs, by = "ID")
write.table(all_covar, file = "results12/covariates_olink_batch12.txt", quote = F, sep = "\t", row.names = F)


#all_covar <- left_join(d_wide_shared[,c("SampleID", "ID", "TP")], collect_date[,c("SampleID", "batch","storage_months")], by = "SampleID")
#all_covar <- left_join(all_covar, pheno_0[,c("ID", "from","Age", "BMI")], by = "ID")

#covar_per_id <- all_covar[,! colnames(all_covar) %in% c("SampleID", "TP")]%>% 
#           distinct(ID, .keep_all = TRUE)
#write.table(covar_per_id, file = paste0(out_basedir, "covariates_per_id_olink_batch12.txt"), quote = F, sep = "\t", row.names = F)

all_covar[rowSums(is.na(all_covar)) > 0,]

d_wide_shared_adj <- regress_covariates_lmm(d_wide_shared, all_covar, covars_longitudinal = T)
pca_per_tp_adj <- run_pca_nipals_per_tp(d_wide_shared_adj)$pca_per_tp

season_kw_adj <-  run_kruskal_test_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","season")])
colnames(season_kw_adj)[4] <- "pval"
storage_lm_adj <-  run_lm_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","storage_months")])
storage_lm_adj$beta <- NULL

covar_kw_res_adj <- run_kruskal_test_each_TP(pca_per_tp_adj, pheno_0[,! colnames(pheno_0) %in% continuous_pheno])
colnames(covar_kw_res_adj)[4] <- "pval"
covar_lm_res_adj <- run_lm_each_TP(pca_per_tp_adj, pheno_cont_0)

combined_covar_res_adj <- rbind(season_kw_adj, storage_lm_adj, covar_kw_res_adj, subset(covar_lm_res_adj, select = -c(beta)))


################################################################################
# Functions
################################################################################


run_kruskal_test_each_TP <- function(pca_per_tp, covariates, num_pcs = 10){
  kw_res <- data.frame(matrix(ncol = 4))
  colnames(kw_res) <- c("TP", "PC", "covariate", "KW_test_pval")
  cnt <- 1
  for (tp in 1:4){
    #print(tp)
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    #pca10$TP <- NULL
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      merged <- left_join(pca10, tmp_covariates, by = c("ID", "TP"))
      #tmp_covariates$TP = NULL
    } else {
      tmp_covariates <- covariates
      merged <- left_join(pca10, tmp_covariates, by = "ID")
    }
    
    
    
    all_covs <- colnames(tmp_covariates)[! colnames(tmp_covariates) %in% c("SampleID", "ID", "TP")]
    all_pcs <- colnames(pca10)[! colnames(pca10) %in% c("SampleID", "ID", "TP")]
    for (pc in all_pcs){
      #print(pc)
      for (cov in all_covs){
        #print(cov)
        subs <- merged[,c(pc,cov) ]
        colnames(subs) <- c("PC", "cov")
        kw <- kruskal.test(PC ~ cov, data=subs)$p.value
        
        kw_res[cnt,] <- c(tp, pc, cov, kw)
        cnt <- cnt + 1
      }
    }
  }
  
  kw_res <- na.omit(kw_res) %>%
    mutate(across(-c(PC, covariate), as.numeric))   
  
  kw_res
}

run_lm_each_TP <- function(pca_per_tp, covariates, num_pcs = 10){
  lm_res <- data.frame(matrix(ncol = 5))
  colnames(lm_res) <- c("TP", "PC", "covariate", "beta", "pval")
  cnt <- 1
  for (tp in 1:4){
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    #pca10$TP <- NULL
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      merged <- left_join(pca10, tmp_covariates, pca10, by = c("ID", "TP"))
      #tmp_covariates$TP = NULL
    } else {
      tmp_covariates <- covariates
      merged <- left_join(pca10, tmp_covariates, pca10, by = "ID")
    }
    
    
    all_covs <- colnames(tmp_covariates)[! colnames(tmp_covariates) %in% c("SampleID", "ID", "TP")]
    all_pcs <- colnames(pca10)[! colnames(pca10) %in% c("SampleID", "ID", "TP")]
    
    for (pc in all_pcs){
      for (cov in all_covs){
        #cat(tp, pc, cov, "\n")
        subs <- merged[,c(pc,cov) ]
        colnames(subs) <- c("PC", "cov")
        lm_fit <- lm(PC ~ cov, data = subs)
        lm_res[cnt,] <- c(tp, pc, cov, summary(lm_fit)$coefficients['cov', 1], summary(lm_fit)$coefficients['cov', 4])
        cnt <- cnt + 1
      }
    }
  }
  
  lm_res <- na.omit(lm_res) %>%
    mutate(across(-c(PC, covariate), as.numeric))   
  
  lm_res
}

plot_PCA_boxplot <- function(pca_per_tp, covariates, pc, cov){
  
  pca_all_tps <- data.frame()
  for (tp in 1:4){
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    pca10$TP <- NULL
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      tmp_covariates$TP = NULL
    } else {
      tmp_covariates <- covariates
    }
    
    merged <- inner_join(tmp_covariates, pca10, by = c("ID"))
    pca_all_tps <- rbind(pca_all_tps, cbind(tp, merged))
  }
  

  subs <- pca_all_tps[c("PC1", "PC2", pc, cov, "tp")]
  colnames(subs) <- c("PC1", "PC2", "PC", "cov", "TP")
  if (length(unique(subs$cov)) < 5){
    p1 <- ggplot(subs, aes(x = cov, y = PC, group = cov)) + geom_boxplot() + theme_bw() + xlab(cov) + ylab(pc) + facet_wrap(~TP)
    p2 <- ggplot(pca_all_tps, aes(x = PC1, y = PC2, color = season)) + geom_point() + scale_color_manual(values = my_colors) + theme_bw()
  } else {
    p1 <- ggplot(subs, aes(x = cov, y = PC)) + geom_point() + geom_smooth(method = 'lm')+ theme_bw() + xlab(cov) + ylab(pc) + facet_wrap(~TP)
    p2 <- ggplot(subs, aes(x = PC1, y = PC2, color = cov)) + geom_point() + theme_bw() + facet_wrap(~TP)
  }
  return(list(p1,p2))
}

regression_olink_pheno <- function(joined_data, prot, ph, scale = F, kw = F){
  d <- na.omit(joined_data[,c("ID", prot, ph)])
  colnames(d) <- c("SampleID", "prot", "pheno")
  is_factor <- length(unique(d$pheno)) < 3
  d <- na.omit(d)
  if (is_factor){
    d$pheno <- as.factor(d$pheno)
  } 
  
  if(scale){
    d$prot <- scale(d$prot)
    if (! is_factor) d$pheno <- scale(d$pheno)
  }
  b <- NA
  pval <- NA
  if (! kw || ! is_factor){
    lm_fit <- lm(prot ~ pheno, data = d)
    pval <- summary(lm_fit)$coefficients[2,4]
    b <- summary(lm_fit)$coefficients[2,1]
    
  } else if (is_factor) {
    pval <- kruskal.test(prot ~ pheno, data=d)$p.value
    b <- NA
  }
  return(list("pval" = pval, "est" = b, "n" = nrow(d)))
}

run_pca_nipals_per_tp <- function(d_wide, nPCs = 70){
  num_pcs_80 <- c()
  pca_per_tp <- data.frame()
  for (tp in 1:4){
    print(tp)
    tmp_wide <- as.data.frame(subset(d_wide[d_wide$TP == tp,], select = -c(TP, SampleID, ID)))
    row.names(tmp_wide) <- d_wide[d_wide$TP == tp, "ID"]
    #tmp_wide <- na.omit(tmp_wide)
    pca <- pcaMethods::pca(tmp_wide, method = 'nipals', nPcs = nPCs, center = T, scale = 'vector')
    cumulative_variance <- cumsum(pca@R2)
    num_pcs_80 <- c(num_pcs_80, which(cumulative_variance >= 0.80)[1])
    
    pca10 <- as.data.frame(pca@scores)[,1:10] %>%
      rownames_to_column(var = 'ID') 
    pca_per_tp <- rbind(pca_per_tp, data.frame(TP = tp, pca10))
  }
  pca_per_tp$SampleID <- paste0(pca_per_tp$ID, "_", pca_per_tp$TP)
  pca_per_tp <- pca_per_tp %>% select(SampleID, ID, TP, everything())
  return(list(pca_per_tp = pca_per_tp, num_pcs_80 = num_pcs_80))
}


run_pca <- function(df, batch_info = NULL, rm_prots_with_na = T, rm_samples_with_na = F){
  if(! "SampleID" %in% colnames(df)) df$SampleID <- paste0(df$ID, "_", df$TP)
  if ("TP" %in% colnames(df)) df$TP <- NULL
  if ("ID" %in% colnames(df)) df$ID <- NULL
  row.names(df) <- df$SampleID
  df$SampleID <- NULL
  
  if(rm_samples_with_na) df <-df[rowSums(is.na(df)) == 0,]
  if(rm_prots_with_na) df <-df[,colSums(is.na(df)) == 0]
  
  pca <- stats::prcomp(df, scale = TRUE, 
                       center = TRUE)
  x <- pca$x[,1:10]
  x <- as.data.frame(x)
  if (! is.null(batch_info)){
    x2 <- cbind(x, batch_info[row.names(x), ])
    return(x2)
  }
  return(x)
}

regress_covariates_lmm <- function(data, covar_data, covars_longitudinal = T){
  
  if (!"SampleID" %in% colnames(covar_data) & covars_longitudinal) {
    covar_data <- cbind(paste0(covar_data$ID, "_",covar_data$TP), covar_data)
    colnames(covar_data)[1] <- "SampleID"
  }
  
  d_adj <- data[,c("SampleID", "ID", "TP")]
  
  data[,"TP"] <- NULL
  covar_data[,"TP"] <- NULL
  
  cat("Adjusting for covariates:\n")
  print(colnames(covar_data)[-1])
  
  cnt <- 1
  for (ph in colnames(data)[3: (ncol(data))]){
    if (ph == 'TST') next
    #print(ph)
    if (covars_longitudinal){
      covar_data$ID <- NULL
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    } else {
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "ID"))
    }
    colnames(subs)[3] <- 'pheno'
    subs$batch <- as.factor(subs$batch)
    fo_lmm <- as.formula(paste("pheno ~ ", paste(colnames(covar_data)[-1], collapse = "+"), "+ (1|ID)"))
    lmm_fit <- lmer(fo_lmm, data = subs)
    subs[,ph] <- subs$pheno - predict(lmm_fit, re.form = NA)
    
    d_adj <- left_join(d_adj, subs[, c("SampleID", ph)], by = "SampleID")
  }
  return(d_adj)
}

simple_lmm_no_covariates <- function(prot_data, pheno_data, prot, ph){
  d_subs <- inner_join(prot_data[,c("SampleID", "ID", "TP", prot)], pheno_data[,c("SampleID" ,ph)], by = c("SampleID"))
  colnames(d_subs) <- c("SampleID", "ID", "TP","prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  model <- lmer(prot ~ pheno + TP + (1|ID), data = d_subs)

  est <- summary(model)$coefficients["pheno", "Estimate"]
  pval <- summary(model)$coefficients["pheno","Pr(>|t|)"]
  return(pval)
}

