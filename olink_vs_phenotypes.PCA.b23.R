source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/utility_functions.R")

setwd("/Users/Dasha/work/Sardinia/W4H/olink/all_batches/data/")

out_basedir = "/Users/Dasha/work/Sardinia/W4H/olink/all_batches/results/"
d_wide <- fread("olink_all_batches.lod150_wide_rm_outliers_4sd.txt", data.table = F)

batch_info <- read.delim("olink_all_batches_batchinfo.txt", sep = "\t", as.is = T, check.names = F)
batch_info = batch_info[batch_info$has_olink_data == T & batch_info$pass_olink_QC == T,]
batch_info$ID = gsub("_.*","", batch_info$SampleID)

# 1.
# take only batch 2 and 3 to avoid sparse data
#
d_wide_b23 <- d_wide[d_wide$SampleID %in% batch_info[batch_info$batch %in% c("batch2", "batch3"),]$SampleID, ]

# exclude proteins with many missing samples
protein_na <- data.frame(
  Protein = names(d_wide_b23)[-1],
  NA_n = colSums(is.na(d_wide_b23[, -1])),
  nonNA_n = colSums(!is.na(d_wide_b23[, -1])),
  NA_pct = colMeans(is.na(d_wide_b23[, -1])) * 100
)
high_miss_prots = protein_na[protein_na$NA_pct > 30, "Protein"]
d_wide_b23 <- d_wide_b23[,! colnames(d_wide_b23) %in% high_miss_prots]


d_wide_b23$ID <- gsub("_.*", "", d_wide_b23$SampleID)
d_wide_b23$TP <- gsub(".*_", "", d_wide_b23$SampleID)
d_wide_b23 <- d_wide_b23 %>%
  dplyr::select(SampleID, ID, TP, everything())

res_pca_b23 <- run_pca_nipals_per_tp(d_wide_b23, nPCs = 80)
cat("Max number of PCs explaining 80 % of variance: ")
max(res_pca_b3$num_pcs_80)
pca_per_tp <- res_pca_b23$pca_per_tp

# 2.
# take only all batches but 630 shared proteins
#

protein_na_b123 <- data.frame(
  Protein = names(d_wide)[-1],
  NA_n = colSums(is.na(d_wide[, -1])),
  nonNA_n = colSums(!is.na(d_wide[, -1])),
  NA_pct = colMeans(is.na(d_wide[, -1])) * 100
)

d_wide_b123 <- d_wide[,c("SampleID",protein_na_b123[protein_na_b123$nonNA_n > 650,]$Protein) ]

d_wide_b123$ID <- gsub("_.*", "", d_wide_b123$SampleID)
d_wide_b123$TP <- gsub(".*_", "", d_wide_b123$SampleID)
d_wide_b123 <- d_wide_b123 %>%
  dplyr::select(SampleID, ID, TP, everything())

res_pca_b123 <- run_pca_nipals_per_tp(d_wide_b123, nPCs = 80)
cat("Max number of PCs explaining 80 % of variance: ")
max(res_pca_b123$num_pcs_80)
# 69
pca_per_tp_b123 <- res_pca_b123$pca_per_tp


# read the questionnaire
pheno <- read.csv("../../../phenotypes/all_batches/clean_questionnaire_250526_last_fixed.csv")
pheno[,c("X.2", "X.1", "X")] <- NULL

################################################################################
# Technical covariates
################################################################################

#
# Season and storage time.
#
library(lubridate)

collect_date <- pheno[pheno$Visit_number != 0,c("Code","datavisitaodierna")]
collect_date$datavisitaodierna <-  as.Date(collect_date$datavisitaodierna, format = "%Y-%m-%d")
collect_date <- collect_date %>%
  group_by(Code) %>%
  slice_max(datavisitaodierna, n = 1, with_ties = FALSE) %>%
  ungroup()

collect_date <- collect_date %>%
  rename(SampleID = Code) %>%
  mutate(TP = as.numeric(gsub(".*_","", SampleID))) %>%
  right_join(batch_info, by = "SampleID")
  
collect_date <- collect_date %>%
  mutate(shipment_date = case_when(
    batch == 'batch1' ~ as.Date("24/09/2024", format = "%d/%m/%Y"),
    batch == 'batch2' ~ as.Date("01/09/2025", format = "%d/%m/%Y"),
    batch == 'batch3' ~ as.Date("24/08/2026", format = "%d/%m/%Y"),
    TRUE ~ as.Date(NA)
  ))
collect_date$storage_months<- interval(collect_date$datavisitaodierna, collect_date$shipment_date) %/% months(1)
collect_date$storage_quarters <- round(collect_date$storage_months / 4)

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
collect_date$season <- getSeason(collect_date$datavisitaodierna)

collect_date$date_collection <- NULL
collect_date$shipment_date <- NULL

median_npx <- unique(b123[,c("SampleID", "sample_median_NPX")])
colnames(median_npx) <- gsub("sample_", "", colnames(median_npx))
median_npx[grepl("^[0-9]",median_npx$SampleID), ]$SampleID <- paste0("T", median_npx[grepl("^[0-9]",median_npx$SampleID), ]$SampleID)
collect_date <- left_join(collect_date, median_npx, by = "SampleID")

write.table(collect_date, file = paste0(out_basedir, "correlations_with_covariates/all_batches_season_storage_time.txt"), quote = F, sep = "\t", row.names = FALSE)

collect_date <- read.delim(paste0(out_basedir, "correlations_with_covariates/all_batches_season_storage_time.txt"), sep = "\t", as.is = T, check.names = F)
collect_date$season <- factor(collect_date$season, levels = c("Winter", "Spring", "Summer", "Autumn"))
collect_date$season_num <- as.numeric(collect_date$season)


# all prots
season_kw <-  run_kruskal_test_each_TP(pca_per_tp, collect_date[,c("ID", "TP","season", "batch")])
storage_lm <-  run_lm_each_TP(pca_per_tp, collect_date[,c("ID", "TP","storage_months","median_NPX")])

write.table(season_kw, file = paste0(out_basedir, "correlations_with_covariates/all_prots_season_KW.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(storage_lm, file = paste0(out_basedir, "correlations_with_covariates/all_prots_storage_lm_per_tp.txt"), quote = F, sep = "\t", row.names = FALSE)

# shared prots
season_kw_b123 <-  run_kruskal_test_each_TP(pca_per_tp_b123, collect_date[,c("ID", "TP","season", "batch")])
storage_lm_b123 <-  run_lm_each_TP(pca_per_tp_b123, collect_date[,c("ID", "TP","storage_months","median_NPX", "season_num")])

write.table(season_kw_b123, file = paste0(out_basedir, "correlations_with_covariates/shared_prots_season_KW.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(storage_lm_b123, file = paste0(out_basedir, "correlations_with_covariates/shared_prots_storage_lm_per_tp.txt"), quote = F, sep = "\t", row.names = FALSE)


################################################################################
# Main covariates
################################################################################

pheno_0 <- pheno[pheno$Visit_number == '0',]
pheno_0[pheno_0 == 'NA'] <- NA
pheno_0[pheno_0 == ''] <- NA

pheno_0$Visit_number = NULL

pheno_0$T2D_relatives <- ifelse(pheno_0$Family_history_T2D == 'Yes' | pheno_0$Second_degree_family_history_T2D == 'Yes', 'Yes', 'No')
pheno_0$hypercholesterolemia_relatives <- ifelse(pheno_0$Family_history_hypercholesterolemia == 'Yes' | pheno_0$Second_degree_family_history_hypercholesterolemia == 'Yes', 'Yes', 'No')

pheno_clean <- pheno_0 %>%
  filter(ID %in% d_wide_b23$ID) %>%
  mutate(across(
    everything(),
    ~ type.convert(.x, as.is = TRUE)
  ))
pheno_clean <- pheno_clean %>%
  mutate(across(
    where(~ n_distinct(.x, na.rm = TRUE) <= 3),
    as.factor
  ))
pheno_flt <- pheno_clean %>%
  select(
    where(~ mean(is.na(.x)) <= 0.90 &
            (!is.factor(.x) ||
               (n_distinct(.x, na.rm = TRUE) > 1 &&
                  all(table(.x, useNA = "no") >= 10))))
  )
pheno_flt <- pheno_flt %>%
  select(-all_of(setdiff(names(pheno_flt)[sapply(pheno_flt, is.character)], "ID")))


selected_covars <-  names(pheno_flt)[!grepl("consumo|freq|quantita|giorni|episodes|amily_history|anni_",names(pheno_flt))]
selected_covars <- selected_covars[! selected_covars %in% c("Number_of_pregnancies", "defecazioni_media", "numero_figli" , "fonte", "Number_partners_three_months", "Abdominal_surgery", "BMI_ranges")]
pheno_flt_sel <- pheno_flt[,selected_covars]

name_convertion <- data.frame(it_name = c("from", "Age", "BMI", "Pill_use", "Waist_circumference", "peso_kg", "season", "storage_months", "Hip_circumference", "numero_gravidanze", "Pregnancy_category", "anni_stop_pillola", "altezza_cm", "diagnosi_covid", "alcol", "WHR", "caffe", "interventi_addome", "colester_famil_secondo", "Smoker", "t2d_famil_secondo"),
                              en_name = c("city_of_collection", "Age", "BMI", "Pill_use", "Waist_circumference", "Weight", "Season", "Storage_months", "Hip_circumference", "Number_of_pregnancies", "Pregnancy_category", "Years_since_pills", "Height", "Covid_diagnosis", "Alcohol_comsumption", "WHR", "Coffee_consumption", "Prior_abdomen_interventions", "high_cholesterol_relatives", "Smoking_status", "T2D_relatives"))


continuous_pheno = names(pheno_flt_sel)[sapply(pheno_flt_sel, is.numeric)]

# all prots
covar_kw_res <- run_kruskal_test_each_TP(pca_per_tp, pheno_flt_sel[,! colnames(pheno_flt_sel) %in% continuous_pheno])

# shared prots
covar_kw_res_b123 <- run_kruskal_test_each_TP(pca_per_tp_b123, pheno_flt_sel[,! colnames(pheno_flt_sel) %in% continuous_pheno])


# Linear regression on continuous phenotypes
pheno_cont_0 <- pheno_flt_sel[,c("ID", continuous_pheno)]
# all prots
covar_lm_res <- run_lm_each_TP(pca_per_tp, pheno_cont_0)
# shared prots
covar_lm_res_b123 <- run_lm_each_TP(pca_per_tp_b123, pheno_cont_0)


# combine for all prots
combined_covar_res <- rbind(covar_kw_res, subset(covar_lm_res, select = -c(beta)),
                            season_kw, subset(storage_lm, select = -c(beta)))

combined_covar_res <- combined_covar_res[combined_covar_res$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]
#combined_covar_res <- combined_covar_res[! combined_covar_res$covariate == "Site",]

num_tests_per_tp <- nrow(combined_covar_res[combined_covar_res$TP == '1',])
combined_covar_res$bonf_sign <- ifelse(combined_covar_res$pval < 0.05/num_tests_per_tp, T, F)

combined_covar_res$logp <- -log10(combined_covar_res$pval)
combined_covar_res$pval <- as.numeric(combined_covar_res$pval)

combined_covar_res <- combined_covar_res[order(combined_covar_res$pval),]


# combine for shared prots
combined_covar_res_b123 <- rbind(covar_kw_res_b123, subset(covar_lm_res_b123, select = -c(beta)),
                            season_kw_b123, subset(storage_lm_b123, select = -c(beta)))

combined_covar_res_b123 <- combined_covar_res_b123[combined_covar_res_b123$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]
#combined_covar_res <- combined_covar_res[! combined_covar_res$covariate == "Site",]

num_tests_per_tp_b123 <- nrow(combined_covar_res_b123[combined_covar_res_b123$TP == '1',])
combined_covar_res_b123$bonf_sign <- ifelse(combined_covar_res_b123$pval < 0.05/num_tests_per_tp_b123, T, F)

combined_covar_res_b123$logp <- -log10(combined_covar_res_b123$pval)
combined_covar_res_b123$pval <- as.numeric(combined_covar_res_b123$pval)

combined_covar_res_b123 <- combined_covar_res_b123[order(combined_covar_res_b123$pval),]


# pheno_visits <- pheno[pheno$Visit_number != '0',]
# 
# pheno_visits <- pheno_visits %>%
#   select(where(~ !all(is.na(.x))))
# 
# colnames(pheno_visits)[! colnames(pheno_visits) %in% colnames(pheno_0)]



# Plot heatmaps for all prots
library(colorRamp2)
library(ComplexHeatmap)
col_order <- unique(combined_covar_res$covariate)
row_order <- paste("PC", seq(1,5), sep = "")

color_palette <- colorRampPalette(c("white", "#E6F2FF", "#99CCFF", "#3399FF", "#0066CC", "#004C99"))(100)
color_palette <- colorRampPalette(c("white", "#66B3FF", "#3399FF", "#0066CC", "#004C99"))(100)
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
pdf(paste0(out_basedir, 'correlations_with_covariates/all_prots_corrplots_PC_vs_covariates.pdf'), width = 20, height  = 5)
draw(plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]], gap = unit(1, "cm"))
dev.off()



# Plot heatmaps for shared prots
col_order <- unique(combined_covar_res_b123$covariate)
row_order <- paste("PC", seq(1,5), sep = "")

color_palette <- colorRampPalette(c("white", "#E6F2FF", "#99CCFF", "#3399FF", "#0066CC", "#004C99"))(100)
color_palette <- colorRampPalette(c("white", "#66B3FF", "#3399FF", "#0066CC", "#004C99"))(100)
breaks <- seq(0, max(combined_covar_res_b123$logp, na.rm = TRUE), length.out = 100)
color_fn <- colorRamp2(breaks, color_palette)

plot_list <- list()
for (tp in 1:4){
  logp_mat <- my_pivot_wider(combined_covar_res_b123[combined_covar_res_b123$TP == tp,c("PC", "covariate", "logp")], row_names = 'PC', names_from = 'covariate', values_from = 'logp')
  pval_mat <- my_pivot_wider(combined_covar_res_b123[combined_covar_res_b123$TP == tp,c("PC", "covariate", "pval")], row_names = 'PC', names_from = 'covariate', values_from = 'pval')
  logp_mat <- logp_mat[row_order, col_order]
  pval_mat <- pval_mat[row_order, col_order]
  labels_matrix <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  labels_matrix[pval_mat < 0.05/num_tests_per_tp_b123] <- "*"
  
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
pdf(paste0(out_basedir, 'correlations_with_covariates/shared_prots_corrplots_PC_vs_covariates.pdf'), width = 20, height  = 5)
draw(plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]], gap = unit(1, "cm"))
dev.off()


################################################################################
# Correct for all covariates and rerun the associations
################################################################################
final_covar_pheno <- pheno_0[,c("ID", "Age", "BMI", "Height_cm")]
all_covar <- left_join(collect_date[,c("SampleID", "ID","batch","storage_months", "season", "city","median_NPX")], final_covar_pheno, by = "ID")

all_covar[rowSums(is.na(all_covar)) > 0,]
all_covar$from <- relevel(as.factor(all_covar$city), ref = "Trieste")
all_covar$city = NULL
all_covar$batch <- relevel(as.factor(all_covar$batch), ref = "batch2")

write.table(all_covar, file = paste0(out_basedir, "covariates_olink_all_batches.extended.txt"), quote = F, sep = "\t", row.names = F)

# check for all prots
d_wide_b23_adj <- regress_covariates_lmm_phase(d_wide_b23, all_covar, covars_longitudinal = T)
pca_per_tp_adj <- run_pca_nipals_per_tp(d_wide_b23_adj,nPCs = 10)$pca_per_tp

season_kw_adj <-  run_kruskal_test_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","season")])
storage_lm_adj <-  run_lm_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","storage_months","median_NPX")])

covar_kw_res_adj <- run_kruskal_test_each_TP(pca_per_tp_adj, pheno_flt_sel[,! colnames(pheno_flt_sel) %in% continuous_pheno])
covar_lm_res_adj <- run_lm_each_TP(pca_per_tp_adj, pheno_cont_0)

combined_covar_res_adj <- rbind(covar_kw_res_adj, subset(covar_lm_res_adj, select = -c(beta)),
                            season_kw_adj, subset(storage_lm_adj, select = -c(beta)))

combined_covar_res_adj <- combined_covar_res_adj[combined_covar_res_adj$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]


num_tests_per_tp_adj <- nrow(combined_covar_res_adj[combined_covar_res_adj$TP == '1',])
combined_covar_res_adj$bonf_sign <- ifelse(combined_covar_res_adj$pval < 0.05/num_tests_per_tp_adj, T, F)

# check for shared prots
d_wide_b123_adj <- regress_covariates_lmm_phase(d_wide_b123, all_covar, covars_longitudinal = T)
pca_per_tp_adj <- run_pca_nipals_per_tp(d_wide_b123_adj,nPCs = 10)$pca_per_tp

season_kw_adj <-  run_kruskal_test_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","season")])
storage_lm_adj <-  run_lm_each_TP(pca_per_tp_adj, collect_date[,c("ID", "TP","storage_months", "season_num")])

covar_kw_res_adj <- run_kruskal_test_each_TP(pca_per_tp_adj, pheno_flt_sel[,! colnames(pheno_flt_sel) %in% continuous_pheno])
covar_lm_res_adj <- run_lm_each_TP(pca_per_tp_adj, pheno_cont_0)

combined_covar_res_adj_b123 <- rbind(covar_kw_res_adj, subset(covar_lm_res_adj, select = -c(beta)),
                                season_kw_adj, subset(storage_lm_adj, select = -c(beta)))

combined_covar_res_adj_b123 <- combined_covar_res_adj_b123[combined_covar_res_adj_b123$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]
num_tests_per_tp_adj <- nrow(combined_covar_res_adj_b123[combined_covar_res_adj_b123$TP == '1',])
combined_covar_res_adj_b123$bonf_sign <- ifelse(combined_covar_res_adj_b123$pval < 0.05/num_tests_per_tp_adj, T, F)


#### Median NPX
median_npx <- unique(b123[,c("SampleID", "sample_median_NPX")])
median_npx[grepl("^[0-9]",median_npx$SampleID), ]$SampleID <- paste0("T", median_npx[grepl("^[0-9]",median_npx$SampleID), ]$SampleID)
median_npx$ID <- gsub("_.*","", median_npx$SampleID)

median_lm_pcs <-  run_lm_each_TP(pca_per_tp, median_npx)
median_lm_pcs_adj <-  run_lm_each_TP(pca_per_tp_adj, median_npx)

################################################################################
# Functions
################################################################################


run_kruskal_test_each_TP <- function(pca_per_tp, covariates, num_pcs = 10){
  kw_res <- data.frame(matrix(ncol = 4))
  colnames(kw_res) <- c("TP", "PC", "covariate", "pval")
  cnt <- 1
  for (tp in 1:4){
    #print(tp)
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,1:(num_pcs+3)]
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
    cat("10 PCs explain", cumulative_variance[10], " of variance\n")
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

