source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/utility_functions.R")

setwd("/Users/Dasha/work/Sardinia/W4H/olink/all_batches/data/")

out_basedir = "/Users/Dasha/work/Sardinia/W4H/olink/all_batches/results/"
d_wide <- fread("olink_all_batches.lod150_wide_rm_outliers_4sd.txt.gz", data.table = F)

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


#d_wide_b23$ID <- gsub("_.*", "", d_wide_b23$SampleID)
#d_wide_b23$TP <- gsub(".*_", "", d_wide_b23$SampleID)
#d_wide_b23 <- d_wide_b23 %>%
#  dplyr::select(SampleID, ID, TP, everything())

res_pca_b23 <- run_pca_nipals(d_wide_b23, nPCs = 120)
cat("number of PCs explaining 80 % of variance: ")
res_pca_b23$num_pcs_80
pca_b23 <- res_pca_b23$pca_res

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

res_pca_b123 <- run_pca_nipals(d_wide_b123, nPCs = 120)
cat("number of PCs explaining 80 % of variance: ")
res_pca_b123$num_pcs_80
# 69
pca_b123 <- res_pca_b123$pca_res


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
season_lmm <- simple_lmm(pca_b23, collect_date[,c("SampleID", "ID","batch", "city", "storage_months", "season", "median_NPX")])
write.table(season_lmm, file = paste0(out_basedir, "correlations_with_covariates/all_prots.all_visits.storage_season.LMM.txt"), quote = F, sep = "\t", row.names = FALSE)

# shared prots
season_lmm_b123 <- simple_lmm(pca_b123, collect_date[,c("SampleID", "ID","batch", "city", "storage_months", "season", "median_NPX")])
write.table(season_lmm_b123, file = paste0(out_basedir, "correlations_with_covariates/shared_prots.all_visits.storage_season.LMM.txt"), quote = F, sep = "\t", row.names = FALSE)

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
  filter(ID %in% gsub("_.*", "", d_wide_b23$SampleID)) %>%
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


selected_covars <-  names(pheno_flt)[!grepl("consumo|freq|quantita|giorni|episodes|amily_history|anni_|COVID|from",names(pheno_flt))]
selected_covars <- selected_covars[! selected_covars %in% c("Number_of_pregnancies", "defecazioni_media", "numero_figli" , "fonte", "Number_partners_three_months", "Abdominal_surgery", "BMI_ranges")]
pheno_flt_sel <- pheno_flt[,selected_covars]

name_convertion <- data.frame(it_name = c("from", "Age", "BMI", "Pill_use", "Waist_circumference", "peso_kg", "season", "storage_months", "Hip_circumference", "numero_gravidanze", "Pregnancy_category", "anni_stop_pillola", "altezza_cm", "diagnosi_covid", "alcol", "WHR", "caffe", "interventi_addome", "colester_famil_secondo", "Smoker", "t2d_famil_secondo"),
                              en_name = c("city_of_collection", "Age", "BMI", "Pill_use", "Waist_circumference", "Weight", "Season", "Storage_months", "Hip_circumference", "Number_of_pregnancies", "Pregnancy_category", "Years_since_pills", "Height", "Covid_diagnosis", "Alcohol_comsumption", "WHR", "Coffee_consumption", "Prior_abdomen_interventions", "high_cholesterol_relatives", "Smoking_status", "T2D_relatives"))

# all prots
covar_lmm <- simple_lmm(pca_b23, pheno_flt_sel)

# shared prots
covar_lmm_b123 <- simple_lmm(pca_b123, pheno_flt_sel)

# combine for all prots
combined_covar_res <- rbind(season_lmm, covar_lmm)
combined_covar_res <- combined_covar_res[combined_covar_res$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]
num_tests <- nrow(combined_covar_res)
combined_covar_res$bonf_sign <- ifelse(combined_covar_res$pval < 0.05/num_tests, T, F)

combined_covar_res$logp <- -log10(combined_covar_res$pval)
combined_covar_res$pval <- as.numeric(combined_covar_res$pval)
combined_covar_res <- combined_covar_res[order(combined_covar_res$pval),]

# combine for shared prots
combined_covar_res_b123 <- rbind(season_lmm_b123, covar_lmm_b123)
combined_covar_res_b123 <- combined_covar_res_b123[combined_covar_res_b123$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]
combined_covar_res_b123$bonf_sign <- ifelse(combined_covar_res_b123$pval < 0.05/num_tests, T, F)

combined_covar_res_b123$logp <- -log10(combined_covar_res_b123$pval)
combined_covar_res_b123$pval <- as.numeric(combined_covar_res_b123$pval)
combined_covar_res_b123 <- combined_covar_res_b123[order(combined_covar_res_b123$pval),]


# Plot heatmaps for all prots
library(colorRamp2)
library(ComplexHeatmap)
col_order <- unique(combined_covar_res$covariate)
row_order <- paste("PC", seq(1,5), sep = "")

color_palette <- colorRampPalette(c("white", "#66B3FF", "#3399FF", "#0066CC", "#004C99"))(100)
breaks <- seq(0, max(combined_covar_res$logp, na.rm = TRUE), length.out = 100)
color_fn <- colorRamp2(breaks, color_palette)

logp_mat <- my_pivot_wider(combined_covar_res[,c("PC", "covariate", "logp")], row_names = 'PC', names_from = 'covariate', values_from = 'logp')
pval_mat <- my_pivot_wider(combined_covar_res[,c("PC", "covariate", "pval")], row_names = 'PC', names_from = 'covariate', values_from = 'pval')
logp_mat <- logp_mat[row_order, col_order]
pval_mat <- pval_mat[row_order, col_order]
labels_matrix <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
labels_matrix[pval_mat < 0.05/num_tests] <- "*"

heatmap_b23 <- Heatmap(
  logp_mat,
  name = "-log10(p)",
  border = TRUE,
  rect_gp = gpar(col = "grey", lwd = 1),
  col = color_fn,
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 12),  # Row font size
  column_names_gp = gpar(fontsize = 12),  # Column font size
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(labels_matrix[i, j] != "") {
      grid.text(labels_matrix[i, j], x, y, gp = gpar(fontsize = 12))
    }
  }  
)

# shared prots
logp_mat <- my_pivot_wider(combined_covar_res_b123[,c("PC", "covariate", "logp")], row_names = 'PC', names_from = 'covariate', values_from = 'logp')
pval_mat <- my_pivot_wider(combined_covar_res_b123[,c("PC", "covariate", "pval")], row_names = 'PC', names_from = 'covariate', values_from = 'pval')
logp_mat <- logp_mat[row_order, col_order]
pval_mat <- pval_mat[row_order, col_order]
labels_matrix <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
labels_matrix[pval_mat < 0.05/num_tests] <- "*"

heatmap_b123 <- Heatmap(
  logp_mat,
  name = "-log10(p)",
  border = TRUE,
  rect_gp = gpar(col = "grey", lwd = 1),
  col = color_fn,
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 12),  # Row font size
  column_names_gp = gpar(fontsize = 12),  # Column font size
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(labels_matrix[i, j] != "") {
      grid.text(labels_matrix[i, j], x, y, gp = gpar(fontsize = 12))
    }
  }
)



pdf(paste0(out_basedir, 'correlations_with_covariates/corrplots_PC_vs_covariates.LMM.pdf'), width = 10, height  = 5)
draw(heatmap_b23 + heatmap_b123, gap = unit(1, "cm"))
dev.off()



################################################################################
# Correct for all covariates and rerun the associations
################################################################################
final_covar_pheno <- pheno_0[,c("ID", "Age", "BMI")]
all_covar <- left_join(collect_date[,c("SampleID", "ID","batch","storage_months", "season","city","median_NPX")], final_covar_pheno, by = "ID")

all_covar[rowSums(is.na(all_covar)) > 0,]
all_covar$from <- relevel(as.factor(all_covar$city), ref = "Trieste")
all_covar$city = NULL
all_covar$batch <- relevel(as.factor(all_covar$batch), ref = "batch2")

write.table(all_covar, file = paste0(out_basedir, "covariates_olink_all_batches.txt"), quote = F, sep = "\t", row.names = F)
write.table(all_covar[,c("SampleID", "ID", "TP", "Age", "BMI", "from")], file = paste0(out_basedir, "biological_covariates_olink_all_batches.txt"), quote = F, sep = "\t", row.names = F)


# check for all prots
d_wide_b23$ID <- gsub("_.*","",d_wide_b23$SampleID)
d_wide_b23$TP <- gsub(".*_","",d_wide_b23$SampleID)
d_wide_b23_adj <- regress_covariates_lmm_phase(d_wide_b23, all_covar, covars_longitudinal = T)
pca_b23_adj <- run_pca_nipals(d_wide_b23_adj,nPCs = 10)$pca_res

storage_b23_adj <-  simple_lmm(pca_b23_adj, collect_date[,c("SampleID", "ID","batch", "city", "storage_months", "season", "median_NPX")])
covar_b23_adj <- simple_lmm(pca_b23_adj, pheno_flt_sel)

combined_covar_res_adj <- rbind(storage_b23_adj, covar_b23_adj)
combined_covar_res_adj <- combined_covar_res_adj[combined_covar_res_adj$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]
num_tests_adj <- nrow(combined_covar_res_adj)
combined_covar_res_adj$bonf_sign <- ifelse(combined_covar_res_adj$pval < 0.05/num_tests_adj, T, F)

# check for shared prots
d_wide_b123$ID <- gsub("_.*","",d_wide_b123$SampleID)
d_wide_b123$TP <- gsub(".*_","",d_wide_b123$SampleID)
d_wide_b123_adj <- regress_covariates_lmm_phase(d_wide_b123, all_covar, covars_longitudinal = T)
pca_b123_adj <- run_pca_nipals(d_wide_b123_adj,nPCs = 10)$pca_res

storage_b123_adj <-  simple_lmm(pca_b123_adj, collect_date[,c("SampleID", "ID","batch", "city", "storage_months", "season", "median_NPX")])
covar_b123_adj <- simple_lmm(pca_b123_adj, pheno_flt_sel)

combined_covar_res_adj_b123 <- rbind(storage_b123_adj, covar_b123_adj)

combined_covar_res_adj_b123 <- combined_covar_res_adj_b123[combined_covar_res_adj_b123$PC %in% c("PC1", "PC2", "PC3", "PC4", 'PC5'),]
combined_covar_res_adj_b123$bonf_sign <- ifelse(combined_covar_res_adj_b123$pval < 0.05/num_tests_adj, T, F)


#### Median NPX
median_npx <- unique(b123[,c("SampleID", "sample_median_NPX")])
median_npx[grepl("^[0-9]",median_npx$SampleID), ]$SampleID <- paste0("T", median_npx[grepl("^[0-9]",median_npx$SampleID), ]$SampleID)
median_npx$ID <- gsub("_.*","", median_npx$SampleID)

median_lm_pcs <-  run_lm_each_TP(pca_per_tp, median_npx)
median_lm_pcs_adj <-  run_lm_each_TP(pca_per_tp_adj, median_npx)


################################################################################
# Create a clean dataset adjusted for technical covariates
################################################################################

d_wide$ID <- gsub("_.*","",d_wide$SampleID)
d_wide$TP <- gsub(".*_","",d_wide$SampleID)

d_wide_adj <- regress_covariates_lmm_phase(d_wide, all_covar[,c("SampleID", "batch", "storage_months", "season", "median_NPX")], covars_longitudinal = T)
write.table(d_wide_adj, file = "olink_all_batches.lod150_wide_rm_outliers_4sd.adj_technical.txt.gz", quote = F, sep = "\t", row.names = F)

d_wide_adj_2 <- regress_covariates_lmm_phase(d_wide, all_covar, covars_longitudinal = T)
write.table(d_wide_adj, file = "olink_all_batches.lod150_wide_rm_outliers_4sd.adj_all_covars.txt.gz", quote = F, sep = "\t", row.names = F)


d_wide_2 <- fread("olink_all_batches.lod150_wide.txt.gz", data.table = F)
d_wide_2$ID <- gsub("_.*","",d_wide_2$SampleID)
d_wide_2$TP <- gsub(".*_","",d_wide_2$SampleID)

d_wide_2_adj <- regress_covariates_lmm_phase(d_wide_2, all_covar[,c("SampleID", "batch", "storage_months", "season", "median_NPX")], covars_longitudinal = T)
write.table(d_wide_2_adj, file = "olink_all_batches.lod150_wide.adj_technical.txt.gz", quote = F, sep = "\t", row.names = F)

################################################################################
# Functions
################################################################################


run_pca_nipals <- function(d_wide, nPCs = 100){
  num_pcs_80 <- c()
  tmp_wide <- d_wide %>%
    column_to_rownames("SampleID") %>%
    select(-any_of(c("ID", "SampleID", "TP")))

  pca <- pcaMethods::pca(tmp_wide, method = 'nipals', nPcs = nPCs, center = T, scale = 'vector')
  cumulative_variance <- cumsum(pca@R2)
  num_pcs_80 <- which(cumulative_variance >= 0.80)[1]
  
  pca10 <- as.data.frame(pca@scores)[,1:10] %>%
    rownames_to_column(var = 'SampleID') 
  
  return(list(pca_res = pca10, num_pcs_80 = num_pcs_80))
}

simple_lmm <- function(pca_res, covar_data) {
  covar_names <- setdiff(colnames(covar_data), c("SampleID", "ID", "TP"))
  pc_names <- setdiff(colnames(pca_res), c("SampleID", "ID", "TP"))
  pca_res$ID <- gsub("_.*", "", pca_res$SampleID)
  
  lmm_results <- data.frame()
  
  for (pc in pc_names){
    for (covar in covar_names) {
      #cat(pc, covar, "\n")
      if ("SampleID" %in% colnames(covar_data)){
        d_subs <- inner_join(pca_res[,c("SampleID", pc)], covar_data[,c("SampleID", "ID",covar)], by = c("SampleID")) %>%
          rename(PC = all_of(pc), covar = all_of(covar))
      } else {
        d_subs <- inner_join(pca_res[,c("ID", pc)], covar_data[,c("ID",covar)], by = c("ID")) %>%
          rename(PC = all_of(pc), covar = all_of(covar))
      }
      model <- lmer(PC ~ covar + (1|ID), data = d_subs)
      coef <- summary(model)$coefficients[,c(1,5)]
      coef <- coef[row.names(coef) != "(Intercept)",,drop =F]
      lmm_results <- rbind(lmm_results, data.frame(PC = pc, covariate = covar, covarlevel = row.names(coef), estimate = coef[,1], pval = coef[,2]))
    }
  }
  row.names(lmm_results) <- NULL
  
  lmm_results$covariate <- ifelse(lmm_results$covarlevel == "covar", 
                                      lmm_results$covariate, 
                                      paste0(lmm_results$covariate, "_", gsub("covar","",lmm_results$covarlevel)))
  
  lmm_results$covarlevel <- NULL
  
  return(lmm_results)
}

