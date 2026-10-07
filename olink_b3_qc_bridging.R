
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(patchwork)
library(data.table)
library(tibble)
set.seed(123)
setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch3/data/")

outfile_base = "olink_batch3."

###
### Batch3 read and QC
###
# read both CVD and INF panel files removing failed samples
b3 <- read_NPX("Q-08153_NPX_2026-09-10.parquet") 

# add a batch3 suffix to control samples
b3$SampleID <- ifelse(b3$SampleType != "SAMPLE", paste0(b3$SampleID, "_batch3"), b3$SampleID)

b3$geo <- "Trieste"
b3[startsWith(b3$SampleID, "S"), "geo"] <- "Sardinia"
b3[startsWith(b3$SampleID, "B"), "geo"] <- "Bologna"

# Generate check log
check_log_b3 <- OlinkAnalyze::check_npx(
  df = b3
)

b3_clean <- OlinkAnalyze::clean_npx(
  df = b3,
  # remove control assays
  remove_control_assay = TRUE,
  remove_control_sample = TRUE,
  # remove assay and QC warnings
  remove_assay_warning = TRUE,
  remove_qc_warning = TRUE,
  check_log = check_log_b3
)

# Generate check log on cleaned data
check_log_b3_clean <- OlinkAnalyze::check_npx(
  df = b3_clean
)

length(unique(b3_clean$SampleID))
# [1] 198
length(unique(b3_clean$Assay))
# [1] 5387

# add LOD
b3_clean <- olink_lod(b3_clean,
                     lod_method = "FixedLOD",
                     lod_file_path = "../../batch12/data/Explore HT_Fixed LOD.csv",
                     check_log = check_log_b3_clean)

get_failed_assays <- function(df, threshold = 0.10) {
  study <- df %>%
    filter(SampleType == "SAMPLE")
  
  n_samples <- n_distinct(study$SampleID)
  
  study %>%
    group_by(OlinkID) %>%
    summarise(
      n_observed = n_distinct(SampleID[!is.na(NPX)]),
      missing_rate = 1 - n_observed / n_samples,
      .groups = "drop"
    ) %>%
    filter(missing_rate > threshold)
}
get_failed_samples <- function(df, threshold = 0.10) {
  study <- df %>%
    filter(SampleType == "SAMPLE")
  
  n_assays <- n_distinct(study$OlinkID)
  
  study %>%
    group_by(SampleID) %>%
    summarise(
      n_observed = n_distinct(OlinkID[!is.na(NPX)]),
      missing_rate = 1 - n_observed / n_assays,
      .groups = "drop"
    ) %>%
    filter(missing_rate > threshold)
}

## Remove assays with high missing rate > 10%
b3_high_missing_assays <- get_failed_assays(b3_clean)
cat("Proteins with missing rate > 10%:\n")
b3_high_missing_assays
b3_clean <- b3_clean %>%
  filter(!OlinkID %in% b3_high_missing_assays$OlinkID)

b3_high_missing_samples <- get_failed_samples(b3_clean)
cat("Samples with missing rate > 10%:\n")
b3_high_missing_samples
b3_clean <- b3_clean %>%
  filter(!SampleID %in% b3_high_missing_samples$SampleID)


cat("Batch 3: number of samples: ", length(unique(b3_clean$SampleID)), "for", length(unique(gsub("_.*","",b3_clean[b3_clean$SampleType == "SAMPLE",]$SampleID))), "unique women\n")
#Batch 3: number of samples:  197 for 75 unique women
cat("Batch 3: number of proteins:", length(unique(b3_clean$Assay)), "\n")
#Batch 3: number of proteins: 5387 

###
### Read Batch2
###
b2 <- read_NPX("../../batch12/data/Q-08152_NPX_2025-09-17.parquet")

check_log_b2 <- OlinkAnalyze::check_npx(
  df = b2
)

# Add LOD
b2 <- OlinkAnalyze::olink_lod(b2,
                              lod_method = "FixedLOD",
                              lod_file_path = "../../batch12/data/Explore HT_Fixed LOD.csv",
                              check_log = check_log_b2)



b2$geo <- "Trieste"
b2[startsWith(b2$SampleID, "S"), "geo"] <- "Sardinia"
b2[startsWith(b2$SampleID, "B"), "geo"] <- "Bologna"

b2_clean <- OlinkAnalyze::clean_npx(
  df = b2,
  # remove control assays
  remove_control_assay = TRUE,
  remove_control_sample = TRUE,
  # remove assay and QC warnings
  remove_assay_warning = TRUE,
  remove_qc_warning = TRUE,
  check_log = check_log_b2
)

# Generate check log on cleaned data
check_log_b2_clean <- OlinkAnalyze::check_npx(
  df = b2_clean
)

length(unique(b2_clean$SampleID))
# [1] 344
length(unique(b2_clean$Assay))
# [1] 5415

## Remove assays with high missing rate > 10%
b2_high_missing_assays <- get_failed_assays(b2_clean)
cat("Proteins with missing rate > 10%:\n")
b2_high_missing_assays
b2_clean <- b2_clean %>%
  filter(!OlinkID %in% b2_high_missing_assays$OlinkID)

b2_high_missing_samples <- get_failed_samples(b2_clean)
cat("Samples with missing rate > 10%:\n")
b2_high_missing_samples
b2_clean <- b2_clean %>%
  filter(!SampleID %in% b2_high_missing_samples$SampleID)


cat("Batch 2: number of samples: ", length(unique(b2_clean$SampleID)), "for", length(unique(gsub("_.*","",b2_clean$SampleID))), "unique women\n")
# Batch 2: number of samples:  344 for 130 unique women
cat("Batch 2: number of proteins:", length(unique(b2_clean$Assay)), "\n")
# Batch 2: number of proteins: 5407 


###
### Look at the overlap between batches:
###

cat("Protein with plate control normalization: ", as.character(unique(b3_clean[b3_clean$Normalization == "Plate control", "Assay"])), "\n")

overlapping_samples <- intersect(b3_clean$SampleID, b2_clean$SampleID)
cat ("Number of overlapping samples usable for bridging:", length(overlapping_samples), "\n")

overlapping_prots <- intersect(unique(b3_clean$OlinkID), b2_clean$OlinkID)
cat ("Number of overlapping proteins:", length(overlapping_prots), "\n")

###
### PCA plots before bridging
###


# Plot PCA to see if bridging samples are outliers
b3_before_br <- b3_clean |>
  mutate(Type = if_else(SampleID %in% overlapping_samples,
                        paste0("batch3_bridge"),
                        paste0("batch3_sample")),
         Batch = "batch3")

b2_before_br <- b2_clean |>
  dplyr::mutate(Type = if_else(SampleID %in% overlapping_samples,
                               paste0("batch2_bridge"),
                               paste0("batch2_sample")),
                Batch = "batch2")

b3_log <- check_npx(b3_before_br)
b2_log <- check_npx(b2_before_br)

### PCA plot and check 
pca_b3 <- olink_pca_plot(df = b3_before_br,
                         color_g = "Type",
                         quiet = TRUE,
                         check_log = b3_log) 
pca_b3_geo <- olink_pca_plot(df = b3_before_br,
                         color_g = "geo",
                         quiet = TRUE,
                         check_log = b3_log) 
pca_b2 <- olink_pca_plot(df = b2_before_br,
                         color_g = "Type",
                         quiet = TRUE,
                         check_log = b2_log)
pca_b2_geo <- olink_pca_plot(df = b2_before_br,
                             color_g = "geo",
                             quiet = TRUE,
                             check_log = b2_log)

outliers = pca_b3[[1]]@data[pca_b3[[1]]@data$PCY > 0.1,]$SampleID

npx_df <- bind_rows(b3_before_br[b3_before_br$OlinkID %in% overlapping_prots,], 
                    b2_before_br[b2_before_br$OlinkID %in% overlapping_prots,])



###
### Run bridging
###



overlap_samples_list <- list("DF1" = overlapping_samples,
                             "DF2" = overlapping_samples)

npx_br_data <- olink_normalization_bridge(project_1_df = b3_before_br,
                                          project_2_df = b2_before_br,
                                          bridge_samples = overlap_samples_list,
                                          project_1_name = "batch3",
                                          project_2_name = "batch2",
                                          project_ref_name = "batch2",
                                          project_1_check_log = b3_log,
                                          project_2_check_log = b2_log)


# Extract from bridged df only the batch 3 samples, and add batch2 from the original df, because otherwise only shared proteins are retained
npx_br_data_b3 <- npx_br_data[npx_br_data$Batch == 'batch3',]

# Check adjustment factor
ggplot(npx_br_data_b3, aes(x = Adj_factor)) + geom_histogram(bins = 50) + theme_minimal()
ggplot(npx_br_data_b3, aes(x = Adj_factor)) + geom_histogram(bins = 50) + theme_minimal() + coord_cartesian(ylim = c(0, 10000))

# proteins with outlying adj factor
adj_factor_outiers = unique(npx_br_data_b3[abs(npx_br_data_b3$Adj_factor) > 2, ]$Assay)
adj_factor_outiers
#[1] "CLBA1"  "DIXDC1" "ECD"    "KPNA7"  "ODAD4"  "SLC1A4"

# plot violin plots for outliers
bridge_samples <- npx_df |>
  dplyr::filter(.data[["SampleID"]] %in% .env[["overlapping_samples"]]) |>
  dplyr::filter(.data[["Assay"]] %in% adj_factor_outiers) |>
  dplyr::mutate(Assay_OID = .data[["Assay"]])

# Generate violin plot for proteins with extreme adj factor
npx_df |>
  dplyr::filter( .data[["Assay"]] %in% adj_factor_outiers) |>
  dplyr::mutate(Assay_OID = .data[["Assay"]]) |>
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["Batch"]], y = .data[["NPX"]])) +
  ggplot2::geom_violin(mapping = ggplot2::aes(fill = .data[["Batch"]])) +
  ggplot2::geom_point(data = bridge_samples,  position = ggplot2::position_jitter(width = 0.1)) +
  ggplot2::theme(legend.position = "none") +
  OlinkAnalyze::set_plot_theme() +
  ggplot2::facet_wrap(facets = ggplot2::vars(.data[["Assay_OID"]]), nrow = 2)

npx_df$PanelDataArchiveVersion = NULL
npx_df$ExploreVersion = NULL

# remove proteins with a large ads adjustment factor
npx_br_b3_clean <- npx_br_data_b3 %>%
  filter(!SampleID %in% overlapping_samples) %>%
  filter(!Assay %in% adj_factor_outiers)

length(unique(npx_br_b3_clean$SampleID))
# Here combine again with batch 2 data
shared_colnames <- intersect(colnames(b2_before_br), colnames(npx_br_b3_clean))
npx_br_clean <- rbind(b2_before_br[,shared_colnames], npx_br_b3_clean[,shared_colnames])
npx_br_clean_log <- check_npx(df = npx_br_clean)

###
### QC plots post bridging
###

# QC plot
library("ggrepel") 
olink_qc_plot(
  df = npx_br_clean,
  color_g = "Batch",
  label_outliers = T,
  check_log = npx_br_clean_log,IQR_outlierDef = 4
)
olink_qc_plot(
  df = npx_br_clean,
  color_g = "geo",
  label_outliers = T,
  check_log = npx_br_clean_log,IQR_outlierDef = 4
)


median_iqr_outliers <- c("B042_3", "B042_4", "B042_1", "B037_1", "S034_3")

# PCA plots
pca1 <- olink_pca_plot(
  df = npx_br_clean,
  color_g = "Type",
  byPanel = TRUE,
  check_log = npx_br_clean_log, outlierDefY = 4, outlierDefX = 4,outlierLines = T
)
pca_data <-pca1$Explore_HT@data
outliers_after_br <- pca_data[pca_data$Outlier == 1, "SampleID"]

pca2 <- olink_pca_plot(
  df = npx_br_clean[!npx_br_clean$SampleID %in% outliers_after_br,],
  color_g = "Type",
  byPanel = TRUE,
  check_log = npx_br_clean_log, outlierDefY = 4, outlierDefX = 4,outlierLines = T
)
pca2_data <-pca2$Explore_HT@data
outliers_after_br2 <- pca2_data[pca2_data$Outlier == 1, "SampleID"]
outliers_after_br2

npx_br_clean_flt <- npx_br_clean[!npx_br_clean$SampleID %in% c(outliers_after_br, outliers_after_br2),]

pca3 <- olink_pca_plot(
  df = npx_br_clean_flt,
  color_g = "Type",
  byPanel = TRUE,
  check_log = npx_br_clean_log,  outlierDefY = 4, outlierDefX = 4,outlierLines = T
)

olink_pca_plot(
  df = npx_br_clean_flt,
  color_g = "geo",
  byPanel = TRUE, 
  check_log = npx_br_clean_log,  outlierDefY = 4, outlierDefX = 4,outlierLines = T
)

if (requireNamespace(package = "ggrepel", quietly = TRUE)) {
  OlinkAnalyze::olink_qc_plot(
    df = npx_br_clean_flt,
    color_g = "Batch",
    label_outliers = T,
    check_log = npx_br_clean_log,IQR_outlierDef = 4
  )
}

# Write the resulting files
fwrite(npx_br_data, file = paste0(outfile_base, "b2+b3.bridged_raw.txt.gz"), quote = F, sep = "\t", row.names = F)
fwrite(npx_br_clean_flt, file = paste0(outfile_base, "b2+b3.filtered.txt"), quote = F, sep = "\t", row.names = F)

rm(npx_br_data)


###
### Combine with batch 1
###
npx_br_clean_flt <- fread("olink_batch3.b2+b3.filtered.txt", data.table = F)

b12 <- fread("../../batch12/data/olink_batch12.intensity.bridged_all_proteins.txt", data.table = F)
b2_samples <- unique(intersect(b12$SampleID, npx_br_clean_flt$SampleID))
length(b2_samples)
b1 <- b12[! b12$SampleID %in% b2_samples,]

# bridge normalized NPX for batch 1 are stored in NPX_normalized. To concatenate correct values, replace the raw NPX with bridged NPX
b1$NPX <- b1$NPX_normalized

length(unique(b1$SampleID))
# 217 samples
length(unique(b1$Assay))
# 689 proteins

shared_colnames = intersect(colnames(npx_br_clean_flt), colnames(b1))

b123 <- rbind(npx_br_clean_flt[,shared_colnames], b1[,shared_colnames])
b123 <- tibble(b123)
b123$UniProt <- gsub("NTproBNP","NT-proBNP", b123$UniProt)


b123_log <- check_npx(df = b123)

b123 <- b123 %>%
  group_by(SampleID) %>%
  mutate(sample_median_NPX = median(NPX, na.rm = TRUE)) %>%
  ungroup()

fwrite(b123, file =  "../../all_batches/data/olink_all_batches.no_flt.txt", quote = F, sep = "\t", row.names = F)

b123_pca <- olink_pca_plot(b123,color_g='geo', check_log = b123_log)
b123_pca <- olink_pca_plot(b123,color_g='Batch', check_log = b123_log)

pca_dat <- b123_pca[[1]]@data
pca_dat <- left_join(pca_dat, unique(b123[,c("SampleID", "sample_median_NPX")]),by = "SampleID")
ggplot(pca_dat, aes(x = PCX, y = PCY, color = sample_median_NPX)) + geom_point() + theme_minimal() + scale_color_viridis_c()

olink_qc_plot(
  df = b123,
  color_g = "Batch",
  label_outliers = T, 
  check_log = b123_log, IQR_outlierDef = 4
)

cat("All batches: number of samples: ", length(unique(b123$SampleID)), "for", length(unique(gsub("_.*","",b123[b123$SampleType == "SAMPLE",]$SampleID))), "unique women\n")
#Batch 3: number of samples:  197 for 75 unique women
cat("Batch 3: number of proteins:", length(unique(b123$Assay)), "\n")
#Batch 3: number of proteins: 5407 

sample_info <- read.delim("../../all_batches/data/olink_all_batches_batchinfo.txt", sep = "\t", as.is = T, check.names = F)
sample_info = sample_info[sample_info$has_olink_data == T & sample_info$pass_olink_QC == T,]
sample_info$ID = gsub("_.*","", sample_info$SampleID)

addmargins(table(sample_info[,c("batch", "city")]))

b123_wide <- as.data.frame(b123[,c("SampleID", "Assay", "NPX")] %>%
                             pivot_wider(names_from = "Assay", values_from = "NPX"))
dim(b123_wide)
b123_wide[grepl("^[0-9]",b123_wide$SampleID), "SampleID"] <- paste0("T", b123_wide[grepl("^[0-9]",b123_wide$SampleID), "SampleID"])

fwrite(b123_wide, file =  "../../all_batches/data/olink_all_batches.no_flt.wide.txt", quote = F, sep = "\t", row.names = F)


# Summary of missing rate
sample_na <- data.frame(
  SampleID = b123_wide$SampleID,
  NA_n = rowSums(is.na(b123_wide[, -1])),
  NA_pct = rowMeans(is.na(b123_wide[, -1])) * 100
)
sample_na <- left_join(sample_na, sample_info[,c("SampleID", "batch")], by = "SampleID")
View(sample_na)
ggplot(sample_na, aes(x = SampleID, y = NA_pct, color = batch)) + geom_point() + theme_minimal()


protein_na <- data.frame(
  Protein = names(b123_wide)[-1],
  NA_n = colSums(is.na(b123_wide[, -1])),
  nonNA_n = colSums(!is.na(b123_wide[, -1])),
  NA_pct = colMeans(is.na(b123_wide[, -1])) * 100
)

View(protein_na)
ggplot(protein_na, aes(x = Protein, y = NA_n)) + geom_point() + theme_minimal() + ylab("number of samples missing")


###
### Remove proteins with many below LOD values
###

num_below_lod <- b123[,c("Assay",  "NPX" ,"LOD", "Count")] %>%
  group_by(Assay) %>%
  mutate(above_lod = NPX >= LOD) %>%
  summarise(
    count_above_lod = sum(above_lod == TRUE, na.rm = T),
    count_below_lod = sum(above_lod == FALSE, na.rm = T),
    mean_abund_count = mean(Count, na.rm = T),
    fraction_above_lod = mean(NPX >= LOD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(count_below_lod)

ggplot(num_below_lod, aes(x = reorder(Assay, count_above_lod), y = count_above_lod)) +
  geom_bar(stat = "identity", fill = 'dodgerblue3') +
  geom_hline(yintercept = 50, color = 'red') + 
  geom_hline(yintercept = 150, color = 'red') + 
  labs(title = "Number of samples above LOD",
       x = "Assay",
       y = "# samples above LOD") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))


prots_to_remove_lod <- na.omit(num_below_lod[num_below_lod$count_above_lod < 150, ]$Assay)
cat("The number of proteins to remove, because they have less than 150 samples above LOD:\n")
length(prots_to_remove_lod)
cat(length(na.omit(num_below_lod[num_below_lod$count_above_lod >= 150, ]$Assay)), "proteins remaining\n")


b123_wide_lod_flt <- b123_wide[,! colnames(b123_wide) %in% prots_to_remove_lod]
dim(b123_wide_lod_flt)
fwrite(b123_wide_lod_flt, file =  "../../all_batches/data/olink_all_batches.lod150_wide.txt", quote = F, sep = "\t", row.names = F)

#prot_summary$pass_LOD_filter <- ifelse(! prot_summary$Assay %in% prots_to_remove_lod, T, F)


###
#### Remove outliers
###

remove_outliers_per_feature <- function(d, sd_cutoff = 4) {
  zscore <- scale(d)  # Compute z-scores
  d[abs(zscore) > sd_cutoff] <- NA  # Replace outliers with NA
  return(d)
}


remove_outliers_dataframe <- function(df, sd_cutoff = 4) {
  # Create a copy of the data frame to store the outlier mask
  outlier_mask <- df %>%
    mutate(across(.cols = -c(SampleID), .fns = ~abs(scale(.)) > sd_cutoff))
  
  # Apply the outlier removal
  df_cleaned <- df %>%
    mutate(across(.cols = -c(SampleID), .fns = ~remove_outliers_per_feature(., sd_cutoff)))
  
  # Return the cleaned data and the outlier mask
  list(cleaned_data = df_cleaned, outlier_mask = outlier_mask)
}



rm_outliers_tmp <- remove_outliers_dataframe(b123_wide_lod_flt,  sd_cutoff = 4)
b123_rm_outliers <- rm_outliers_tmp$cleaned_data


write.table(b123_rm_outliers, file = "../../all_batches/data/olink_all_batches.lod150_wide_rm_outliers_4sd.txt", quote = F, sep = "\t", row.names = FALSE)

num_outliers <- data.frame(colSums(rm_outliers_tmp$outlier_mask[,2:ncol(rm_outliers_tmp$outlier_mask)], na.rm = T)) %>%
  rownames_to_column(var = "prot")
colnames(num_outliers)[2] <- "num_outliers_4sd"

rm_outliers_tmp_6 <- remove_outliers_dataframe(b123_wide_lod_flt,  sd_cutoff = 6)
num_outliers_6 <- data.frame(colSums(rm_outliers_tmp_6$outlier_mask[,2:ncol(rm_outliers_tmp_6$outlier_mask)], na.rm = T)) %>%
  rownames_to_column(var = "prot")
colnames(num_outliers_6)[2] <- "num_outliers_6sd"

num_outliers <- left_join(num_outliers, num_outliers_6, by = 'prot')
num_outliers <- num_outliers[order(num_outliers$num_outliers_4sd, decreasing = T),]



ggplot(num_outliers[num_outliers$num_outliers_4sd > 3,], aes(x = reorder(prot, -num_outliers_4sd), y = num_outliers_4sd)) +
  geom_col(aes(fill = "4SD Outliers"), width = 0.7) +
  geom_col(aes(y = num_outliers_6sd, fill = "6SD Outliers"), width = 0.7) +
  scale_fill_manual(values = c("4SD Outliers" = "lightblue", "6SD Outliers" = "darkblue")) + 
  labs(title = "Outliers by Protein",
       x = "Protein",
       y = "Number of Outliers",
       fill = "Outlier Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3))


write.table(num_outliers, file = "../../all_batches/data/num_outliers_4sd.txt", quote = F, sep = "\t")




