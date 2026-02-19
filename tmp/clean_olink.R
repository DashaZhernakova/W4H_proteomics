library(dplyr)
library(ggplot2)

library(tidyr)
library(pheatmap)
library(patchwork)
library(tibble)
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")

setwd("/Users/Dasha/work/Sardinia/W4H/olink/data/")


################################################################################
# 1. Remove technical samples and assays from the data table
################################################################################


inf <- read.delim("INF_Q-14695_NPX_2024-09-28.txt", as.is = T, check.names = F, sep = "\t")
cvd <- read.delim("CVD_Q-08150_NPX_2024-09-28.txt", as.is = T, check.names = F, sep = "\t")

filter_table <- function(d){
  d <- d[,-1]
  d <- d[d$AssayType != "ext_ctrl",]
  d_clean <- d[d$SampleType == "SAMPLE" & d$AssayType == "assay",]
  d_clean$Timepoint <- gsub(".*_", "", d_clean$SampleID)
  d_clean$ID <- gsub("_.*", "", d_clean$SampleID)
  d_clean <- d_clean[d_clean$Normalization != "EXCLUDED",]
  cat("Total number of samples:", length(unique(d_clean$SampleID)), "\n\tNumber of individuals:", length(unique(d_clean$ID)), "\n\tNumber of proteins: ", length(unique(d_clean$Assay)), "\n")
  return(d_clean)
}

inf2 <- filter_table(inf)
cvd2 <- filter_table(cvd)

################################################################################
# 2. Remove 100_2 - the sample that failed CVD completely and partially INF
################################################################################

inf2 <- inf2[inf2$SampleID != '100_2',]
cvd2 <- cvd2[cvd2$SampleID != '100_2',]


################################################################################
# 3. Combine INF and CVD panels
################################################################################

merged_long <- rbind(inf2[,c("SampleID", "Assay", "NPX", "Panel")], cvd2[,c("SampleID", "Assay", "NPX", "Panel")])
merged_long$Assay2 <- paste0(merged_long$Assay, "_", merged_long$Panel)
merged_wide <-merged_long[,c("SampleID", "Assay2", "NPX")] %>%
  pivot_wider(names_from = Assay2, values_from = NPX)


# Combine clean without panel name
# For the overlapping proteins keep the INF version as default and add a suffix to the CVD version:
overlap <- intersect(unique(inf2$Assay), unique(cvd$Assay))
merged_long2 <- merged_long
merged_long2[merged_long2$Assay2 %in% c(paste0(overlap, "_Cardiometabolic")), "Assay"] <- merged_long2[merged_long2$Assay2 %in% c(paste0(overlap, "_Cardiometabolic")), "Assay2"]
merged_long2$Assay2 <- NULL
merged_long2$Panel <- NULL
merged_long2$Timepoint <- gsub(".*_","",merged_long2$SampleID)
merged_long2$ID <- gsub("_.*","",merged_long2$SampleID)
merged_long2$Assay <- gsub("_Cardiometabolic", "_CVD", merged_long2$Assay)
merged_wide2 <-merged_long2[,c("SampleID", "Assay", "NPX")] %>%
  pivot_wider(names_from = Assay, values_from = NPX)

num_cols_with_na <- sum(colSums(is.na(merged_wide2)) > 0)
num_rows_with_na <- sum(rowSums(is.na(merged_wide2)) > 0)

num_rows_with_na
num_cols_with_na

################################################################################
# 4. Write combined cleaned tables
################################################################################

ID <- gsub("_.*", "", merged_wide2$SampleID)
TP <- gsub(".*_", "", merged_wide2$SampleID)

merged_wide2 <- cbind(ID, TP, merged_wide2)

write.table(merged_wide2, file = "olink_clean_CVD+INF.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(merged_long2, file = "olink_clean_CVD+INF_long.txt", quote = F, sep = "\t", row.names = FALSE)


################################################################################
# 5. Remove below LOD values
################################################################################
merged_data <- read.delim("olink_clean_CVD+INF.txt", check.names =  F, sep = "\t", as.is = T, colClasses = c(ID = "character"))

lod_values <- read.delim("Explore 3072_Fixed LOD_2024-12-19.csv", sep = ";",check.names =  F, as.is = T)
lod_values <- lod_values[lod_values$DataAnalysisRefID %in% c("E70006", "E50007"),]
lod_values[lod_values$Assay == "TNF" & lod_values$Panel == "Cardiometabolic", "Assay"] <- "TNF_CVD"
lod_values[lod_values$Assay == "IL6" & lod_values$Panel == "Cardiometabolic", "Assay"] <- "IL6_CVD"
lod_values[lod_values$Assay == "CXCL8" & lod_values$Panel == "Cardiometabolic", "Assay"] <- "CXCL8_CVD"

lod_map <- setNames(lod_values$LODNPX, lod_values$Assay)

all_prots <- colnames(merged_data)[4:ncol(merged_data)]

# Apply the function to each protein column
cleaned_data_rm_below_lod <- merged_data
num_below_lod <- data.frame(protein = all_prots, num_samples = 0, num_below_lod = 0)
row.names(num_below_lod) <- all_prots

for (protein in all_prots) {
  lod_value <- lod_map[protein]
  column <- cleaned_data_rm_below_lod[[protein]]
  num_below_lod[protein, "num_samples"] <- sum(!is.na(column))
  num_below_lod[protein, "num_below_lod"] <- length(column[!is.na(column) & column < lod_value])
  num_below_lod[protein, "num_valid_samples"] <- num_below_lod[protein, "num_samples"] - num_below_lod[protein, "num_below_lod"]
  num_below_lod[protein, "num_valid_indiv"] <- length(unique(cleaned_data_rm_below_lod[!is.na(column) & column > lod_value,"ID"]))
  column[column < lod_value] <- NA  
  cleaned_data_rm_below_lod[[protein]] <- column
}
num_below_lod$perc_below_lod <- num_below_lod$num_below_lod / num_below_lod$num_samples
num_below_lod = left_join(num_below_lod, lod_values[,c("Assay", "LODNPX")], by = c("protein" = "Assay"))

pdf("num_below_lod.pdf",height = 7, width = 20)
ggplot(num_below_lod, aes(x = fct_reorder(protein, -num_below_lod), y = num_below_lod)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 263, color = 'grey') + 
  theme(axis.text.x = element_text(hjust = 1, size = 4, angle = 45),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlab("protein")
dev.off()

write.table(cleaned_data_rm_below_lod, file = "olink_clean_CVD+INF_rm_below_lod_keep_NA.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(num_below_lod, file = "num_below_lod.txt", quote = F, sep = "\t", row.names = FALSE)


# remove proteins that have more than half of the samples below LOD
prots_half_below_lod <- num_below_lod[num_below_lod$perc_below_lod > 0.5,"protein"]
length(prots_half_below_lod)
prots_less_30_indiv <- num_below_lod[num_below_lod$num_valid_indiv < 30, "protein"]
prots_less_100_samples <- num_below_lod[num_below_lod$num_valid_samples < 100, "protein"]

cleaned_data_rm_below_lod_flt <- cleaned_data_rm_below_lod

#cleaned_data_rm_below_lod_flt[,prots_half_below_lod] <- NULL
#write.table(cleaned_data_rm_below_lod_flt, file = "olink_clean_CVD+INF_rm_below_lod_prot_with_more_half_samples.txt", quote = F, sep = "\t", row.names = FALSE)

#cleaned_data_rm_below_lod_flt[,prots_less_30_indiv] <- NULL
#write.table(cleaned_data_rm_below_lod_flt, file = "olink_clean_CVD+INF_rm_below_lod_more_30_indiv.txt", quote = F, sep = "\t", row.names = FALSE)

cleaned_data_rm_below_lod_flt[,prots_less_100_samples] <- NULL
write.table(cleaned_data_rm_below_lod_flt, file = "olink_clean_CVD+INF_rm_below_lod_more_100_samples.txt", quote = F, sep = "\t", row.names = FALSE)


################################################################################
# 6. Set extreme outliers > 4 SDs from the mean to NA
################################################################################

remove_outliers("olink_clean_CVD+INF_rm_below_lod_more_100_samples.txt")
remove_outliers("olink_clean_CVD+INF_rm_below_lod_more_30_indiv.txt")
remove_outliers("olink_clean_CVD+INF_rm_below_lod_prot_with_more_half_samples.txt")


remove_outliers_per_feature <- function(d, sd_cutoff = 4, iqr_cutoff = 3, method = 'zscore') {
  if (method == 'zscore'){
    zscore <- scale(d)  # Compute z-scores
    d[abs(zscore) > sd_cutoff] <- NA  # Replace outliers with NA
  } else if (method == 'IQR'){
    q <- quantile(d, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr <- diff(q)
    d[d < (q[1] - iqr_cutoff * iqr) | d > (q[2] + iqr_cutoff * iqr)] <- NA
  } else {
    stop("Wrong method, should be zscore or IQR")
  }
  return(d)
}

remove_outliers_dataframe <- function(df, sd_cutoff = 5) {
  # Create a copy of the data frame to store the outlier mask
  outlier_summary <- df %>%
     select(-c(SampleID, ID, TP)) %>%  # Exclude specific columns
     summarise(across(everything(), ~any(abs(scale(.)) > sd_cutoff, na.rm = TRUE))) %>%
     pivot_longer(cols = everything(), names_to = "column", values_to = "has_outliers")
  
  # Apply the outlier removal
  df_cleaned <- df %>%
    mutate(across(.cols = -c(SampleID, ID, TP), .fns = ~remove_outliers_per_feature(., sd_cutoff)))
  
  # Return the cleaned data and the outlier mask
  list(cleaned_data = df_cleaned, outlier_mask = outlier_summary)
}

plot_features_with_outliers <- function(df, features_with_outliers, cutoff = 4, output_dir = NULL, method = 'zscore') {
  
  # Iterate over each feature and create a plot
  plots <- list()
  for (feature in features_with_outliers) {
    feature_values <- df[,feature]
    if (method == 'zscore'){
      zscore <- scale(feature_values)
      threshold_upper <- mean(feature_values, na.rm = TRUE) + cutoff * sd(feature_values, na.rm = TRUE)
      threshold_lower <- mean(feature_values, na.rm = TRUE) - cutoff * sd(feature_values, na.rm = TRUE)
    } else {
      q <- quantile(feature_values, probs = c(0.25, 0.75), na.rm = TRUE)
      iqr <- diff(q)
      
      threshold_lower <- q[1] - cutoff * iqr
      threshold_upper <- q[2] + cutoff * iqr
    }
    # Create the plot
    p <- ggplot(df, aes(x = SampleID, y = .data[[feature]], color = ID)) +
      geom_point() +
      geom_hline(yintercept = c(threshold_lower, threshold_upper), color = "red", linetype = "dashed") +
      labs(
        title = paste("Outliers in", feature),
        x = "SampleID",
        y = feature
      ) +
      theme_minimal() +
      theme(legend.position="none")
    
    # Save the plot if output_dir is specified
    if (!is.null(output_dir)) {
      ggsave(filename = file.path(output_dir, paste0(feature, "_outliers.png")), plot = p)
    }
    
    # Store the plot in a list
    plots[[feature]] <- p
  }
  
  return(plots)
}

remove_outliers <- function(fname,  sd_cutoff = 4, make_plots = F) {
  cleaned_data_rm_below_lod <- read.delim(fname, check.names =  F, sep = "\t", as.is = T, colClasses = c(ID = "character"))
  
  res <- remove_outliers_dataframe(cleaned_data_rm_below_lod,  sd_cutoff)
  cleaned_data <- res$cleaned_data
  outliers <- res$outlier_mask[res$outlier_mask$has_outliers == T, ]$column
  length(outliers)
  write.table(cleaned_data, file = paste0(gsub(".txt$","",fname), "_rm_outliers_", sd_cutoff, "sd.txt"), quote = F, sep = "\t", row.names = FALSE)

  if (make_plots) {
    plots <- plot_features_with_outliers(cleaned_data_rm_below_lod, outliers, cutoff = 4)
    return(plots)
  }
}