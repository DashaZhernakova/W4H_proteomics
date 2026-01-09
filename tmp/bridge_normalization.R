library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12/data/")

b1 <-  list.files(path = "./",
                  pattern = "2024-09-28.parquet$",
                  full.names = TRUE) |>
  lapply(OlinkAnalyze::read_NPX)  |>
  dplyr::bind_rows()


# Remove samples that failed the batch 1 QC:
b1_failed_qc <- c("100_2", "018_4", "102_3", "102_2", "091_2", "099_4", "109_1")
b1 <- b1[! b1$SampleID %in% b1_failed_qc,]

# remove Assays that failed QC in batch 1
table(b1[,c("AssayQC", "Normalization")])
b1 <- b1[b1$Normalization != 'EXCLUDED',]

# Extract only real samples
table(b1[,c("SampleType","SampleQC")])
b1 <- b1[b1$SampleType == 'SAMPLE',]

b1$geo <- "Trieste"

cat("Batch 1: number of samples: ", length(unique(b1$SampleID)), "for", length(unique(gsub("_.*","",b1$SampleID))), "unique women\n")
cat("Batch 1: number of proteins:", nrow(unique(b1[b1$AssayType == 'assay', "Assay"])), "\n")

# leave only non-control assays
b1 <- b1[b1$AssayType == 'assay',]


b2 <- read_NPX("Q-08152_NPX_2025-09-17.parquet")

b2$geo <- "Trieste"
b2[startsWith(b2$SampleID, "S"), "geo"] <- "Sardinia"
b2[startsWith(b2$SampleID, "B"), "geo"] <- "Bologna"


table(b2[,c("AssayQC", "AssayType", "Normalization")])

cat("Number of proteins with QC status = WARN:", nrow(unique(b2[b2$AssayQC == 'WARN', "Assay"])), "\n")
unique(b2[b2$AssayQC == 'WARN', "Assay"])

cat("Removing these 9 proteins!\n")

b2 <- b2[b2$AssayQC != 'WARN',]

# Extract only real samples
table(b2[,c("SampleType","SampleQC")])
b2 <- b2[b2$SampleType == 'SAMPLE',]


## Remove assays with high missing rate
b2_wide <- b2[b2$AssayType == "assay",c("SampleID", "Assay", "NPX")] %>%
  pivot_wider(names_from = "Assay", values_from = "NPX")

b2_wide <- as.data.frame(b2_wide)
row.names(b2_wide) <- b2_wide$SampleID
b2_wide$SampleID <- NULL


# Identify proteins with missing rate > 10%
missing_rates <- colMeans(is.na(b2_wide)) 
high_missing_proteins <- names(missing_rates[missing_rates > 0.1])

cat("Removing proteins with high missing rate:\n")
print(missing_rates[high_missing_proteins])

b2 <- b2[! b2$Assay %in% high_missing_proteins,]
b2_wide <- b2_wide[,!colnames(b2_wide) %in% high_missing_proteins]

cat("Batch 1: number of samples: ", length(unique(b2$SampleID)), "for", length(unique(gsub("_.*","",b2$SampleID))), "unique women\n")
cat("Batch 1: number of proteins:", nrow(unique(b2[b2$AssayType == 'assay', "Assay"])), "\n")


# leave only non-control assays
b2 <- b2[b2$AssayType == 'assay',]


# check overlapping bridging sample ids and overlapping proteins
overlapping_samples <- intersect(b1$SampleID, b2$SampleID)
cat ("Number of overlapping samples usable for bridging:", length(overlapping_samples), "\n")

overlapping_prots <- intersect(unique(b1[b1$AssayType == 'assay',]$Assay), b2[b2$AssayType == "assay", ]$Assay)
cat ("Number of overlapping proteins:", length(overlapping_prots), "\n")


# Plot PCA to see if bridging samples are outliers
b1_before_br <- b1 |>
  mutate(Type = if_else(SampleID %in% overlapping_samples,
                        paste0("batch1_bridge"),
                        paste0("batch1_sample")),
         Batch = "batch1")

b2_before_br <- b2 |>
  dplyr::mutate(Type = if_else(SampleID %in% overlapping_samples,
                               paste0("batch2_bridge"),
                               paste0("batch2_sample")),
                Batch = "batch2")

### PCA plot
pca_b1 <- olink_pca_plot(df = b1_before_br,
                         color_g = "Type",
                         quiet = TRUE) 
pca_b2 <- olink_pca_plot(df = b2_before_br,
                         color_g = "Type",
                         quiet = TRUE)
pca_b2_geo <- olink_pca_plot(df = b2_before_br,
                             color_g = "geo",
                             quiet = TRUE)

pca_b1[[1]]

pca_b2[[1]]

pca_b2_geo[[1]]


# PCA myself

b2_pca <- stats::prcomp(b2_wide, scale = TRUE, 
              center = TRUE)
x <- b2_pca$x[,1:5]
x <- as.data.frame(x)
x$geo = "Trieste"
x[startsWith( row.names(x), "S"), "geo"] <- "Sardinia"
x[startsWith(row.names(x), "B"), "geo"] <- "Bologna"

ggplot(x, aes(PC1, PC2, color = geo)) + geom_point() + theme_minimal()

# Get the loadings (rotations) for PC1
pc1_loadings <- b2_pca$rotation[, "PC1"]
# Sort by absolute value to see most influential proteins
pc1_importance <- sort(abs(pc1_loadings), decreasing = TRUE)
write.table(pc1_importance, file = "PC1_loadings_abs.txt")




npx_df <- bind_rows(b1_before_br[b1_before_br$Assay %in% overlapping_prots,], 
                    b2_before_br[b2_before_br$Assay %in% overlapping_prots,])


# Plot NPX density before bridging normalization
npx_df %>%
  mutate(Panel = gsub("Olink ", "", Panel)) %>%
  ggplot(aes(x = NPX, fill = Batch)) +
  geom_density(alpha = 0.4) +
  olink_fill_discrete(coloroption = c("red", "darkblue")) +
  set_plot_theme() +
  ggtitle("Before bridging normalization: NPX distribution") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "top")


# get number of shared samples for each protein

npx_df_subs <- npx_df[,c("SampleID", "Assay", "NPX","Batch")]
npx_df_subs <- npx_df_subs[npx_df_subs$SampleID %in% overlapping_samples,]
res <- data.frame()
for (p in unique(npx_df_subs$Assay)){
  b1_samples <- npx_df_subs[npx_df_subs$Assay == p & npx_df_subs$Batch == "batch1", "SampleID"]
  b2_samples <- npx_df_subs[npx_df_subs$Assay == p & npx_df_subs$Batch == "batch2", "SampleID"]
  overlap <- length(intersect(b1_samples$SampleID, b2_samples$SampleID))
  res <- rbind(res, c(p, overlap))
}

colnames(res) <- c("protein", "num_overlapping_samples")
not_enough_samples <- res[res$num_overlapping_samples < 40,]

cat("Proteins with less than 40 bridging samples:\n")
not_enough_samples


## Perform normalization
overlap_samples_list <- list("DF1" = overlapping_samples,
                             "DF2" = overlapping_samples)

npx_br_data <- olink_normalization_bridge(project_1_df = b1_before_br,
                                          project_2_df = b2_before_br,
                                          bridge_samples = overlap_samples_list,
                                          project_1_name = "batch1",
                                          project_2_name = "batch2",
                                          project_ref_name = "batch2")




## check if including or removing technical assays matters


npx_br_data_no_tech <- olink_normalization_bridge(project_1_df = b1_before_br[b1_before_br$AssayType == 'assay',],
                                          project_2_df = b2_before_br[b2_before_br$AssayType == 'assay',],
                                          bridge_samples = overlap_samples_list,
                                          project_1_name = "batch1",
                                          project_2_name = "batch2",
                                          project_ref_name = "batch2")

