library(OlinkAnalyze)
library(stringr)
library(dplyr)
setwd('/Users/Dasha/work/Sardinia/W4H/olink/batch12//')

b2 <- read_NPX("data/Q-08152_NPX_2025-09-17.parquet")
# Add LOD
b2 <- OlinkAnalyze::olink_lod(b2,
                              lod_method = "FixedLOD",
                              lod_file_path = "data/Explore HT_Fixed LOD.csv")

b2$city <- "T"
b2[startsWith(b2$SampleID, "B"),"city"] <- "B"
b2[startsWith(b2$SampleID, "S"),"city"] <- "S"

set.seed(150)
set.seed(111)

b2 <- b2 %>%
  mutate(IndividualID = gsub("_.*", "", SampleID))
b2 <- b2[b2$AssayQC != 'WARN',]

b2_high_missing_assays <- b2 %>%
  group_by(Assay) %>%
  summarise(missing_rate = 1 - sum(!is.na(NPX)) / num_samples ) %>% 
  filter(missing_rate > 0.1) %>% 
  arrange(desc(missing_rate))

cat("Proteins with missing rate > 10%:\n")
b2_high_missing_assays
if (nrow(b2_high_missing_assays) > 0) b2 <- b2[! b2$Assay %in% b2_high_missing_assays$Assay,]

# read the filtered dataset
x <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.txt", as.is = T, check.names = F, sep = "\t")
# select only batch 2 samples
x <- x[x$SampleID %in% b2$SampleID,]
dim(x)
missing_per_sample <- data.frame(SampleID = x$SampleID, num_missing = rowSums(is.na(x)))
table(missing_per_sample$num_missing)
samples_low_missing <- x$SampleID[rowSums(is.na(x)) < 4]
length(samples_low_missing)


### Check missing in filtered data
b2 <- left_join(b2, missing_per_sample, by = "SampleID")
b2$num_missing <- as.factor(b2$num_missing)
p1 <- b2[b2$SampleType == 'SAMPLE',] %>% 
  mutate(Bridge = ifelse(SampleID %in% bridge_samples_unique, "Bridge", "Sample")) %>% 
  olink_pca_plot(color_g = "num_missing")

dat <- p1[[1]]$data

dat$missing <- as.numeric(as.character(dat$colors))
sample_right <- dat[dat$missing < 4 & dat$PCX > 0.1 & dat$PCY > 0,]$SampleID

# Select 100 bridge samples using Olink's olink_bridgeselector:
# Remove outliers, samples that fail QC, randomly select 100 samples that span the whole NPX range 
bridge_samples_all<- b2[b2$SampleID %in% samples_low_missing,] %>% 
  olink_bridgeselector(sampleMissingFreq = 0.6,
                       n = 60)

# Randomly pick 1 timepoint per individual for the selected bridge samples
b2_unique <- b2[b2$SampleID %in% bridge_samples_all$SampleID,] %>%
  group_by(IndividualID) %>%
  slice_sample(n = 1) %>%          
  ungroup()

# randomly subset 32 samples
bridge_samples_unique <- unique(c(sample_right, sample(b2_unique$SampleID, 29)))
length(bridge_samples_unique)
pdf("data/bridge_samples_pca_unique_30.pdf")
b2[b2$SampleType == 'SAMPLE',] %>% 
  mutate(Bridge = ifelse(SampleID %in% bridge_samples_unique, "Bridge", "Sample")) %>% 
  olink_pca_plot(color_g = "Bridge")
dev.off()
write.table(bridge_samples_unique, file = "data/bridge_samples_unique_30.txt", quote = F, sep = "\t", row.names = FALSE)


### Check missing in filtered data
b2$num_missing <- as.factor(b2$num_missing)
b2[b2$SampleType == 'SAMPLE',] %>% 
     mutate(Bridge = ifelse(SampleID %in% bridge_samples_unique, "Bridge", "Sample")) %>% 
     olink_pca_plot(color_g = "num_missing")

b2[b2$SampleType == 'SAMPLE',] %>% 
  mutate(Bridge = ifelse(SampleID %in% bridge_samples_unique, "Bridge", "Sample")) %>% 
  olink_pca_plot(color_g = "city")
