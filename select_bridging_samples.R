library(OlinkAnalyze)
library(stringr)
library(dplyr)
setwd('/Users/Dasha/work/Sardinia/W4H/olink//')

inf_data <- read_NPX("data/INF_Q-14695_NPX_2024-09-28.txt")
cvd_data <- read_NPX("data/CVD_Q-08150_NPX_2024-09-28.txt")

set.seed(123)

combined_data <- rbind(inf_data, cvd_data)
combined_data <- combined_data %>%
  mutate(IndividualID = gsub("_.*", "", SampleID))

# Select 100 bridge samples using Olink's olink_bridgeselector:
# Remove outliers, samples that fail QC, randomly select 150 samples that span the whole NPX range 
bridge_samples_all<- combined_data[combined_data$SampleType == 'SAMPLE',] %>% 
  olink_bridgeselector(sampleMissingFreq = 0.1,
                       n = 100)

# Randomly pick 1 timepoint per individual for the selected bridge samples
combined_data_unique <- combined_data[combined_data$SampleID %in% bridge_samples_all$SampleID,] %>%
  group_by(IndividualID) %>%
  slice_sample(n = 1) %>%          
  ungroup()

# randomly subset 50 samples
bridge_samples_unique <- sample(combined_data_unique$SampleID, 40)

pdf("data/bridge_samples_pca_unique_40.pdf")
combined_data[combined_data$SampleType == 'SAMPLE',] %>% 
  mutate(Bridge = ifelse(SampleID %in% bridge_samples_unique, "Bridge", "Sample")) %>% 
  olink_pca_plot(color_g = "Bridge")
dev.off()
write.table(bridge_samples_unique, file = "data/bridge_samples_unique_40.txt", quote = F, sep = "\t", row.names = FALSE)



