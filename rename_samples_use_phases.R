

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12/")

fname = "data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.txt"
fname = "../../phenotypes/batch12/cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.txt"
fname = "results12/covariates_olink_batch12.txt"
fname = "mb_covariates.txt"

phases <- read.delim("../../phenotypes/batch12/phases_251125.csv", as.is = T, check.names = F, sep = ",")

phases <- unique(phases) %>%
  mutate(
    phase_name = case_when(
      Phase == 1 ~ "F",
      Phase == 2 ~ "O",
      Phase == 3 ~ "EL",
      Phase == 4 ~ "LL",
      .default = "other"
    ),
    ID = gsub("_.*","", Code),
    SampleID_new = paste0(ID, "_", phase_name)
  )


d <- read.delim(fname, as.is = T, check.names = F, sep = "\t")
if ("Code" %in% colnames(d)) colnames(d) <- gsub("Code", "SampleID", colnames(d))
d[grepl("^[0-9]",d$SampleID), "SampleID"] <- paste0("X", d[grepl("^[0-9]",d$SampleID), "SampleID"])


nrow(d)
nrow(phases)
length(intersect(d$SampleID, phases$Code))
length(unique(gsub("_.*","",(d$SampleID))))

d_phases <- inner_join(phases[,c("Code", "SampleID_new")], d, by = c("Code" = "SampleID"))
d$SampleID[!d$SampleID %in% d_phases$Code]
d_phases$Code <- NULL
d_phases$ID = NULL

if (!covariates_file) {
  d_avg <- d_phases %>%
    group_by(SampleID_new) %>%
    summarise(
      across(where(is.numeric), ~ mean(., na.rm = TRUE)),  # Average numeric columns
      .groups = 'drop'
    )
} else {
d_avg <- d_phases %>%
  group_by(SampleID_new) %>%
  summarise(
    across(where(is.numeric), ~ mean(., na.rm = TRUE)),  # Average numeric columns
    batch = ifelse(any(batch == "batch2"), "batch2", first(batch)),  # Prefer batch2 if exists
    from = first(from),  # Take first from value (or use similar logic as batch if needed)
    .groups = 'drop'
  )
}
colnames(d_avg)[1] <- "SampleID"
length(unique(gsub("_.*","",(d_avg$SampleID))))
colnames(d_avg) <- gsub("^BES$", "17BES", colnames(d_avg))
d_avg$TP <- NULL
write.table(d_avg, file = paste0(gsub(".txt$","", fname), ".phase_avg.txt"), quote = F, sep = "\t", row.names = FALSE)

