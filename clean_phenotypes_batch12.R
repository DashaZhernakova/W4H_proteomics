library(dplyr)
library(ggplot2)
library(gridExtra)

fname = "../../phenotypes/batch12/cleaned_phenotypes_251125_uniformed_adjusted.csv"

d <- read.delim(fname, as.is = T, check.names = F, sep = ",", row.names =1)
d <-d[,!grepl("_AOU$", colnames(d))]
d$X.1 <- NULL
d$X <- NULL
d$cycle_index <- NULL

# remove duplicate rows that have NAs for phenotypes
#d_full <- d %>%
#  group_by(Code) %>%
#  arrange(desc(!is.na(GL))) %>% 
#  dplyr::slice(1) %>%  
#  ungroup()

homa<- read.delim("../../phenotypes/batch12/blood_hormones_18102025_log_adj_storage_batch_withHOMA.txt", as.is = T, check.names = F, sep = "\t")
d <- left_join(d, homa[c("SampleID", "HOMA_IR", "HOMA_B")], by = c("Code" = "SampleID"))

phenos <- c("GL", "AST", "ALT", "TRI", "COL", "HDL", "LDL")
#d <- d[,c("Code", "GL", "AST", "ALT", "TRI", "COL", "HDL", "LDL")]
#d[d == ""] <- NA

d$cycle_index <-NULL
d$ID <- NULL
d$X <- NULL
d$Visit_number <- NULL

plot_list <- list()
for (ph in phenos){
  ph_name = sym(ph)
  plot_list[[ph]] <- ggplot(d, aes(x = !!ph_name)) + geom_density() + theme_minimal()
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  

pheno_to_log <- c("AST", "ALT", "TRI")

d_log <- d
d_log[,pheno_to_log] <- log(d[,pheno_to_log])

plot_list <- list()
for (ph in phenos){
  ph_name = sym(ph)
  plot_list[[ph]] <- ggplot(d_log, aes(x = !!ph_name)) + geom_density() + theme_minimal()
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  

write.table(d_log, file = paste0(fname, "withHOMA.log_some.txt"), quote = F, sep = "\t", row.names = F)


# add phases
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


d_log_phases <- inner_join(phases[,c("Code", "SampleID_new")], d_log, by = "Code")
d_log$Code[!d_log$Code %in% d_log_phases$Code]

d_log_avg <- d_log_phases %>%
  group_by(SampleID_new) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)))

d_log_avg$Code = NULL
colnames(d_log_avg) <- gsub("SampleID_new", "SampleID", colnames(d_log_avg))

write.table(d_log_avg, file = paste0(fname, "withHOMA.log_some.phase_avg.txt"), quote = F, sep = "\t", row.names = F)
