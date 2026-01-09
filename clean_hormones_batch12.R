my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
my_colors <- c("#b71f57", "#96d1aa", "#099197", "#112f2c")
library(dplyr)
library(lme4)
library(lubridate)
library(lmerTest)
library(patchwork)
library(gridExtra)

set.seed(123)
setwd("/Users/Dasha/work/Sardinia/W4H/phenotypes/batch12/")

d <- read.delim("ALLSAMPLES_hormone_levels_original_file_13102025.csv", 
                sep = ",", as.is = T, check.names = F)
colnames(d) <- gsub("Record ID", "ID", colnames(d))
colnames(d) <- gsub("Numero Visita", "TP", colnames(d))
colnames(d) <- gsub("data analisi profilo ormonale", "date_of_analysis", colnames(d))
colnames(d) <- gsub(" .*","",colnames(d))
d$ID <- gsub("_", "", d$ID)
d$ID <- gsub("ID", "X", d$ID)
d[d == "-"] <- NA
d <- d[!is.na(d$date_of_analysis),]
colnames(d) <- gsub("^17BES", "X17BES", colnames(d) )

d$batch <- 1
d[endsWith(d$date_of_analysis, "24"), "batch"] <- 2
d[endsWith(d$date_of_analysis, "25"), "batch"] <- 3

d$date_of_analysis <- as.Date(d$date_of_analysis, format = "%d/%m/%y")

dates <- read.delim("date_of_collection_clean.txt", sep = "\t", as.is = T, check.names = F)
dates[grepl("^[0-9]",dates$SampleID), "SampleID"] <- paste0("X", dates[grepl("^[0-9]",dates$SampleID), "SampleID"])

d <- left_join(d, dates, by = c("ID" = "SampleID", "TP" = "visit_number"))

d$storage <- interval(d$date_collection, d$date_of_analysis) / months(1)

d[is.na(d$storage),]

write.table(d[,c("ID", "TP", "storage")], file = "hormone_storage_months.txt", sep = "\t", quote = F, row.names = F)


d$Note <- NULL
d$date_of_analysis <- NULL
d$date_collection <- NULL

d$SampleID <- paste0(d$ID, "_",d$TP)


d <- d %>% dplyr::select(SampleID, ID, TP, storage, batch, CITY, everything())
  
d <- d %>%
  mutate(across(-c(SampleID, ID, TP, storage, batch, CITY), as.numeric))

pdf("storage_distribution.pdf", width = 5, height = 5)
ggplot(d, aes(x = storage, color = CITY, group = CITY)) + geom_density() + theme_minimal() + xlab("storage time in months")
dev.off()

#number of samples
table(d$CITY)
#number of women
table(unique(d[,c("ID", "CITY")])$CITY)
#number of samples per TP
table(d[,c("TP", "CITY")])

d_clean <- d
#
# Calculate HOMA
#
glucose <- read.delim("cleaned_phenotypes_251125_uniformed_adjusted.csv", as.is = T, check.names = F, sep = ",")
glucose <- na.omit(glucose[,c("Code", "GL")])

d <- left_join(d, glucose, by = c("SampleID" = "Code"))
d$HOMA_IR <- d$GL * d$INS / 405
d$HOMA_B <- 360 * d$INS / (d$GL - 63)
d$GL <- NULL

d[!is.na(d$HOMA_B) & d$HOMA_B < 0, ]$HOMA_IR <- NA
d[!is.na(d$HOMA_B) & d$HOMA_B < 0, ]$HOMA_B <- NA

all_hormones <- c("INS", "PROG", "TSH", "FT4", "LH", "FSH", "X17BES", "TST", "PRL", "HOMA_B", "HOMA_IR")

#
# Distributions
#

plot_list = list()
for (ph in all_hormones){
  plot_list[[ph]] <- ggplot(d, aes(x = !!ensym(ph))) + geom_density() + theme_bw()
  plot_list[[paste0(ph, "_log")]] <- ggplot(d, aes(x = log( !!ensym(ph) ) )) + geom_density() + theme_bw()
}
pdf("plots/hormone_distributions_batch12.pdf", height = 20, width = 20)
grid.arrange(grobs = plot_list, ncol = 5, nrow = 5)  
dev.off()




#
# Log-transform all hormones
#
d_log <- d
d_log[, all_hormones] <- log(d[, all_hormones])



#
# Differences between cities
#

d_log4pca <- d_log[,7:ncol(d_log)]
row.names(d_log4pca) <- d_log$SampleID
d_log4pca[,c("SampleID", "TSH", "FT4", "TST")] <- NULL
pca <- stats::prcomp(d_log4pca, scale = TRUE, 
                     center = TRUE)

x <- pca$x[,1:5]
x <- as.data.frame(x)
x$SampleID <- row.names(x)
x <- left_join(x, d_log[,c("SampleID", "CITY", "storage", "batch")], by = "SampleID")

pdf("plots/hormone_differences_pca.pdf", width = 10, height = 15)
p1 <- ggplot(x, aes(PC1, PC2, color = CITY)) + geom_point(alpha = 0.4) + theme_minimal() 
p2 <- ggplot(x, aes(PC3, PC4, color = CITY)) + geom_point(alpha = 0.4) + theme_minimal() 

p3 <- ggplot(x, aes(PC1, PC2, color = batch)) + geom_point(alpha = 0.8) + theme_minimal() 
p4 <- ggplot(x, aes(PC3, PC4, color = batch)) + geom_point(alpha = 0.8) + theme_minimal() 

print((p1 + p2) / (p3 + p4)/ plot_spacer())

plot_list = list()
for (ph in colnames(d)[7: ncol(d_log)]){
  plot_list[[ph]] <- ggplot(d_log, aes(x = CITY, y = !!ensym(ph))) + geom_violin(trim = F) + geom_boxplot(width = 0.2) + theme_minimal() 
}
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
dev.off()


#
# Batch effects
#
d$batch <- relevel(as.factor(d$batch), ref = "3")
d_log$batch <- as.factor(d_log$batch)

batch_effects <- data.frame(matrix(nrow = ncol(d) - 6, ncol = 3))
colnames(batch_effects) <- c("pheno", "batch3_1_pval", "batch3_2_pval")
cnt <- 1
for (ph in all_hormones){
  cat(ph, "\n")
  
  subs <- na.omit(d[,c("ID", "TP", ph, "batch")])
  #subs$batch <- as.factor(subs$batch)
  colnames(subs)[3] <- 'pheno'
  if (! ph %in% c('TSH', 'FT4', 'TST')) {
    lmm_fit <- lmer(pheno ~ TP + batch + (1|ID), data = subs)
    pval1 <- summary(lmm_fit)$coefficients['batch1', 5]
    pval2 <- summary(lmm_fit)$coefficients['batch2', 5]
  } else {
    lm_fit <- lm(pheno ~ batch , data = subs)
    pval1 <- summary(lm_fit)$coefficients['batch1', 4]
    pval2 <- summary(lm_fit)$coefficients['batch2', 4]
  }
  batch_effects[cnt,] <- c(ph, pval1, pval2)
  cnt <- cnt + 1
}
batch_effects$batch3_1_pval <- as.numeric(batch_effects$batch3_1_pval)
batch_effects$batch3_2_pval <- as.numeric(batch_effects$batch3_2_pval)
View(na.omit(batch_effects))
write.table(batch_effects, file = "hormone_batch_effects_2.txt", sep = "\t", quote = F)


plot_list <- list()
for (ph in all_hormones){
  ph_name = sym(ph)
  plot_list[[ph]] <- ggplot(d_log, aes(x = batch, y = !!ph_name, group = batch)) + facet_wrap(~TP) + geom_violin(trim = F) + geom_boxplot(width = 0.3) + theme_bw()
}
pdf("plots/hormones_log_batch_effect_boxplots.pdf", height = 15, width = 15)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 3)  
dev.off()

plot_list <- list()
for (ph in all_hormones){
  ph_name = sym(ph)
  plot_list[[ph]] <- ggplot(d_log, aes( x = !!ph_name, group = batch, fill = batch)) +  geom_density(alpha = 0.5) + theme_bw() + facet_wrap(~TP) + xlab(paste0("log ", ph))
}
pdf("plots/hormones_batch_effect_distributions.pdf", height = 15, width = 20)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 3)  
dev.off()




#
# Correct for batch effect by matching the means of the affected hormones and regress out storage time
#

# get the effect of storage time

storage_pvals <- data.frame()
for (ph in all_hormones){
  if (! ph %in% c('TSH', 'FT4', 'TST')) {
    fo <- as.formula(paste(ph, "~ storage + (1|ID)"))
    lmm_fit <- lmer(fo, data = d_log)
    pval <- summary(lmm_fit)$coefficients["storage", 5]
  } else {
    fo <- as.formula(paste(ph, "~ storage"))
    lm_fit <- lm(fo, data = d_log)
    pval <- summary(lm_fit)$coefficients["storage", 4]
  }
  storage_pvals <- rbind(storage_pvals, c(ph, pval))
}
colnames(storage_pvals) <- c("hormone", "storage_effect_pval")
storage_pvals$storage_effect_pval <- as.numeric(storage_pvals$storage_effect_pval)
View(storage_pvals)


# regress storage time
d_log_adj <- regress_covariates_lmm(data = d_log[, c("SampleID", "ID", "TP", all_hormones)], covar_data = d_log[,c("SampleID", "storage"),drop =F])
d_log_adj$ID <- NULL
d_log_adj$TP <- NULL
d_log_adj <- full_join(d_log[,! colnames(d_log) %in% all_hormones], d_log_adj, by = "SampleID")

# correct for batch effect
pheno_to_adj_batch <- all_hormones
d_matched_means <- d_log_adj[,! colnames(d_log_adj) %in% pheno_to_adj_batch]

for (ph in pheno_to_adj_batch){
  d_matched_means <- full_join(d_matched_means, make_equal_means(d_log_adj, ph), by = "SampleID")
}
d_matched_means$batch <- NULL
d_matched_means$storage <- NULL
d_matched_means$CITY <- NULL
write.table(d_matched_means, file = "blood_hormones_18102025_log_adj_storage_batch_withHOMA.txt", quote = F, sep = "\t", row.names = FALSE)




### Only correct for the batch effect
d_matched_means_2 <- d_log[,! colnames(d_log) %in% pheno_to_adj_batch]

for (ph in pheno_to_adj_batch){
  d_matched_means_2 <- full_join(d_matched_means_2, make_equal_means(d_log, ph), by = "SampleID")
}
d_matched_means_2$batch <- NULL
d_matched_means_2$storage <- NULL
d_matched_means_2$CITY <- NULL
write.table(d_matched_means_2, file = "blood_hormones_18102025_log_adj_batch_withHOMA.txt", quote = F, sep = "\t", row.names = FALSE)



# Check if there is batch effect left
d_matched_means$batch <- as.factor(d_matched_means$batch)

batch_effects_after <- data.frame(matrix(nrow = ncol(d_matched_means) - 6, ncol = 3))
colnames(batch_effects_after) <- c("pheno", "batch1_2_pval", "batch1_3_pval")
cnt <- 1
for (ph in all_hormones){
  cat(ph, "\n")
  
  subs <- na.omit(d_matched_means[,c("ID", "TP", ph, "batch")])
  subs$batch <- as.factor(subs$batch)
  colnames(subs)[3] <- 'pheno'
  if (! ph %in% c('TSH', 'FT4', 'TST')) {
    lmm_fit <- lmer(pheno ~ TP + batch + (1|ID), data = subs)
    pval2 <- summary(lmm_fit)$coefficients['batch2', 5]
    pval3 <- summary(lmm_fit)$coefficients['batch3', 5]
  } else {
    lm_fit <- lm(pheno ~ batch , data = subs)
    pval2 <- summary(lm_fit)$coefficients['batch2', 4]
    pval3 <- summary(lm_fit)$coefficients['batch3', 4]
  }
  batch_effects_after[cnt,] <- c(ph, pval2, pval3)
  cnt <- cnt + 1
}
batch_effects_after$batch1_2_pval <- as.numeric(batch_effects_after$batch1_2_pval)
batch_effects_after$batch1_3_pval <- as.numeric(batch_effects_after$batch1_3_pval)
View(na.omit(batch_effects_after))

plot_list = list()
for (ph in all_hormones){
 ph_name = sym(ph)
 plot_list[[ph]] <- ggplot(d_matched_means, aes(x = batch, y = !!ph_name, group = batch)) + facet_wrap(~TP) + geom_violin(trim = F) + geom_boxplot(width = 0.3) + theme_bw()
}
pdf("plots/hormones_log_batch_effect_boxplots_after_correction.pdf", height = 15, width = 15)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 3)  
dev.off()

plot_list = list()
for (ph in all_hormones){
  plot_list[[ph]] <- ggplot(d_matched_means, aes(x = CITY, y = !!ensym(ph))) + geom_violin(trim = F) + geom_boxplot(width = 0.2) + theme_minimal() 
}

pdf("plots/hormone_city_differences_after_batch_correction.pdf", width = 10, height = 15)
grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
dev.off()


d_matched_means$batch <- NULL
d_matched_means$storage <- NULL
d_matched_means$CITY <- NULL

d_log$batch <- NULL
d_log$storage <- NULL
d_log$CITY <- NULL


d_log_rm_outliers4 <- remove_outliers_dataframe(d_log,  sd_cutoff = 4)$cleaned_data
rm_outliers_tmp <- remove_outliers_dataframe(d_matched_means,  sd_cutoff = 4)
d_matched_means4 <- rm_outliers_tmp$cleaned_data


#
# Write final tables 
#
write.table(d, file = "blood_hormones_18102025.txt", quote = F, sep = "\t", row.names = FALSE)
write.table(d_log, file = "blood_hormones_18102025_log.txt", quote = F, sep = "\t", row.names = FALSE)




make_equal_means <- function(pheno_with_batch, ph){
  b1 <- pheno_with_batch[pheno_with_batch$batch == '1', c("SampleID",  "batch", ph)]
  b2 <- pheno_with_batch[pheno_with_batch$batch == '2', c("SampleID",  "batch", ph)]
  b3 <- pheno_with_batch[pheno_with_batch$batch == '3', c("SampleID",  "batch", ph)]
  
  ref_mean <- mean(b3[,ph], na.rm = T)
  m1 <- mean(b1[,ph], na.rm = T)
  m2 <- mean(b2[,ph], na.rm = T)
  mean_dif1 <- m1 - ref_mean
  mean_dif2 <- m2 - ref_mean
  
  #if (m1 > m2) stop("Mean of batch 1 is larger than the mean of batch2")
  
  b2[,ph] <- b2[,ph] - mean_dif2
  b1[,ph] <- b1[,ph] - mean_dif1
  
  adjusted <- rbind(b1, b2, b3)
  
  pheno_with_batch$batch <- as.factor(pheno_with_batch$batch)
  pheno_name = sym(ph)
  p1 <- ggplot(pheno_with_batch, aes (x = !!pheno_name, group = batch, fill = batch)) + 
    geom_density(alpha = 0.3) + 
    geom_vline(xintercept = m1, color = my_colors[1]) + 
    geom_vline(xintercept = m2, color = my_colors[2]) + 
    geom_vline(xintercept = ref_mean, color = my_colors[3]) + 
    scale_fill_manual(values = my_colors) +
    theme_minimal()
  
  p2 <- ggplot(adjusted, aes (x = !!pheno_name, group = batch, fill = batch)) + 
    geom_density(alpha = 0.3) + geom_vline(xintercept = ref_mean)  + 
    theme_minimal() 
  
  print(p1 + p2)
  adjusted$batch <- NULL
  
  return (adjusted)
}




regress_covariates_lmm <- function(data, covar_data, covars_longitudinal = T){
  
  if (!"SampleID" %in% colnames(covar_data) & covars_longitudinal) {
    covar_data <- cbind(paste0(covar_data$ID, "_",covar_data$TP), covar_data)
    colnames(covar_data)[1] <- "SampleID"
  }
  
  d_adj <- data[,c("SampleID", "ID", "TP")]
  
  data[,"TP"] <- NULL
  covar_data[,"TP"] <- NULL
  
  cnt <- 1
  for (ph in colnames(data)[3: (ncol(data))]){
    if (covars_longitudinal){
      covar_data$ID <- NULL
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    } else {
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "ID"))
    }
    colnames(subs)[3] <- 'pheno'
    
    if (! ph %in% c('TSH', 'FT4', 'TST')) {
      fo_lmm <- as.formula(paste("pheno ~ ", paste(colnames(covar_data)[-1], collapse = "+"), "+ (1|ID)"))
      lmm_fit <- lmer(fo_lmm, data = subs)
      subs[,ph] <- subs$pheno - predict(lmm_fit, re.form = NA)
    } else {
      fo_lm <- as.formula(paste("pheno ~ ", paste(colnames(covar_data)[-1], collapse = "+")))
      lm_fit <- lm(fo_lm, data = subs)
      subs[,ph] <- residuals(lm_fit)
    }
    
    d_adj <- left_join(d_adj, subs[, c("SampleID", ph)], by = "SampleID")
  }
  return(d_adj)
}


remove_outliers_per_feature <- function(d, sd_cutoff = 4) {
  zscore <- scale(d)  # Compute z-scores
  d[abs(zscore) > sd_cutoff] <- NA  # Replace outliers with NA
  return(d)
}


remove_outliers_dataframe <- function(df, sd_cutoff = 4) {
  # Create a copy of the data frame to store the outlier mask
  outlier_mask <- df %>%
    mutate(across(.cols = -c(SampleID, ID, TP), .fns = ~abs(scale(.)) > sd_cutoff))
  
  # Apply the outlier removal
  df_cleaned <- df %>%
    mutate(across(.cols = -c(SampleID, ID, TP), .fns = ~remove_outliers_per_feature(., sd_cutoff)))
  
  # Return the cleaned data and the outlier mask
  list(cleaned_data = df_cleaned, outlier_mask = outlier_mask)
}

