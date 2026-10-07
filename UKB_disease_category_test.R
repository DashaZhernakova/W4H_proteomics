library(dplyr)
library(readr)
library(stringr)
library(tidyr)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")

# read association results from both batch 1+2 and batch2
d1 <- read.delim("results12/intensity_shared_prots_261125/prot_vs_hormones.gam.spline.shared_prots.txt", as.is = T, check.names  = F, sep = "\t")
d2 <- read.delim("results12/intensity_batch2_prots_261125/prot_vs_hormones.gam.spline.shared_prots.txt", as.is = T, check.names  = F, sep = "\t")
d1$batch = "batch12"
d2$batch <- "batch2"
d <- rbind(d1[,c("prot", "pheno", "pval", "estimate","batch","BH_pval")], d2[,c("prot", "pheno", "pval", "estimate","batch","BH_pval")])
d$abs_estimate <- abs(d$estimate)


#### UKB disease associations
df_ukb <- read_csv("data/all_tested_2453_olink_annot_UKB.csv")
dis_cat <- read_csv("/Users/Dasha/work/resources/phecode_definitions1.2.csv")

df_long <- df_ukb[,c("Gene", "UKB Disease Risk")] %>%
  mutate(disease_string = str_split(`UKB Disease Risk`, ", ")) %>%  # split by comma+space
  unnest(disease_string) %>%
  mutate(
    # Remove HR, p-value etc. to get clean disease name
    disease_clean = str_remove(disease_string, "\\s+HR\\s*=\\s*[0-9\\.eE+-]+\\s*\\|\\s*p\\s*=\\s*[0-9\\.eE+-]+"),
    disease_clean = str_trim(disease_clean)
  )

dis_cat_clean <- dis_cat %>%
  mutate(phenotype_clean = str_trim(tolower(phenotype))) %>%
  select(phenotype_clean, category)

# add the few diseases missing in the phewas category DB
additional_dis_cat <- data.frame(matrix(c("all-cause mortality", "all-cause mortality",
                                          "cancer of bronchus lung", "neoplasms",
                                          "congestive heart failure nonhypertensive", "circulatory system",
                                          "esophagitis", "digestive",
                                          "gerd and related diseases","digestive",
                                          "nephritis nephrosis renal sclerosis", "genitourinary",
                                          "pulmonary collapse interstitial and compensatory emphysema", "respiratory",
                                          "recurrent seizures", "neurological"), byrow = T, ncol = 2))
colnames(additional_dis_cat) <- c("phenotype_clean", "category")

dis_cat_clean <- rbind(dis_cat_clean, additional_dis_cat)

df_with_cat <- df_long %>%
  mutate(disease_clean_lower = str_trim(tolower(disease_clean))) %>%
  left_join(dis_cat_clean, by = c("disease_clean_lower" = "phenotype_clean")) %>%
  select(Gene, category) %>%
  unique() 

df_with_cat_wide <- df_with_cat %>%
  group_by(Gene) %>%
  summarise(
    category_merged = toString(category),
    .groups = "drop"
  )

d_with_cat_wide <- df_with_cat %>%
  group_by(Gene) %>%
  summarise(
    category_merged = toString(category),
    .groups = "drop"
  )

# Are estimates higher for proteins associated with any specific UKB disease category?
d_annot <- left_join(d, df_with_cat_wide[, c("Gene", "category_merged")], by = c("prot" = "Gene"))
all_categories <- na.omit(unique(df_with_cat$category))

res_table_UKB <- data.frame(Disease_category = all_categories, KS_twosided = NA, KS_greater = NA, KS_less = NA)

for (p in all_categories){
  d_annot$category_bin <- grepl(pattern = p, x = d_annot$category_merged)
  cat(p, "\n")
  twosided <- ks.test(d_annot[d_annot$category_bin == T,"abs_estimate"], d_annot[d_annot$category_bin == F,"abs_estimate"])
  greater <- ks.test(d_annot[d_annot$category_bin == T,"abs_estimate"], d_annot[d_annot$category_bin == F,"abs_estimate"], alternative = "greater")
  less <- ks.test(d_annot[d_annot$category_bin == T,"abs_estimate"], d_annot[d_annot$category_bin == F,"abs_estimate"], alternative = "less")
  
  res_table_UKB[res_table_UKB$Disease_category == p,2:4] <- c(twosided$p.value, greater$p.value, less$p.value)
}
View(res_table_UKB)


res_table_UKB_per_hormone <- expand.grid(
  Disease_category = all_categories,
  hormone = unique(d_annot$pheno), 
  n_true = NA, n_false = NA,
  twosided = NA, greater = NA, less = NA, 
  stringsAsFactors = FALSE)


for (h in unique(d_annot$pheno)){
  d_h <- d_annot[d_annot$pheno == h,]
  for (p in all_categories){
    d_h$category_bin <- grepl(pattern = p, x = d_h$category_merged)
    cat(p, "\n")
    twosided <- wilcox.test(d_h[d_h$category_bin == T,"abs_estimate"], d_h[d_h$category_bin == F,"abs_estimate"])
    greater <- wilcox.test(d_h[d_h$category_bin == T,"abs_estimate"], d_h[d_h$category_bin == F,"abs_estimate"], alternative = "greater")
    less <- wilcox.test(d_h[d_h$category_bin == T,"abs_estimate"], d_h[d_h$category_bin == F,"abs_estimate"], alternative = "less")
    
    n_true <- length(d_h[d_h$category_bin == T, "prot"])
    n_false <- length(d_h[d_h$category_bin == F, "prot"])
    res_table_UKB_per_hormone[res_table_UKB_per_hormone$Disease_category == p & res_table_UKB_per_hormone$hormone == h,3:7]  <- c(n_true, n_false, twosided$p.value, greater$p.value, less$p.value)
  }
}
View(res_table_UKB_per_hormone)


p_cutoff <- max(d1[d1$BH_pval < 0.05, "pval"])
d_annot$signif <- d_annot$pval < p_cutoff

calculate_enrichment <- function(p, df, colname = "category_merged") {
  x <- data.frame(disease = grepl(pattern = p, x = df[,colname]), associated_hormones = df$signif)
  contingency_table <- table(
    x$disease == TRUE,         # In panel?
    x$associated_hormones == 1  # Associated with Hormone?
  )
  if (all(dim(contingency_table) == 2)) {
    
    test <- fisher.test(contingency_table)
    data.frame(
      gene_set = p,
      p_value = test$p.value,
      odds_ratio = test$estimate,
      count_in_set = sum(x$disease == T),
      overlap_count = contingency_table[2, 2] # Intersection of both
    )
  } else {
    data.frame(
      gene_set = p,
      p_value = NA,
      odds_ratio = NA,
      count_in_set = sum(x$pathway == T),
      overlap_count = NA # Intersection of both
    )
  }
}

calculate_enrichment_per_hormone <- function(p, df, h,colname = "category_merged") {
  df <- df[df$pheno == h,]
  x <- data.frame(disease = grepl(pattern = p, x = df[,colname]), associated_hormones = df$signif)
  contingency_table <- table(
    x$disease == TRUE,         # In panel?
    x$associated_hormones == 1  # Associated with Hormone?
  )
  if (all(dim(contingency_table) == 2)) {
    
    test <- fisher.test(contingency_table)
    data.frame(
      gene_set = p,
      hormone = h,
      p_value = test$p.value,
      odds_ratio = test$estimate,
      count_in_set = sum(x$disease == T),
      overlap_count = contingency_table[2, 2] # Intersection of both
    )
  } else {
    data.frame(
      gene_set = p,
      hormone = h,
      p_value = NA,
      odds_ratio = NA,
      count_in_set = sum(x$pathway == T),
      overlap_count = NA # Intersection of both
    )
  }
}


enrichment_results <- lapply(all_categories, calculate_enrichment, d_annot) %>%
  bind_rows() %>%
  mutate(
    fdr = ifelse(overlap_count > 0, p.adjust(p_value, method = "BH"), NA)
  )


dis_clean <- df_long %>%
  group_by(Gene) %>%
  summarise(
    disease_merged = toString(disease_clean),
    .groups = "drop"
  )
d_annot_dis <- left_join(d, dis_clean[, c("Gene", "disease_merged")], by = c("prot" = "Gene"))
all_diseases <- na.omit(unique(df_long$disease_clean))
d_annot_dis$signif <- d_annot_dis$pval < p_cutoff

enrichment_results_per_disease <- lapply(all_diseases, calculate_enrichment, d_annot_dis, "disease_merged") %>%
  bind_rows() %>%
  mutate(
    fdr = ifelse(overlap_count > 0, p.adjust(p_value, method = "BH"), NA)
  )



