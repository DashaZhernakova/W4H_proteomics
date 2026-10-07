library(dplyr)

### Are there more associated proteins in CVD and INF panels than in other panels?

df2 <- read_csv("data/all_tested_2453_olink_annot.csv")
all_pathways <- unique(unlist(str_split(df2$`Human Pathways`, ", ")))
all_pathways <- all_pathways[!is.na(all_pathways)]
all_panels <- na.omit(unique(df2$Panel))
all_panels_category <- na.omit(unique(gsub("_.*","",df2$Panel)))


calculate_enrichment_human_pathways <- function(p, df2) {
  x <- data.frame(pathway = grepl(pattern = p, x = df2$`Human Pathways`), associated_hormones = df2$`associated_p1e-3`)
  contingency_table <- table(
    x$pathway == TRUE,         # In pathway?
    x$associated_hormones == 1  # Associated with Hormone?
  )
  if (all(dim(contingency_table) == 2)) {

  test <- fisher.test(contingency_table)
  data.frame(
    gene_set = p,
    p_value = test$p.value,
    odds_ratio = test$estimate,
    count_in_set = sum(x$pathway == T),
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

calculate_enrichment_panel <- function(p, df2) {
  x <- data.frame(panel = grepl(pattern = p, x = df2$Panel), associated_hormones = df2$`associated_p1e-3`)
  contingency_table <- table(
    x$panel == TRUE,         # In panel?
    x$associated_hormones == 1  # Associated with Hormone?
  )
  if (all(dim(contingency_table) == 2)) {
    
    test <- fisher.test(contingency_table)
    data.frame(
      gene_set = p,
      p_value = test$p.value,
      odds_ratio = test$estimate,
      count_in_set = sum(x$panel == T),
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

enrichment_results_pathways <- lapply(all_pathways, calculate_enrichment_human_pathways, df2) %>%
  bind_rows() %>%
  mutate(
    fdr = ifelse(overlap_count > 0, p.adjust(p_value, method = "BH"), NA)
  )

enrichment_results_panel <- lapply(all_panels_category, calculate_enrichment_panel, df2) %>%
  bind_rows() %>%
  mutate(
    fdr = ifelse(overlap_count > 0, p.adjust(p_value, method = "BH"), NA)
  )


#### Are estimates larger for proteins in INF and CVD than for the rest
d1 <- read.delim("results12/intensity_shared_prots_261125/prot_vs_hormones.gam.spline.shared_prots.txt", as.is = T, check.names  = F, sep = "\t")
d2 <- read.delim("results12/intensity_batch2_prots_261125/prot_vs_hormones.gam.spline.shared_prots.txt", as.is = T, check.names  = F, sep = "\t")
d1$batch = "batch12"
d2$batch <- "batch2"
d <- rbind(d1[,c("prot", "pheno", "pval", "estimate","batch","BH_pval")], d2[,c("prot", "pheno", "pval", "estimate","batch","BH_pval")])
d$abs_estimate <- abs(d$estimate)

df2 <- read_csv("data/all_tested_2453_olink_annot.csv")
all_panels_category <- na.omit(unique(gsub("_.*","",df2$Panel)))

d_annot <- left_join(d, df2[, c("Gene", "Panel")], by = c("prot" = "Gene"))

res_table_panels <- data.frame(Panel = all_panels_category, KS_twosided = NA, KS_greater = NA, KS_less = NA, Wilcox_twosided = NA)
for (p in all_panels_category){
  d_annot$panel_bin <- grepl(pattern = p, x = d_annot$Panel)
  cat(p, "\n")
  twosided <- ks.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"])
  greater <- ks.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"], alternative = "greater")
  less <- ks.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"], alternative = "less")
  wilcoxon <- wilcox.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"])
  
  res_table_panels[res_table_panels$Panel == p,2:5] <- c(twosided$p.value, greater$p.value, less$p.value, wilcoxon$p.value)
}


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
  left_join(df2, by = "Gene") %>%
  select(Gene, `associated_p1e-3`, category) %>%
  unique() 

df_with_cat_wide <- df_with_cat %>%
  group_by(Gene, `associated_p1e-3`) %>%
  summarise(
    Panel = toString(category),
    .groups = "drop"
  )



enrichment_results_UKB <- lapply(na.omit(unique(df_with_cat$category)), calculate_enrichment_panel, df_with_cat_wide) %>%
  bind_rows() %>%
  mutate(
    p_bonf = ifelse(overlap_count > 0, p.adjust(p_value, method = "bonferroni"), NA)
  )

df_with_dis <- df_long %>%
  mutate(disease_clean_lower = str_trim(tolower(disease_clean))) %>%
  left_join(dis_cat_clean, by = c("disease_clean_lower" = "phenotype_clean")) %>%
  left_join(df2, by = "Gene") %>%
  select(Gene, `associated_p1e-3`, disease_clean) %>%
  unique() 

df_with_dis_wide <- df_with_dis %>%
  group_by(Gene, `associated_p1e-3`) %>%
  summarise(
    Panel = toString(disease_clean),
    .groups = "drop"
  )
enrichment_results_UKB_dis <- lapply(unique(df_with_dis$disease_clean), calculate_enrichment_panel, df_with_dis_wide) %>%
  bind_rows() %>%
  mutate(
    qval_BH = ifelse(overlap_count > 0, p.adjust(p_value, method = "BH"), NA)
  )


# Are estimates higher for proteins associated with any specific UKB disease category?
d_annot <- left_join(d, df_with_cat_wide[, c("Gene", "Panel")], by = c("prot" = "Gene"))
all_panels_category <- na.omit(unique(df_with_cat$category))
res_table_UKB <- data.frame(Panel = all_panels_category, KS_twosided = NA, KS_greater = NA, KS_less = NA, Wilcox_twosided = NA)

for (p in all_panels_category){
  d_annot$panel_bin <- grepl(pattern = p, x = d_annot$Panel)
  cat(p, "\n")
  twosided <- ks.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"])
  greater <- ks.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"], alternative = "greater")
  less <- ks.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"], alternative = "less")
  wilcoxon <- wilcox.test(d_annot[d_annot$panel_bin == T,"abs_estimate"], d_annot[d_annot$panel_bin == F,"abs_estimate"])
  
  res_table_UKB[res_table_UKB$Panel == p,2:5] <- c(twosided$p.value, greater$p.value, less$p.value, wilcoxon$p.value)
}
View(res_table_UKB)

# per hormone separately
res_table_UKB_per_hormone <- expand.grid(
  Panel = all_panels_category,
  hormone = unique(d_annot$pheno), 
  KS_twosided = NA, KS_greater = NA, KS_less = NA, Wilcox_twosided = NA,
  stringsAsFactors = FALSE)

for (ph in unique(d_annot$pheno)){
  d_subs <- d_annot[d_annot$pheno == ph,]
  for (p in all_panels_category){
    d_subs$panel_bin <- grepl(pattern = p, x = d_subs$Panel)
    cat(p, "\n")
    twosided <- ks.test(d_subs[d_subs$panel_bin == T,"abs_estimate"], d_subs[d_subs$panel_bin == F,"abs_estimate"])
    greater <- ks.test(d_subs[d_subs$panel_bin == T,"abs_estimate"], d_subs[d_subs$panel_bin == F,"abs_estimate"], alternative = "greater")
    less <- ks.test(d_subs[d_subs$panel_bin == T,"abs_estimate"], d_subs[d_subs$panel_bin == F,"abs_estimate"], alternative = "less")
    wilcoxon <- wilcox.test(d_subs[d_subs$panel_bin == T,"abs_estimate"], d_subs[d_subs$panel_bin == F,"abs_estimate"])
    
    res_table_UKB_per_hormone[res_table_UKB_per_hormone$Panel == p & res_table_UKB_per_hormone$hormone == ph,3:6] <- c(twosided$p.value, greater$p.value, less$p.value, wilcoxon$p.value)
  }
}
View(res_table_UKB_per_hormone)
