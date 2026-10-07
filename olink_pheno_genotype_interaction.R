my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/utility_functions.R")

library(ggplot2)
library(ggExtra)

set.seed(123)

out_basedir <- "results12/intensity_all_prots_220526/"

# READ PROTEINS
d_wide <- read.delim(paste0(out_basedir, "olink_batch12.all_proteins.phase_avg.adj_batch_storage.txt"), as.is = T, check.names = F, sep = "\t")

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$phase <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    dplyr::select(SampleID, ID, phase, everything())
}

d_wide$phase <- relevel(factor(d_wide$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
d_wide$TP <- as.numeric(d_wide$phase)

all_phases <- c("F", "O", "EL", "LL")
all_prots <- colnames(d_wide)[! colnames(d_wide) %in% c("SampleID", "ID", "TP","phase")]


# READ COVARIATES
covariates <- read.delim("results12/covariates_olink_batch12.phase_avg.txt", sep = "\t", check.names = F, as.is = T)

if (! "ID" %in% colnames(covariates)){
  covariates$ID <- gsub("_.*", "", covariates$SampleID)
  covariates$phase <- gsub(".*_", "", covariates$SampleID)
  covariates <- covariates %>%
    dplyr::select(SampleID, ID, phase, everything())
}
covariates$phase <- relevel(factor(covariates$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
covariates$TP <- as.numeric(covariates$phase)

covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})


covariate_names <- c("Age","BMI", "from")
covariates$from <- relevel(as.factor(covariates$from), ref = "X")
covariates$storage_months = covariates$batch  = NULL

# READ PHENOTYPES
pheno <- read.delim("../../phenotypes/batch12/cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.phase_avg.txt", as.is = T, check.names = F, sep = "\t")            

all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
all_phenos <- c("INS", "HOMA_B", "HOMA_IR", "GL", "AST", "ALT", "TRI", "COL", "HDL", "LDL")
all_phenos_combined <- c(all_hormones, all_phenos)
covariate_names_pheno <- c("from", "Age", "BMI")

if (! "ID" %in% colnames(pheno)){
  pheno$ID <- gsub("_.*", "", pheno$SampleID)
  pheno$phase <- gsub(".*_", "", pheno$SampleID)
  pheno <- pheno %>%
    dplyr::select(SampleID, ID, phase, everything())
}

pheno$phase <- relevel(factor(pheno$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
pheno$TP <- as.numeric(pheno$phase)

# READ PRS
prs <- read.delim("data/merged_protein_PRS.tsv", sep = "\t", as.is = T, check.names = F)
prs[grepl("^[0-9]",prs$IID), "IID"] <- paste0("X", prs[grepl("^[0-9]",prs$IID), "IID"])

prs <- prs[, colSums(is.na(prs)) <= 100]

# READ SIGNIFICANT ASSOCIATIONS
gam_res_hormones <- read.delim(paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.txt"), as.is = T, sep = "\t", check.names = FALSE)
gam_res_phenos<- read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.txt"), as.is = T, sep = "\t", check.names = FALSE)

#signif_pairs <- rbind(
#  gam_res_phenos[gam_res_phenos$BH_pval < 0.05, c("prot", "pheno")],
#  gam_res_hormones[gam_res_hormones$BH_pval < 0.05, c("prot", "pheno")]
#)
signif_pairs = gam_res_hormones[gam_res_hormones$BH_pval < 0.05, c("prot", "pheno")]

signif_pairs <- signif_pairs[signif_pairs$prot %in% colnames(prs),]
signif_pairs$gam_inter_p <- NA
signif_pairs$lmm_inter_p <- NA

# RUN PROTEIN - HORMONE/PHENO ASSOCIATIONS

for (i in 1:nrow(signif_pairs)){
  prot = signif_pairs[i, ]$prot
  ph = signif_pairs[i, ]$pheno
  cur_covars <- covariate_names
  
  d_subs <- d_wide[,c("SampleID", "TP", "ID", prot)] %>%
    inner_join(pheno[,c("SampleID", "TP", "ID", ph)], by = c("SampleID", "TP", "ID")) %>%
    inner_join(covariates, by = c("SampleID", "TP", "ID")) %>%
    inner_join(prs %>% transmute(ID = IID, prs = as.numeric(scale(!!sym(prot)))), by = "ID") %>%
    rename(
      prot = !!prot, 
      pheno = !!ph
    ) %>%
    mutate(
      prot = as.numeric(scale(prot)),
      pheno = as.numeric(scale(pheno)),
      ID = as.factor(ID)
    ) %>%
    drop_na()
 if (length(unique(d_subs$batch)) == 1) cur_covars <- cur_covars[cur_covars != "batch"]
  
  fo_gam <- as.formula(paste("prot ~ pheno*prs + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(cur_covars, collapse = "+")))
  fo_lmm <- as.formula(paste("prot ~ pheno*prs + poly(TP, 3) +", paste(cur_covars, collapse = "+"), "+ (1|ID)"))
  
  gam_fit <- gam(fo_gam, data = d_subs, method = 'REML')
  lmm_fit <- lmer(fo_lmm, data = d_subs)

  signif_pairs[i,]$gam_inter_p = summary(gam_fit)$p.table["pheno:prs","Pr(>|t|)"]
  signif_pairs[i,]$lmm_inter_p <- summary(lmm_fit)$coefficients["pheno:prs","Pr(>|t|)"]
}
write.table(signif_pairs, file = paste0(out_basedir, "genetics/prot_vs_phenotypes.interaction_PRS.txt"), , sep = "\t", quote = F, row.names = F)

prot = 'MLN'
ph = 'TRI'

d_subs <- d_wide[,c("SampleID", "TP", "ID", prot)] %>%
  inner_join(pheno[,c("SampleID", "TP", "ID", ph)], by = c("SampleID", "TP", "ID")) %>%
  inner_join(covariates, by = c("SampleID", "TP", "ID")) %>%
  inner_join(prs %>% transmute(ID = IID, prs = as.numeric(scale(!!sym(prot)))), by = "ID") %>%
  rename(
    prot = !!prot, 
    pheno = !!ph
  ) %>%
  mutate(
    prot = as.numeric(scale(prot)),
    pheno = as.numeric(scale(pheno)),
    ID = as.factor(ID)
  ) %>%
  drop_na()

d_subs <- d_subs %>%
  mutate(
    prs_cat = ntile(prs, 3), # This returns 1, 2, or 3 for every row
    prs_cat = factor(prs_cat, labels = c("Low", "Middle", "High"))
  )
ggplot(d_subs, aes(x = pheno, y = prot, color = prs_cat, group = prs_cat)) + geom_point() + geom_smooth(method = 'lm') + theme_minimal() + xlab("TRI") + ylab("MLN")

####
## PROTEIN - PHASE ASSOCIATIONS
#####
# READ SIGNIFICANT ASSOCIATIONS
gam_res <- read.delim(paste0(out_basedir, "prot_vs_tp_gam.txt"), as.is = T, sep = "\t", check.names = FALSE)

signif_prots <- gam_res[gam_res$gam_BH_pval < 0.05, "prot", drop = F]
signif_prots <- signif_prots[signif_prots$prot %in% colnames(prs),, drop = F]
signif_prots$lmm_inter_p <- NA

# RUN INTERACTION ANALYSIS

for (i in 1:nrow(signif_prots)){
  prot = signif_prots[i, ]$prot
  cur_covars <- covariate_names
  
  d_subs <- d_wide[,c("SampleID", "TP", "ID", prot)] %>%
    inner_join(covariates, by = c("SampleID", "TP", "ID")) %>%
    inner_join(prs %>% transmute(ID = IID, prs = as.numeric(scale(!!sym(prot)))), by = "ID") %>%
    rename(
      prot = !!prot
    ) %>%
    mutate(
      prot = as.numeric(scale(prot)),
      ID = as.factor(ID)
    ) %>%
    drop_na()
  if (length(unique(d_subs$batch)) == 1) cur_covars <- cur_covars[cur_covars != "batch"]
  
  fo_lmm <- as.formula(paste("prot ~ TP*prs + ", paste(cur_covars, collapse = "+"), "+ (1|ID)"))
  lmm_fit <- lmer(fo_lmm, data = d_subs)
  
  signif_prots[i,]$lmm_inter_p <- summary(lmm_fit)$coefficients["TP:prs","Pr(>|t|)"]

}
write.table(signif_prots, file = paste0(out_basedir, "genetics/prot_vs_phase.lmm_linear_interaction_PRS.txt"), , sep = "\t", quote = F, row.names = F)



# TEST ALL PQTLS

pqtls <- read.delim("../../genotypes/batch14/all_clumped_snp_proteins.with_genotypes.txt.gz", sep = "\t", as.is = T, check.names = F)
pqtls <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/E2_prot_qtls.ESR1_peaks.with_genotypes.txt.gz", sep = "\t", as.is = T, check.names = F)

prot = "PROK1"
ph = "PROG"
cur_covars <- c(covariate_names, "pqtl*pheno")

res_table2 <- data.frame(matrix(nrow = nrow(signif_pairs) * nrow(pqtls[pqtls$prot %in% signif_pairs$prot,]), ncol = 5))
cnt <- 1

signif_pairs_subs <- signif_pairs[signif_pairs$prot %in% pqtls$prot & signif_pairs$pheno %in% c("PROG", "17BES"),]
for (i in 1:nrow(signif_pairs_subs)){
  prot = signif_pairs_subs[i, ]$prot
  ph = signif_pairs_subs[i, ]$pheno

  genos <- pqtls[pqtls$prot == prot,]
  genos <- genos[order(genos$pval),]
  row.names(genos) <- genos$SNP
  genos[,c("SNP","prot", "effect_allele", "other_allele", "beta","pval", "CHR","POS")] <-NULL
  genos <- as.data.frame(t(genos)) %>%
    rownames_to_column("ID")
  genos$ID <- gsub("T", "X", genos$ID)
  
  for (pqtl in colnames(genos)[2:ncol(genos)]){
    if (min(table(genos[,pqtl])) < 10) next
    
    d_subs <- d_wide[,c("SampleID", "TP", "ID", prot)] %>%
      inner_join(pheno[,c("SampleID", "TP", "ID", ph)], by = c("SampleID", "TP", "ID")) %>%
      inner_join(covariates, by = c("SampleID", "TP", "ID")) %>%
      inner_join(genos[,c("ID", pqtl)], by = "ID") %>%
      rename(
        prot = !!prot, 
        pheno = !!ph,
        pqtl = !!pqtl
      ) %>%
      mutate(
        prot = as.numeric(scale(prot)),
        pheno = as.numeric(scale(pheno)),
        ID = as.factor(ID)
      ) %>%
      drop_na()
    
    #fo_gam <- as.formula(paste("prot ~ pheno + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(cur_covars, collapse = "+")))
    #fo_lmm <- as.formula(paste("prot ~ pheno + TP +", paste(cur_covars, collapse = "+"), "+ (1|ID)"))
    
    fo_gam <- as.formula(paste("prot ~ pheno + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(cur_covars, collapse = "+")))
    fo_lmm <- as.formula(paste("prot ~ pheno + TP +", paste(cur_covars, collapse = "+"), "+ (1|ID)"))
    
    gam_fit <- gam(fo_gam, data = d_subs, method = 'REML')
    lmm_fit <- lmer(fo_lmm, data = d_subs)
    
    gam_pval <- summary(gam_fit)$p.table["pheno:pqtl", "Pr(>|t|)"]
    lmm_pval <- summary(lmm_fit)$coefficients["pheno:pqtl","Pr(>|t|)"]
    
    res_table2[cnt,] <- c(ph, prot, pqtl, gam_pval, lmm_pval)
    cnt <- cnt + 1
  }
}
colnames(res_table2) <- c("pheno", "prot", "pQTL", "inter_gam_pval", "inter_lmm_pval")

res_table2 <- res_table2 %>%
  na.omit(res_table2) %>%
  mutate(across(-c(pheno, prot, pQTL), as.numeric)) 
write.table(res_table2, file = paste0(out_basedir, "pqtl_interaction/tmp_pqtl_hormone_interaction_p4_e2_minAC10.txt"), quote = F, sep = "\t", row.names = FALSE)



## PLOTs

prot = "PROK1"
ph = 'PROG'
pqtl <- 'rs4839011'

genos <- pqtls[pqtls$prot == prot,]
#genos <- genos[order(genos$pval),]
row.names(genos) <- genos$SNP
genos[,c("SNP","prot", "effect_allele", "other_allele", "beta","pval", "CHR","POS")] <-NULL
genos <- as.data.frame(t(genos)) %>%
  rownames_to_column("ID")
genos$ID <- gsub("T", "X", genos$ID)


d_subs <- d_wide[,c("SampleID", "TP", "ID", prot)] %>%
    inner_join(pheno[,c("SampleID", "TP", "ID", ph)], by = c("SampleID", "TP", "ID")) %>%
    inner_join(covariates, by = c("SampleID", "TP", "ID")) %>%
    inner_join(genos[,c("ID", pqtl)], by = "ID") %>%
    rename(
      prot = !!prot, 
      pheno = !!ph,
      pqtl = !!pqtl
    ) %>%
    mutate(
      prot = as.numeric(scale(prot)),
      pheno = as.numeric(scale(pheno)),
      ID = as.factor(ID),
      pqtl_factor = as.factor(pqtl)
    ) %>%
    drop_na()

ggplot(d_subs, aes(x = prot, y = pheno, color = pqtl_factor)) + 
  geom_point() + geom_smooth(method = 'lm') + theme_minimal() +
  xlab(prot) + ylab(ph)


##################


#### Remap + Fabian SNPs

# E2
pqtls_tmp <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/E2_prot_qtls._with_pos.txt", sep = "\t", as.is = T, check.names = F, col.names = c("prot", "SNP", "chr", "pos"))
pqtls_tmp$chrpos <- paste0(pqtls_tmp$chr, ":", pqtls_tmp$pos)

esr1_peaks <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/remap_overlap.ESR1.my.txt", sep = "\t", as.is = T, check.names = F, , col.names = c("SNP", "ESR1_celltype"))
esr2_peaks <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/remap_overlap.ESR2.my.txt", sep = "\t", as.is = T, check.names = F , col.names = c("SNP", "ESR2_celltype"))

esr1_peaks <- esr1_peaks %>%
  mutate(ESR1_celltype = str_remove(ESR1_celltype, "^ESR1:")) %>%
  group_by(SNP) %>%
  summarise(ESR1_celltype = paste(unique(ESR1_celltype), collapse = ","), .groups = "drop")

esr2_peaks <- esr2_peaks %>%
  mutate(ESR2_celltype = str_remove(ESR2_celltype, "^ESR2:")) %>%
  group_by(SNP) %>%
  summarise(ESR2_celltype = paste(unique(ESR2_celltype), collapse = ","), .groups = "drop")

esr_peaks <- full_join(esr1_peaks, esr2_peaks, by = "SNP")

fabian <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/W4H_batch14.snps.b38.E2_prots_snps.fabian_res.tsv", sep = "\t", as.is = T, check.names = F, 
                     col.names = c("variant", "tf", "model_id",	"database", "ref_score", "alt_score", "start_ref",	"end_ref","start_alt","end_alt",	"strand_ref",	"strand_alt",	"prediction",	"score"))
fabian$max_ref_alt_score <- max(fabian$ref_score, fabian$alt_score)
fabian$abs_score <- abs(fabian$score)
fabian$chrpos <- gsub("[ACGT]>.*","", fabian$variant)

fabian_esr1_snps <- fabian[fabian$max_ref_alt_score > 0.8 & fabian$abs_score > 0.6 & fabian$tf == 'ESR1', "chrpos"]
fabian_esr2_snps <- fabian[fabian$max_ref_alt_score > 0.8 & fabian$abs_score > 0.6 & fabian$tf == 'ESR2', "chrpos"]

pqtls_annot <- left_join(pqtls_tmp, esr_peaks, by = "SNP")
pqtls_annot$ESR1_motif_change <- ifelse(pqtls_annot$chrpos %in% fabian_esr1_snps, T, NA)
pqtls_annot$ESR2_motif_change <- ifelse(pqtls_annot$chrpos %in% fabian_esr2_snps, T, NA)

genotypes <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/W4H_batch14.snps.b38.E2_prots_snps.raw.genotypes.txt.gz", sep = "\t", as.is = T, check.names = F)

pqtls_annot[,c("chr", "pos", "chrpos")] <- NULL
pqtls_annot <- pqtls_annot %>%
  filter(!if_all(3:6, is.na))

pqtls <- left_join(pqtls_annot[,c(1,2)], genotypes, by = "SNP")


# P4
pqtls_tmp <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/P4_assoc_prots/P4_prot_qtls._with_pos.txt", sep = "\t", as.is = T, check.names = F, col.names = c("prot", "SNP", "chr", "pos"))
pqtls_tmp$chrpos <- paste0(pqtls_tmp$chr, ":", pqtls_tmp$pos)

pgr_peaks <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/P4_assoc_prots/remap_overlap.PGR.my.txt", sep = "\t", as.is = T, check.names = F, , col.names = c("SNP", "PGR_celltype"))

pgr_peaks <- pgr_peaks %>%
  mutate(PRG_celltype = str_remove(PGR_celltype, "^PRG:")) %>%
  group_by(SNP) %>%
  summarise(PGR_celltype = paste(unique(PGR_celltype), collapse = ","), .groups = "drop")


fabian <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/P4_assoc_prots/W4H_batch14.snps.b38.P4_prots_snps.fabian_res.tsv", sep = "\t", as.is = T, check.names = F, 
                     col.names = c("variant", "tf", "model_id",	"database", "ref_score", "alt_score", "start_ref",	"end_ref","start_alt","end_alt",	"strand_ref",	"strand_alt",	"prediction",	"score"))
fabian$max_ref_alt_score <- max(fabian$ref_score, fabian$alt_score)
fabian$abs_score <- abs(fabian$score)
fabian$chrpos <- gsub("[ACGT]>.*","", fabian$variant)

fabian_pgr_snps <- fabian[fabian$max_ref_alt_score > 0.8 & fabian$abs_score > 0.6 & fabian$tf == 'PGR', "chrpos"]

pqtls_annot <- left_join(pqtls_tmp, pgr_peaks, by = "SNP")
pqtls_annot$PGR_motif_change <- ifelse(pqtls_annot$chrpos %in% fabian_pgr_snps, T, NA)

genotypes <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/HRE/P4_assoc_prots/W4H_batch14.snps.b38.P4_prots_snps.raw.genotypes.txt.gz", sep = "\t", as.is = T, check.names = F)

pqtls_annot[,c("chr", "pos", "chrpos")] <- NULL
pqtls_annot <- pqtls_annot %>%
  filter(!if_all(3:4, is.na))

pqtls <- left_join(pqtls_annot[,c(1,2)], genotypes, by = "SNP")



cur_covars <- c(covariate_names, "pqtl*pheno")

signif_pairs_subs <- signif_pairs[signif_pairs$prot %in% pqtls$prot & signif_pairs$pheno == "PROG",]

res_table_p4 <- data.frame(matrix(nrow = nrow(signif_pairs_subs) * nrow(pqtls[pqtls$prot %in% signif_pairs_subs$prot,]), ncol = 5))
cnt <- 1

for (i in 1:nrow(signif_pairs_subs)){
  prot = signif_pairs_subs[i, ]$prot
  ph = signif_pairs_subs[i, ]$pheno
  
  if (! prot %in% pqtls$prot) next
  
  genos <- pqtls[pqtls$prot == prot,]
  
  row.names(genos) <- genos$SNP
  genos[,c("SNP","prot", "effect_allele", "other_allele", "beta","pval", "CHR","POS")] <-NULL
  genos <- as.data.frame(t(genos)) %>%
    rownames_to_column("ID")
  genos$ID <- gsub("T", "X", genos$ID)
  
  for (pqtl in colnames(genos)[2:ncol(genos)]){
    if (min(table(genos[,pqtl])) < 10) next
    
    d_subs <- d_wide[,c("SampleID", "TP", "ID", prot)] %>%
      inner_join(pheno[,c("SampleID", "TP", "ID", ph)], by = c("SampleID", "TP", "ID")) %>%
      inner_join(covariates, by = c("SampleID", "TP", "ID")) %>%
      inner_join(genos[,c("ID", pqtl)], by = "ID") %>%
      rename(
        prot = !!prot, 
        pheno = !!ph,
        pqtl = !!pqtl
      ) %>%
      mutate(
        prot = as.numeric(scale(prot)),
        pheno = as.numeric(scale(pheno)),
        ID = as.factor(ID)
      ) %>%
      drop_na()
    
    #fo_gam <- as.formula(paste("prot ~ pheno + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(cur_covars, collapse = "+")))
    #fo_lmm <- as.formula(paste("prot ~ pheno + TP +", paste(cur_covars, collapse = "+"), "+ (1|ID)"))
    
    fo_gam <- as.formula(paste("prot ~ pheno + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(cur_covars, collapse = "+")))
    fo_lmm <- as.formula(paste("prot ~ pheno + TP +", paste(cur_covars, collapse = "+"), "+ (1|ID)"))
    
    gam_fit <- gam(fo_gam, data = d_subs, method = 'REML')
    lmm_fit <- lmer(fo_lmm, data = d_subs)
    
    gam_pval <- summary(gam_fit)$p.table["pheno:pqtl", "Pr(>|t|)"]
    lmm_pval <- summary(lmm_fit)$coefficients["pheno:pqtl","Pr(>|t|)"]
    
    res_table_p4[cnt,] <- c(ph, prot, pqtl, gam_pval, lmm_pval)
    cnt <- cnt + 1
  }
}
colnames(res_table_p4) <- c("pheno", "prot", "pQTL", "inter_gam_pval", "inter_lmm_pval")

res_table_p4 <- res_table_p4 %>%
  na.omit(res_table_p4) %>%
  mutate(across(-c(pheno, prot, pQTL), as.numeric)) 
write.table(res_table_p4, file = paste0(out_basedir, "pqtl_interaction/PGR_pqtl_hormone_interaction.remap_fabian_snps.txt"), quote = F, sep = "\t", row.names = FALSE)



