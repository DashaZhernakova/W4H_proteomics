d1 <- read.delim("results12/intensity_shared_prots_261125/prot_vs_hormones.gam.spline.shared_prots.txt", as.is = T, check.names  = F, sep = "\t")
d2 <- read.delim("results12/intensity_batch2_prots_261125/prot_vs_hormones.gam.spline.shared_prots.txt", as.is = T, check.names  = F, sep = "\t")
d1$batch = "batch12"
d2$batch <- "batch2"
d <- rbind(d1[,c("prot", "pheno", "pval", "batch","BH_pval")], d2[,c("prot", "pheno", "pval", "batch","BH_pval")])
d$BH_pval_joint <- p.adjust(d$pval, method = 'BH')

d$sign <- ifelse(d$BH_pval < 0.05, T, F)
d$sign_joint <- ifelse(d$BH_pval_joint < 0.05, T, F)

table(d[,c("sign", "sign_joint", "batch")])


d[d$sign == T & d$sign_joint == F,]
