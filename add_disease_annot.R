annot <- read.csv("data/all_proteins_annotated.csv")
prot_names <- read.delim("data/batch12_protein_names.txt",  sep = "\t", as.is = T, check.names = F)

annot <- left_join(prot_names[,c("Assay", "UniProt")], annot, by = c("Assay" = "hgnc_symbol"))

mr <- read.delim("../papers_clean/Su_MR.txt", sep = "\t", as.is = T, check.names = F)
mr_joined <-mr %>%
  group_by(UniProt) %>%
  summarise(MR_outcomes = paste(unique(`Outcome name`), collapse = ", "))

mr_cut <- inner_join(annot[,c("Assay", "UniProt")], mr[,c("UniProt", "Outcome name", "Beta")], by = "UniProt")
colnames(mr_cut)[3] <- "MR_outcomes"
write.table(mr_cut, file = "../papers_clean/Su_MR.subset_cut.txt", quote = F, sep = "\t", row.names = FALSE)


coloc <- read.delim("../papers_clean/Pietzner_coloc.txt", sep = "\t", as.is = T, check.names = F)
coloc$EntrezGeneSymbol <- gsub("^'", "", coloc$EntrezGeneSymbol)
coloc_joined <-coloc %>%
  separate_rows(EntrezGeneSymbol, sep = "\\s+") %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(coloc_traits = paste(unique(`mapped_trait`), collapse = ", "))

annot2 <- left_join(annot, mr_joined, by = "UniProt")
annot2 <- left_join(annot2, coloc_joined, by = c("Assay" = "EntrezGeneSymbol"))

annot2$combined_annotations <- NULL
annot2$ensembl_gene_id <- NULL
annot2$chromosome_name <- NULL

write.table(annot2, file = "data/all_proteins_annotated_diseases.txt", quote = F, sep = "\t", row.names = FALSE)

