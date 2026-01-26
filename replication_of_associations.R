library(readxl)
data <- read_excel("/Users/Dasha/work/Sardinia/W4H/olink/batch12/results12/intensity_shared_prots_261125/replication_BH_pval.xlsx", sheet = 1,  skip = 1)
data <- data[!is.na(data$p.adj_PROG),]
data$estimate_Dordevic <- NULL
data$p.adj_Dordevic <- NULL
data$estimate_MHT <- -1*data$estimate_MHT

est_columns <- c(colnames(data)[grepl("estimate", colnames(data))], "logFC_Tarca")
est_columns <- est_columns[! est_columns %in% c("estimate_PROG", "estimate_17BES")]

pval_columns <- colnames(data)[grepl("p.adj", colnames(data))]
pval_columns <- pval_columns[! pval_columns %in% c("p.adj_PROG", "p.adj_17BES")]


# Progesterone
subs <- data[data$p.adj_PROG < 0.05,]
nrow(subs)

subs$replicated_pval <- rowSums(subs[, pval_columns] < 0.05, na.rm = TRUE) > 0

subs$replicated_pval_and_direction <- FALSE

# Loop through each study
for(i in c(4, 5,6)) {
  # Get current study's p-value and estimate columns
  pval_col <- pval_columns[i]
  est_col <- est_columns[i]
  
  # Check: p-value < 0.05 AND estimate has same sign as estimate_PROG
  condition <- !is.na(subs[[pval_col]]) & 
    subs[[pval_col]] < 0.05 & 
    sign(subs[[est_col]]) == sign(subs$estimate_PROG)
  subs[,paste0("repl_dir_", i)] <- condition
  # Update replicated_pval: TRUE if any study meets condition
  subs$replicated_pval_and_direction <- subs$replicated_pval_and_direction | condition
}

table(subs$replicated_pval)
table(subs$replicated_pval_and_direction)


# Estrogen
subs <- data[data$p.adj_17BES < 0.05,]
nrow(subs)

subs$replicated_pval <- rowSums(subs[, pval_columns] < 0.05, na.rm = TRUE) > 0

subs$replicated_pval_and_direction <- FALSE

# Loop through each study
#for(i in seq_along(pval_columns)) {
i = 3
  # Get current study's p-value and estimate columns
pval_col <- pval_columns[i]
est_col <- est_columns[i]

# Check: p-value < 0.05 AND estimate has same sign as estimate_PROG
condition <- !is.na(subs[[pval_col]]) & 
  subs[[pval_col]] < 0.05 & 
  sign(subs[[est_col]]) == sign(subs$estimate_17BES)
subs[,paste0("repl_dir_", i)] <- condition
# Update replicated_pval: TRUE if any study meets condition
subs$replicated_pval_and_direction <- subs$replicated_pval_and_direction | condition
#}

table(subs$replicated_pval)
table(subs$replicated_pval_and_direction)
