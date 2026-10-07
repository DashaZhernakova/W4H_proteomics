
setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
source("/Users/Dasha/work/Sardinia/W4H/olink/scripts/utility_functions.R")

set.seed(123)

out_basedir <- "results12/intensity_all_prots_220526/"

d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.phase_avg.txt", as.is = T, check.names = F, sep = "\t")

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$phase <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    dplyr::select(SampleID, ID, phase, everything())
}

d_wide$phase <- relevel(factor(d_wide$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

covariates <- read.delim("results12/covariates_olink_batch12.phase_avg.txt", sep = "\t", check.names = F, as.is = T)

if (! "ID" %in% colnames(covariates)){
  covariates$ID <- gsub("_.*", "", covariates$SampleID)
  covariates$phase <- gsub(".*_", "", covariates$SampleID)
  covariates <- covariates %>%
    dplyr::select(SampleID, ID, phase, everything())
}
covariates$phase <- relevel(factor(covariates$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})


# Make a dataframe with proteins adjusted for all covariates per visit
covariate_names <- c("Age","BMI","batch", "from", "storage_months")
covariates$from <- relevel(as.factor(covariates$from), ref = "X")


all_phases <- c("F", "O", "EL", "LL")
all_prots <- colnames(d_wide)[! colnames(d_wide) %in% c("SampleID", "ID", "TP","phase")]

joined_data <- full_join(covariates, d_wide, by = c("SampleID", "ID", "phase"))
d_wide_adj_covar <- read.delim(paste0(out_basedir, "olink_batch12.all_proteins.phase_avg.adj_all_covariates.txt"), as.is = T, check.names = F, sep = "\t")
d_wide_adj_covar$phase <- relevel(factor(d_wide_adj_covar$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

# mean protein abundance
# calculate mean abundance per phase for all proteins
mean_abund_per_phase_adj_covar_long <- d_wide_adj_covar %>%
  pivot_longer(cols = all_of(all_prots), 
               names_to = "Protein", 
               values_to = "Level") %>%
  group_by(Protein, phase) %>%
  summarise(mean_prot = mean(Level, na.rm = TRUE), .groups = "drop")

mean_abund_per_phase_adj_covar <- mean_abund_per_phase_adj_covar_long %>%
  pivot_wider(names_from = phase, values_from = mean_prot)


# limma
protein_batch_info <- lapply(all_prots, function(prot) {
  observed_batches <- unique(joined_data$batch[!is.na(joined_data[[prot]])])
  data.frame(protein = prot, n_batches = length(observed_batches),  measured_batches = paste(observed_batches, collapse = ","))
}) %>%
  bind_rows()
covariates_no_batch <- setdiff(covariate_names, "batch")

table(protein_batch_info$measured_batches)

protein_groups <- split(protein_batch_info$protein, protein_batch_info$measured_batches)
base_covariates <- setdiff(covariate_names, "batch")

phase_comb <- t(combn(all_phases, 2))

all_results <- list()
result_index <- 1

for (i in seq_len(nrow(phase_comb))) {
  
  tp1 <- phase_comb[i, 1]
  tp2 <- phase_comb[i, 2]
  
  comparison_results <- list()
  group_index <- 1
  
  for (pattern in names(protein_groups)) {
    
    target_proteins <- protein_groups[[pattern]]
    batches <- strsplit(pattern, ",", fixed = TRUE)[[1]]
    
    # Proteins measured in every batch used by this model
    fit_proteins <- protein_batch_info$protein[
      vapply(
        protein_batch_info$measured_batches,
        function(x) {
          measured <- strsplit(x, ",", fixed = TRUE)[[1]]
          all(batches %in% measured)
        },
        logical(1)
      )
    ]
    
    analysis_data <- joined_data %>%
      filter(batch %in% batches) %>%
      mutate(batch = droplevels(factor(batch)))
    
    model_covariates <- if (length(batches) > 1) {
      c(base_covariates, "batch")
    } else {
      base_covariates
    }
    
    limma_res <- run_limma(
      joined_data = analysis_data,
      tp1 = tp1,
      tp2 = tp2,
      all_prots = fit_proteins,
      covariate_names = model_covariates
    ) %>%
      rownames_to_column(var = "prot") %>%
      filter(prot %in% target_proteins) %>%
      mutate(measured_batches = pattern)
    
    comparison_results[[group_index]] <- limma_res
    group_index <- group_index + 1
  }
  
  comparison_results <- bind_rows(comparison_results) %>%
    mutate(
      tp1 = tp1,
      tp2 = tp2,
      phase1_phase2 = paste(tp1, tp2, sep = "_"),
      # Recalculate BH across all proteins for this comparison
      adj.P.Val = p.adjust(P.Value, method = "BH")
    )
  
  all_results[[result_index]] <- comparison_results
  result_index <- result_index + 1
}

limma_res_all <- bind_rows(all_results)
limma_res_all$adj.P.Val <- p.adjust(limma_res_all$P.Value, method = 'BH')
limma_res_all$sign <- ifelse(limma_res_all$adj.P.Val < 0.05,T,F)
table(limma_res_all$sign)

limma_res_with_means <- limma_res_all %>%
  separate(phase1_phase2, into = c("phase1", "phase2"), sep = "_", remove = FALSE) %>%
  left_join(mean_abund_per_phase_adj_covar_long, by = c("prot" = "Protein", "phase1" = "phase")) %>%
  rename(phase1_mean = mean_prot) %>%
  left_join(mean_abund_per_phase_adj_covar_long, by = c("prot" = "Protein", "phase2" = "phase")) %>%
  rename(phase2_mean = mean_prot) %>%
  select(-phase1, -phase2)

length(unique(limma_res_all[limma_res_all$sign == T, "prot"]))
write.table(limma_res_with_means, file = paste0(out_basedir, "limma_DEPs_withmeans_fixed_paired_samples.txt"), quote = F, sep = "\t", row.names = FALSE)



limma_res_with_means = read.delim(paste0(out_basedir, "limma_DEPs_withmeans_fixed_paired_samples.txt"), as.is = T, check.names = F, sep = "\t")


# stacked barplot
results <- limma_res_with_means %>%
  mutate(
    Regulation = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < 0 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

summary_data <- results %>%
  filter(Regulation != "Not significant") %>%  # Exclude non-significant proteins
  group_by(phase1_phase2, Regulation) %>%
  summarize(Count = n(), .groups = "drop") %>%
  mutate(phase1_phase2 = factor(phase1_phase2, 
                                levels = c("F_O", "F_EL", "F_LL", "O_EL", "O_LL", "EL_LL")))


barplot <- ggplot(summary_data, aes(x = phase1_phase2, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Phase Comparison",
    y = "Number of DEP",
    fill = "Effect direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors[c(3,2)])

# CSF3

csf3_plot <- ggplot(d_wide_adj_covar, aes(x = phase, y = CSF3, group = phase)) + 
  geom_boxplot(width = 0.3, color = my_colors[3], outliers = F) + 
  geom_jitter(alpha = 0.2, width = 0.1, color = my_colors[3]) + 
  theme_minimal() +
  ylab("CSF3 adjusted levels")

pdf(paste0(out_basedir, "plots/limma_DEPs_barplot_CSF3_fixed_limma.pdf"), width = 3, height = 8)
print((barplot + theme(legend.position="bottom"))/csf3_plot)
dev.off()


prot_subs <- unique(limma_res_with_means[limma_res_with_means$sign == T, "prot"])

# heatmap of mean values
mat_mean_abund_per_phase_adj_covar <- mean_abund_per_phase_adj_covar %>%
  filter(Protein %in% prot_subs) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

max_val <- max(abs(min(mat_mean_abund_per_phase_adj_covar[!row.names(mat_mean_abund_per_phase_adj_covar) %in% c('PROK1', 'CXCL13'),])), max(mat_mean_abund_per_phase_adj_covar[!row.names(mat_mean_abund_per_phase_adj_covar) %in% c('PROK1','CXCL13'),]))
mat_mean_abund_per_phase_adj_covar[mat_mean_abund_per_phase_adj_covar > max_val] <- max_val
breaksList = seq(-max_val, max_val, by = 0.01)
if(!0 %in% breaksList) breaksList <- sort(c(breaksList, 0))

full_palette <- rev(brewer.pal(n = 11, name = "RdYlBu"))
full_palette[ceiling(length(full_palette)/2)] <- "#FFFFFF"
colorList <- colorRampPalette(full_palette)(length(breaksList))

mat_mean_abund_per_phase_adj_covar <- mat_mean_abund_per_phase_adj_covar[,c("F", "O", "EL", "LL")]

signif_labels <- my_pivot_wider(limma_res_with_means[limma_res_with_means$prot %in% prot_subs,], "prot", "phase1_phase2", "adj.P.Val")
signif_labels <- ifelse(signif_labels < 0.05, "*", "")
signif_labels[is.na(signif_labels)] <- ""
signif_labels <- signif_labels[protein_order,]


split_at <- ceiling(length(protein_order) / 2)
protein_blocks <- list(
  protein_order[seq_len(split_at)],
  protein_order[seq.int(split_at + 1L, length(protein_order))]
)

comparison_order <- colnames(signif_labels)

plot_heatmap_significance <- function(clust_method = 'complete', num_k = 5){
  # Abundance panel
  hm_abund <- pheatmap(
    mat_mean_abund_per_phase_adj_covar[protein_order, , drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    clustering_method = clust_method,
    show_rownames = FALSE,
    color = abund_colors,
    breaks = breaksList,
    cellwidth = 12,
    cellheight = 6,
    angle_col = "90",
    silent = TRUE,
    cutree_rows = num_k
  )
  legend_scale <- 0.5  # 50% of the original height
  
  i <- which(hm_abund$gtable$layout$name == "legend")
  leg <- hm_abund$gtable$grobs[[i]]
  
  for (j in seq_along(leg$children)) {
    child <- leg$children[[j]]
    
    # Compress vertical positions, keeping the top fixed
    child$y <- grid::unit(1, "npc") -
      (grid::unit(1, "npc") - child$y) * legend_scale
    
    if (inherits(child, "rect")) {
      child$height <- child$height * legend_scale
    }
    
    leg$children[[j]] <- child
  }
  
  hm_abund$gtable$grobs[[i]] <- leg
  protein_order <- hm_abund$tree_row$labels[hm_abund$tree_row$order]
  clusters <- as.data.frame(cutree(hm_abund$tree_row,num_k)) %>%
    rownames_to_column("prot")
  colnames(clusters)[2] <- "cluster"
  row.names(clusters) <- clusters$prot
  clusters_ordered <- clusters[hm_abund$tree_row$labels[hm_abund$tree_row$order], ]
  # Positions where the cluster changes
  gaps_row <- which(clusters_ordered$cluster[-1] != clusters_ordered$cluster[-nrow(clusters_ordered)])
  
  # White cells, containing significance stars
  hm_signif <- pheatmap(
    matrix(
      0,
      nrow = length(protein_order),
      ncol = length(comparison_order),
      dimnames = list(protein_order, comparison_order)
    ),
    display_numbers = signif_labels[protein_order, comparison_order, drop = FALSE],
    number_color = "black",
    fontsize_number = 8,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    color = "white",
    breaks = c(-0.5, 0.5),
    border_color = "grey85",
    fontsize_row = 6,
    fontsize_col = 7,
    cellwidth = 6,
    cellheight = 6,
    legend = FALSE,
    silent = TRUE, 
    gaps_row = gaps_row
  )
  
  # Align the rows of the two panels exactly
  shared_heights <- unit.pmax(
    hm_abund$gtable$heights,
    hm_signif$gtable$heights
  )
  
  hm_abund$gtable$heights <- shared_heights
  hm_signif$gtable$heights <- shared_heights
  
  arrangeGrob(
    hm_abund$gtable,
    hm_signif$gtable,
    ncol = 2,
    widths = unit.c(
      sum(hm_abund$gtable$widths),
      sum(hm_signif$gtable$widths)
    )
  )
}


pdf(paste0(out_basedir, "plots/limma_DEPs_heatmap_fixed_limma_abundance_significance_split_ward_k5.pdf"),  width = 10,  height = 20)
combined <- plot_heatmap_significance(num_k = 5, clust_method = 'ward.D2')
grid::grid.newpage()
grid::grid.draw(combined)
dev.off()

