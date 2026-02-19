library(factoextra)
all_prots_traj <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_93_prots.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)
all_prots_traj <- all_prots_traj[,round(seq(1,100,length.out = 20))]

all_prots_traj_lin <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_70_linear_prots.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)
all_prots_traj_lin <- all_prots_traj_lin[,round(seq(1,100,length.out = 20))]

all_prots_traj_b2 <- read.delim(paste0(out_basedir, "../intensity_batch2_prots_261125/trajectories_gam/protein_trajectories_gam_batch2.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)
signif_b2 <- c("NOS3","AMDHD1","CTSV","KLB","KLK4","CXCL13","FGFR4","NELL1","OXT","WFDC2","FCRL5","NPTX2","STC2","FCN1","AMOT","ADAMTS15","MSMP","B3GNT7","PAMR1","FGF21","CTSF","FOLR2","CYTL1","CD8A","C19orf12","SFRP4","SERPINA9","TREH","INHBB","ABI3BP","VNN2","PRL","FOLR3","ITIH4","ECHDC3","VWC2","INHBC","HS3ST3B1","NPY","HAO1","GCNT1","CILP","MYL3","ASGR1","TGFBR2","LGALS3BP","HSPA13","GAST","CRELD1","PIK3IP1","CLEC7A","CES2","TNFRSF10A","CKB","SH3D19","SHBG","KRT18","EDA2R","FGF23","KLK1","MAMDC4","PRTG","TNFRSF6B","SORD","DCXR","ST3GAL1","CPB2","MATN3","PRSS53","CCN4","KIRREL2","DEFB1","NCR3LG1","KHK","IL7R","AMBP","CRISP3","THY1","SOD2","WFIKKN1","GOLM2","PVR","TMSB10","CD86","TNFAIP6","CD248","IFNAR2")
all_prots_traj_b2 <- all_prots_traj_b2[signif_b2,1:20]
colnames(all_prots_traj_b2) <- colnames(all_prots_traj)

b2_linear <- c("NOS3","AMDHD1","KLB","FCN1","PAMR1","CTSF","CYTL1","CD8A","SERPINA9","TREH","ABI3BP","VNN2","PRL","FOLR3","ITIH4","VWC2","HAO1","TGFBR2","GAST","PIK3IP1","SH3D19","SHBG","KLK1","PRTG","TNFRSF6B","DCXR","ST3GAL1","CPB2","KIRREL2","DEFB1","IL7R","AMBP","CRISP3","WFIKKN1","GOLM2","TNFAIP6","FOLR1","MDGA1","IGF2R","PRG3")
all_prots_traj <- rbind(all_prots_traj, all_prots_traj_lin, all_prots_traj_b2)

n_points = 20

# linear
fitted_matrix_lin = all_prots_traj[row.names(all_prots_traj) %in% c(b2_linear, row.names(all_prots_traj_lin)),1:n_points]
#scaled_matrix <- t(apply(fitted_matrix, 1, scale))
#scaled_matrix <- fitted_matrix
scaled_matrix_lin <- t(apply(fitted_matrix_lin, 1, function(x) (x - min(x)) / (max(x) - min(x))))

colnames(scaled_matrix_lin) <- colnames(all_prots_traj)
# Correlation distance: (1 - Pearson Correlation)
dist_mat_lin <- as.dist(1 - cor(t(scaled_matrix_lin)))

hc_lin <- hclust(dist_mat_lin, method = "ward.D2")
clusters_lin <- cutree(hc_lin, k = 2)

# non-linear

fitted_matrix = all_prots_traj[! row.names(all_prots_traj) %in% c(b2_linear, row.names(all_prots_traj_lin)),1:n_points]
#scaled_matrix <- t(apply(fitted_matrix, 1, scale))
#scaled_matrix <- fitted_matrix
scaled_matrix <- t(apply(fitted_matrix, 1, function(x) (x - min(x)) / (max(x) - min(x))))

colnames(scaled_matrix) <- colnames(all_prots_traj)
# Correlation distance: (1 - Pearson Correlation)
dist_mat <- as.dist(1 - cor(t(scaled_matrix)))

hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 6)

clusters <- clusters + 2
clusters_combined <- c(clusters_lin, clusters)
fitted_matrix_combined <- rbind(fitted_matrix, fitted_matrix_lin)

clusters_combined <- as.data.frame(clusters_combined) %>%
  rownames_to_column("protein")
fitted_matrix_combined <- as.data.frame(fitted_matrix_combined) %>%
  rownames_to_column("protein")

colnames(clusters_combined)[2] <- "clusters"

write.table(as.data.frame(clusters_combined), file = paste0(out_basedir, "trajectories_gam/gam_linear_clustering_inv_correl_dist.combined.b1+b2.txt"),quote = F, sep = "\t")
write.table(as.data.frame(fitted_matrix_combined), file = paste0(out_basedir, "trajectories_gam/gam_trajectories.b1+b2.txt"),quote = F, sep = "\t", row.names = TRUE, col.names = NA)

plot_data <-  left_join(fitted_matrix_combined, clusters_combined, by = "protein") %>%
  tidyr::pivot_longer(cols = -c(protein, clusters), 
                      names_to = "time_point", 
                      values_to = "value") %>%
  mutate(time_point = as.numeric(gsub("X", "", time_point)))

pdf(paste0(out_basedir, "trajectories_gam/gam_clustering_inv_correl_dist.combined.b1+b2.pdf"))
# Plot trajectories by cluster
ggplot(plot_data, aes(x = time_point, y = value, group = protein,)) +
  geom_line(alpha = 0.3, color = 'dodgerblue3') +
  stat_summary(aes(group = clusters), fun = mean, geom = "line", size = 1.5, color = 'dodgerblue4') +
  facet_wrap(~ clusters) +
  labs(x = "Phase", y = "Scaled GAM fitted protein trajectories") +
  theme_minimal()


hclust_cor <- function(x, k) {
  dist_mat <- as.dist(1 - cor(t(x))) 
  hc <- hclust(dist_mat, method = "ward.D2")
  list(data = x, cluster = cutree(hc, k = k))
}

pam_cor <- function(x, k) {
  dist_mat <- as.dist(1 - cor(t(x))) 
  pam_res <- cluster::pam(dist_mat, diss = T, k = k)
  list(cluster = pam_res$clustering)
}



# 2. Run Silhouette method using this custom function
set.seed(123)
p1 <- fviz_nbclust(scaled_matrix, 
                   FUNcluster = hclust_cor, 
                   method = "silhouette")

p2 <- fviz_nbclust(scaled_matrix, 
                   FUNcluster = hclust_cor, 
                   method = "wss")


# 2. Visualize
p3 <- fviz_cluster(hclust_cor(scaled_matrix, 6), 
                   geom = "point",
                   ellipse.type = "convex", 
                   ggtheme = theme_minimal(),
                   main = "Protein Trajectory Clusters (PCA projection)")

print(p1+p2)
print(p3)
dev.off()





# cluster linear and non-linear together
all_prots_traj <- read.delim(paste0(out_basedir, "trajectories_gam/gam_linear_trajectories_163signif.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)

n_points = 100
fitted_matrix = all_prots_traj [,1:n_points]
#scaled_matrix <- t(apply(fitted_matrix, 1, scale))
#scaled_matrix <- fitted_matrix
scaled_matrix <- t(apply(fitted_matrix, 1, function(x) (x - min(x)) / (max(x) - min(x))))

colnames(scaled_matrix) <- seq(1,4, length.out=n_points)
# Correlation distance: (1 - Pearson Correlation)
dist_mat <- as.dist(1 - cor(t(scaled_matrix)))

hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 5)

#pam <- cluster::pam(dist_mat, diss = T, k = 5)
#clusters = pam$clustering


plot_data <- as.data.frame(fitted_matrix) %>%
  mutate(protein = rownames(fitted_matrix),
         cluster = as.factor(clusters)) %>%
  tidyr::pivot_longer(cols = -c(protein, cluster),
                      names_to = "time_point",
                      values_to = "value") %>%
  mutate(time_point = as.numeric(gsub("X", "", time_point)))

ggplot(plot_data, aes(x = time_point, y = value, group = protein)) +
  geom_line(alpha = 0.3, color = 'dodgerblue3') +
  stat_summary(aes(group = cluster), fun = mean, geom = "line", size = 1.5, color = 'dodgerblue4') +
  facet_wrap(~ cluster) +
  labs(title = "Hierarchical clustering (Ward D2) based on inverse correlation distance",
       x = "Phase", y = "Scaled GAM fitted protein trajectories") +
  theme_minimal()




# combine with limma signif results


limma_res_all$cluster <- clusters[limma_res_all$prot]
phase_levels <- c("F_O", "O_EL", "EL_LL") 
midpoints    <- c(1.5, 2.5, 3.5) 
names(midpoints) <- phase_levels

stats_data <- limma_res_all %>%
  filter(signif_direction != 0) %>%
  filter(phase1_phase2 %in% phase_levels) %>%
  group_by(cluster, phase1_phase2, signif_direction) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(x_pos = midpoints[phase1_phase2]) %>%
  mutate(plot_count = ifelse(signif_direction == -1, -count, count)) %>%
  mutate(Direction = factor(ifelse(signif_direction == 1, "Up", "Down"), levels = c("Up", "Down")))

p1 <- ggplot(plot_data, aes(x = time_point, y = value, group = protein)) +
  geom_line(alpha = 0.3, color = 'dodgerblue3') +
  stat_summary(aes(group = cluster), fun = mean, geom = "line", size = 1.5, color = 'dodgerblue4') +
  facet_wrap(~cluster, ncol = 3) +                 
  theme_minimal() +
  labs(title = "Protein Trajectories per Cluster", x = NULL, y = "Scaled Intensity") +
  theme(axis.text.x = element_blank())   

p2 <- ggplot(stats_data, aes(x = x_pos, y = plot_count, fill = Direction)) +
  geom_col(width = 0.6) + # width controls bar thickness
  facet_wrap(~cluster, ncol = 3) +
  scale_fill_manual(values = c("Up" = my_colors[2], "Down" = my_colors[3])) +
  geom_hline(yintercept = 0, color = "black", size = 0.2) + # Zero line
  scale_x_continuous(breaks = 1:4, labels = c("F", "O", "EL", "LL"), limits = c(1, 4)) +
  labs(y = "Count (Sig)", x = "Phases") +
  theme_minimal() +
  theme(
    strip.text = element_blank(), # Hide facet titles for the bottom plot (redundant)
    legend.position = "bottom"
  )

combined_plot <- p1 / p2 + 
  plot_layout(heights = c(3, 1)) # Trajectories are 3x taller than bars

print(combined_plot)




plot_data <- as.data.frame(fitted_matrix) %>%
  mutate(protein = rownames(fitted_matrix),
         cluster = as.factor(clusters)) %>%
  tidyr::pivot_longer(cols = -c(protein, cluster), 
                      names_to = "time_point", 
                      values_to = "value") %>%
  mutate(time_point = as.numeric(gsub("X", "", time_point)))

pdf(paste0(out_basedir, "trajectories_gam/gam_linear_clustering_inv_correl_dist.k6.pdf"), height = 4, width = 4)
# Plot trajectories by cluster
ggplot(plot_data, aes(x = time_point, y = value, group = protein,)) +
  geom_line(alpha = 0.3, color = 'dodgerblue3') +
  stat_summary(aes(group = cluster), fun = mean, geom = "line", size = 1.5, color = 'dodgerblue4') +
  facet_wrap(~ cluster) +
  labs(title = "Hierarchical clustering (Ward D2) based on inverse correlation distance",
       x = "Phase", y = "Scaled GAM fitted protein trajectories") +
  theme_minimal()

dev.off()










### test clustering on logFC
# limma_res_wide_tmp <- my_pivot_wider(limma_res_all[limma_res_all$prot != 'PROK1',], "prot", "phase1_phase2", "logFC")
# limma_res_wide_tmp <- limma_res_wide_tmp[,c("F_O", "O_EL", "EL_LL")]
# 
# k <- 6 # Choose number of clusters
# 
# #hc <- hclust(as.dist(dist_matrix), method = "ward.D2")
# #clusters <- cutree(hc, k = k)
# 
# kmeans_result <- kmeans(limma_res_wide_tmp, centers = k)
# #kmeans_result <- cluster::pam(limma_res_wide_tmp, k, diss = F)
# clusters <- as.data.frame(kmeans_result$cluster) %>%
#   rownames_to_column("prot")
# colnames(clusters)[2] <- "cluster"
# 
# 
# median_by_phase_long <- d_wide %>%
#   mutate(across(where(is.numeric), scale)) %>%
#   pivot_longer(cols = where(is.numeric), names_to = "prot", values_to = "value") %>%
#   group_by(phase, prot) %>%
#   summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")
# 
# plot_data <- inner_join(median_by_phase_long, clusters, by = "prot")
# 
# # Plot trajectories by cluster
# ggplot(plot_data, aes(x = phase, y = median_value, group = prot)) +
#   geom_line(alpha = 0.2, color = 'dodgerblue4') +
#   facet_wrap(~ cluster) +
#   labs(title = "K-means clustering based on limma DE ",
#        x = "Phase", y = "Mean protein level") +
#   stat_summary(aes(group = cluster), fun = mean, geom = "line", size = 1, color = 'black') +
#   theme_minimal()
# 
