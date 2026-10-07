library(factoextra)
all_prots_traj <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_93_prots.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)

n_points = 100
fitted_matrix = all_prots_traj [,1:n_points]
#scaled_matrix <- t(apply(fitted_matrix, 1, scale))
#scaled_matrix <- fitted_matrix
scaled_matrix <- t(apply(fitted_matrix, 1, function(x) (x - min(x)) / (max(x) - min(x))))

colnames(scaled_matrix) <- seq(1,4, length.out=n_points)
# Correlation distance: (1 - Pearson Correlation)
dist_mat <- as.dist(1 - cor(t(scaled_matrix)))

hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 6)

#pam <- cluster::pam(dist_mat, diss = T, k = 5)
#clusters = pam$clustering

write.table(as.data.frame(clusters), file = paste0(out_basedir, "trajectories_gam/gam_clustering_inv_correl_dist.k6.txt"),quote = F, sep = "\t")


# Linear trajectories

all_prots_traj_lin <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_70_linear_prots.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)

n_points = 100
fitted_matrix_lin = all_prots_traj_lin [,1:n_points]
#scaled_matrix <- t(apply(fitted_matrix, 1, scale))
#scaled_matrix <- fitted_matrix
scaled_matrix_lin <- t(apply(fitted_matrix_lin, 1, function(x) (x - min(x)) / (max(x) - min(x))))

colnames(scaled_matrix_lin) <- seq(1,4, length.out=n_points)
# Correlation distance: (1 - Pearson Correlation)
dist_mat_lin <- as.dist(1 - cor(t(scaled_matrix_lin)))

hc_lin <- hclust(dist_mat_lin, method = "ward.D2")
clusters_lin <- cutree(hc_lin, k = 2)


write.table(as.data.frame(clusters_lin), file = paste0(out_basedir, "trajectories_gam/gam_linear_clustering_inv_correl_dist.k2.txt"),quote = F, sep = "\t")

clusters <- clusters + 2
clusters_combined <- c(clusters_lin, clusters)

fitted_matrix_combined <- rbind(fitted_matrix, fitted_matrix_lin)

write.table(as.data.frame(clusters_combined), file = paste0(out_basedir, "trajectories_gam/gam_linear_clustering_inv_correl_dist.combined.txt"),quote = F, sep = "\t")
write.table(as.data.frame(fitted_matrix_combined), file = paste0(out_basedir, "trajectories_gam/gam_linear_trajectories_163signif.txt"),quote = F, sep = "\t", row.names = TRUE, col.names = NA)

clusters_combined <- as.data.frame(clusters_combined) %>%
  rownames_to_column("protein")
fitted_matrix_combined <- as.data.frame(fitted_matrix_combined) %>%
  rownames_to_column("protein")

plot_data <-  left_join(fitted_matrix_combined, clusters_combined, by = "protein") %>%
  tidyr::pivot_longer(cols = -c(protein, clusters_combined), 
                      names_to = "time_point", 
                      values_to = "value") %>%
  mutate(time_point = as.numeric(gsub("X", "", time_point)))

pdf(paste0(out_basedir, "trajectories_gam/gam_clustering_inv_correl_dist.combined.pdf"))
# Plot trajectories by cluster
ggplot(plot_data, aes(x = time_point, y = value, group = prot,)) +
  geom_line(alpha = 0.3, color = 'dodgerblue3') +
  stat_summary(aes(group = clusters_combined), fun = mean, geom = "line", size = 1.5, color = 'dodgerblue4') +
  facet_wrap(~ clusters_combined) +
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
             FUNcluster = pam_cor, 
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



#####
#clustering of individual trajectories
library(factoextra)
set.seed(123) 

tmp_complete <- d_wide_adj_covar %>%
  filter(phase %in% all_phases) %>%
  group_by(ID) %>%
  filter(n_distinct(phase) == 4) %>%
  ungroup()
tmp_complete$phase <- factor(tmp_complete$phase, levels = all_phases)


prot_name = 'PROK1'
prot_name = 'REN'
prot_name = 'LEP'
prot_name = 'MMP7'
prot_name = 'ADA2'
p = cluster_individual_trajectories_per_prot(tmp_complete, prot_name, n_clusters = 3, all_prots_traj[prot_name,], center =T)
p$p_traj
p$p_clust


cluster_individual_trajectories_per_prot <- function(tmp_complete, prot, n_clusters = 3, gam_pred = NULL, center = T){
    
  wide_data <- tmp_complete %>%
    select(ID, phase, all_of(prot)) %>%
    pivot_wider(names_from = phase, values_from = all_of(prot)) %>%
    select(ID, all_of(all_phases))
  
  # Extract the numeric matrix (exclude ID column)
  mat <- as.matrix(wide_data[, -1])
  row.names(mat) <- wide_data$ID
  # center per ID
  if (center) {
    mat_centered <- t(apply(mat, 1, function(row) row - mean(row)))
  } else {
    mat_centered = mat
  }
  
  # Elbow + silhouette + gap statistic in one plot
  p1 <- fviz_nbclust(mat_centered, kmeans, method = "wss")   # elbow
  p2 <- fviz_nbclust(mat_centered, kmeans, method = "silhouette", print.summary = T)
  #fviz_nbclust(mat_centered, kmeans, method = "gap_stat", nboot = 50)
  
  nclusters = n_clusters
  km_centered <- kmeans(mat_centered, centers = nclusters, nstart = 25)
  mat_centered <- as.data.frame(mat_centered)
  mat_centered$cluster <- as.numeric(km_centered$cluster)
  
  plot_data <- mat_centered %>%
    rownames_to_column("ID") %>%
    select(ID, F, O, EL, LL, cluster) %>%
    pivot_longer(cols = c(F, O, EL,LL),
                 names_to = "phase",
                 values_to = "value") %>%
    mutate(phase = factor(phase, levels = all_phases), phase_num = as.numeric(phase))
  
  # Mean trajectory per cluster
  cluster_means <- plot_data %>%
    group_by(cluster, phase_num) %>%
    summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")
  
  if (!is.null(gam_pred)) {
    gam_pred = as.data.frame(t(gam_pred)) %>%
      rownames_to_column("phase_num") %>%
      mutate(phase_num = as.numeric(phase_num))
    colnames(gam_pred)[2] <- "GAM_curve"
    if (center) gam_pred$GAM_curve <- gam_pred$GAM_curve - mean(gam_pred$GAM_curve)
    p_traj <- ggplot(plot_data, aes(x = phase_num, y = value, group = ID, color = factor(cluster))) +
      geom_line(alpha = 0.4) +
      geom_line(data = cluster_means, aes(x = phase_num, y = mean, group = cluster, color = factor(cluster)),
                linewidth = 1, inherit.aes = FALSE) +
      geom_line(data = gam_pred, aes(x = phase_num, y = GAM_curve),
                linewidth = 1, inherit.aes = FALSE) +
      labs(title = paste0(prot, " trajectories by cluster"),
           color = "Cluster") +
      ylab("adjusted centered abundance") +
      theme_minimal() +
      scale_x_continuous(
        breaks = c(1, 2, 3, 4),
        labels = c("F", "O", "EL", "LL")
      )
  } else {
  p_traj <- ggplot(plot_data, aes(x = phase, y = value, group = ID, color = factor(cluster))) +
    geom_line(alpha = 0.4) +
    geom_line(data = cluster_means, aes(x = phase, y = mean, group = cluster, color = factor(cluster)),
               linewidth = 1, inherit.aes = FALSE) +
    labs(title = paste0(prot, " trajectories by cluster"),
         color = "Cluster") +
    ylab("adjusted centered abundance") +
    theme_minimal()
  }
  return(list(p_traj = p_traj, p_clust = (p1 + p2)))
}







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
