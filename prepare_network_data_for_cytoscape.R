library(igraph)
library(dplyr)
library(readr)


setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12/results12/intensity_shared_prots_261125/network")
edges_data <- read.delim("network.spline.edges.causality2.with_pheno-pheno.txt", sep = "\t", check.names = F, as.is = T)  
nodes_data <- read.delim("network.spline.nodes.txt", sep = "\t", check.names = F, as.is = T) 
annot <- read.delim("Su_MR.subset_cut.txt", sep = "\t", check.names = F, as.is = T)

# A. Process edges
edges_processed <- edges_data %>%
  # Join the data with itself, swapping prot and pheno to find reverse matches
  left_join(
    edges_data %>% select(prot, pheno, has_direction),
    by = c("prot" = "pheno", "pheno" = "prot"),
    suffix = c("", "_rev")
  ) %>%
  mutate(
    # Create the new column: TRUE if a reverse edge exists AND it has direction
    reverse_direction = !is.na(has_direction_rev) & has_direction_rev == TRUE
  ) %>%
  filter(!reverse_direction | (reverse_direction & prot < pheno)) %>%
  select(-has_direction_rev)


edges_prep <- edges_processed %>%
  mutate(
    from = prot,
    to = pheno,
    strong_assoc = abs(estimate) > 0.15,
    color = case_when(
      estimate > 0 ~ "#FF0000", # Red
      estimate < 0 ~ "#3498DB"  # Blue
    ),
    lty = 1, 
    arrow.mode = ifelse(has_direction == TRUE, 2, 0),  # 2 = arrow to 'to', 0 = no arrow
    rev_arrow.mode = ifelse(reverse_direction == TRUE, 2, 0)
  ) %>%
  select(-prot, -pheno,  -estimate, -has_direction, -reverse_direction)


# B. Process Nodes
nodes_prep <- nodes_data %>%
  filter(feature %in% c(edges_prep$from, edges_prep$to)) %>%
  mutate(
    id = feature,
    color = case_when(
      type == "protein" ~ "#97C2FC",
      type == "phenotype" ~ "#FFB6C1",
      type == "hormone" ~ "#98FB98",
      TRUE ~ "#DDDDDD"
    ),
    shape = "circle", 
    size = 6,
    label.color = "black"
  )



# C. Add Disease Logic 
SHOW_DISEASES <- TRUE 

if (SHOW_DISEASES) {
  visible_proteins <- nodes_prep$id
  annot_filtered <- annot %>% filter(Assay %in% visible_proteins)
  
  if (nrow(annot_filtered) > 0) {
    disease_edges <- data.frame(
      from = annot_filtered$Assay,
      to = annot_filtered$MR_outcomes,
      color = ifelse(annot_filtered$Beta > 0, "#FF0000", "#3498DB"),
      lty = 2, # Dashed line for diseases
      arrow.mode = 2, 
      rev_arrow.mode = 0
    )
    
    unique_diseases <- unique(annot_filtered$MR_outcomes)
    disease_nodes <- data.frame(
      id = unique_diseases,
      feature = unique_diseases,
      type = "disease",
      color = "#FFA500", # Orange
      shape = "square",  # Box shape
      size = 8,
      label.color = "black",
      nodes_to_select = TRUE
    )
    
    # Combine
    edges_prep <- bind_rows(edges_prep, disease_edges)
    nodes_prep <- bind_rows(nodes_prep, disease_nodes)
  }
}

edges_prep <- unique(edges_prep)

write.table(nodes_prep, file = "nodes_all.txt", quote = F, sep = "\t", row.names = F)
write.table(edges_prep, file = "edges_all.txt", quote = F, sep = "\t", row.names = F)

# strong
edges_prep_strong <- edges_prep[edges_prep$strong_assoc == T,]
write.table(edges_prep_strong, file = "edges_strong.txt", quote = F, sep = "\t", row.names = F)

# no leafs, no diseases
node_counts <- table(c(edges_prep$from, edges_prep$to))
node_sel <- names(node_counts[node_counts > 1])
node_sel <- node_sel[! node_sel %in% nodes_prep[nodes_prep$type == 'disease', "feature"]]
edges_prep_multi <- edges_prep[edges_prep$from %in% node_sel & edges_prep$to %in% node_sel,]
write.table(edges_prep_multi, file = "edges_mulitple.txt", quote = F, sep = "\t", row.names = F)


# strong & no leafs
node_counts <- table(c(edges_prep_strong$from, edges_prep_strong$to))
node_sel <- names(node_counts[node_counts > 1])

edges_prep_strong_multi <- edges_prep_strong[edges_prep_strong$from %in% node_sel & edges_prep_strong$to %in% node_sel,]
write.table(edges_prep_strong_multi, file = "edges_strong_mulitple.txt", quote = F, sep = "\t", row.names = F)



### shared between hormones and phenotypes

edges_with_types <- edges_prep %>%
  left_join(nodes_data %>% select(feature, type_from = type), by = c("from" = "feature")) %>%
  left_join(nodes_data %>% select(feature, type_to = type), by = c("to" = "feature")) %>%
  mutate(
    edge_type_cat = paste(pmin(type_from, type_to), "-", pmax(type_from, type_to))
  ) 


prots_h <- c(edges_with_types[edges_with_types$type_from == 'protein' & edges_with_types$type_to == 'hormone', "from"], edges_with_types[edges_with_types$type_to == 'protein' & edges_with_types$type_from == 'hormone', "to"])
prots_ph <- c(edges_with_types[edges_with_types$type_from == 'protein' & edges_with_types$type_to != 'hormone', "from"], edges_with_types[edges_with_types$type_to == 'protein' & edges_with_types$type_from != 'hormone', "to"])
shared_prots <- intersect(prots_h, prots_ph)

edges_prep <- edges_prep[edges_prep$from %in% shared_prots | edges_prep$to %in% shared_prots,]
nodes_prep <- nodes_prep[nodes_prep$feature %in% edges_prep$to | nodes_prep$feature %in% edges_prep$from,]

edges_prep[,c("prot", "pheno", "estimate", "strong_assoc", "has_direction")] <- NULL

visible_proteins <- nodes_prep$id
annot_filtered <- annot %>% filter(Assay %in% visible_proteins)

if (nrow(annot_filtered) > 0) {
  disease_edges <- data.frame(
    from = annot_filtered$Assay,
    to = annot_filtered$MR_outcomes,
    color = ifelse(annot_filtered$Beta > 0, "#FF0000", "#3498DB"),
    lty = 2, # Dashed line for diseases
    width = 1,
    arrow_size = 0.2,
    arrow.mode = 2
  )
  
  unique_diseases <- unique(annot_filtered$MR_outcomes)
  disease_nodes <- data.frame(
    feature = unique_diseases,
    type = "disease",
    id = unique_diseases,
    color = "#FFA500", # Orange
    shape = "square",  # Box shape
    size = 8,
    label.color = "black"
  )
  
  # Combine
  edges_prep <- bind_rows(edges_prep, disease_edges)
  nodes_prep <- bind_rows(nodes_prep, disease_nodes)
}


write.table(nodes_prep, file = "nodes_shared.txt", quote = F, sep = "\t", row.names = F)
write.table(unique(edges_prep), file = "edges_shared.txt", quote = F, sep = "\t", row.names = F)

