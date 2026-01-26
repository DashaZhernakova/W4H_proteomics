library(shiny)
library(visNetwork)
library(dplyr)
library(readr)
library(bslib)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12/results12/intensity_shared_prots_261125/network")
edges_data <- read.delim("network.spline.edges.causality2.with_pheno-pheno.txt", sep = "\t", check.names = F, as.is = T)  
nodes_data <- read.delim("network.spline.nodes.txt", sep = "\t", check.names = F, as.is = T) 
annot <- read.delim("Su_MR.subset_cut.txt", sep = "\t", check.names = F, as.is = T)

node_counts <- table(c(edges_data$prot, edges_data$pheno))
node_sel <- names(node_counts[node_counts > 1])
nodes_data$nodes_to_select <- ifelse(nodes_data$feature %in% node_sel, TRUE, FALSE)

# Process edges data
edges_vis <- edges_data %>%
  mutate(
    from = prot,
    to = pheno,
    strong_assoc = if_else(abs(estimate) > 0.15, TRUE, FALSE),
    title = paste0(prot, " → ", pheno, 
                   "<br>Estimate: ", round(estimate, 3),
                   "<br>Directed: ", has_direction),
    color = case_when(
      estimate > 0 ~ "#FF0000",
      estimate < 0 ~ "#3498DB",
    ),
    arrows = ifelse(has_direction, "to", ""),
    width = 1 ,
    smooth = F,
  )

# Process nodes data
nodes_vis <- nodes_data %>%
  filter(feature %in% c(edges_vis$from, edges_vis$to)) %>%
  mutate(
    id = feature,
    label = feature,  
    title = paste0("Type: ", type, "<br>ID: ", feature),
    color = case_when(
      type == "protein" ~ "#97C2FC",
      type == "phenotype" ~ "#FFB6C1",
      type == "hormone" ~ "#98FB98"
    ),
    node_with_two_edges = nodes_to_select,
    shape = 'circle',
    size = 25,
    borderWidth = 2,
    borderWidthSelected = 4
  )

edges_with_types <- edges_vis %>%
  left_join(nodes_data %>% select(feature, type_from = type), by = c("from" = "feature")) %>%
  left_join(nodes_data %>% select(feature, type_to = type), by = c("to" = "feature")) %>%
  mutate(
    edge_type_cat = paste(pmin(type_from, type_to), "-", pmax(type_from, type_to))
  )

# MR disease annot
annot_flt <- annot %>% filter(Assay %in% nodes_vis$feature)

disease_edges_all <- data.frame(
  from = annot_flt$Assay,
  to = annot_flt$MR_outcomes,
  color = case_when(
    annot_flt$Beta > 0 ~ "#FF0000",
    annot_flt$Beta < 0 ~ "#3498DB",
  ),      
  dashes = TRUE,
  width = 1.5,
  arrows = 'to',
  smooth = FALSE,
  title = paste("Protein:", annot_flt$Assay, "<br>Disease:", annot_flt$MR_outcomes),
  edge_type_cat = "protein - disease",
  strong_assoc = FALSE
)

unique_diseases_all <- unique(annot_flt$MR_outcomes)

disease_nodes_all <- data.frame(
  id = unique_diseases_all,
  feature = unique_diseases_all,
  label = unique_diseases_all,
  title = paste("Disease (MR Outcome):", unique_diseases_all),
  type = "disease",
  color = "#FFA500",
  shape = "box",
  size = 40,
  borderWidth = 2,
  borderWidthSelected = 4,
  widthConstraint = 150,
  node_with_two_edges = NA
)

ui <- page_sidebar(
  
  sidebar = sidebar(
    checkboxGroupInput("edge_filters", 
                       "Show Edge Types:",
                       choices = c("hormone - hormone",
                                   "phenotype - phenotype",
                                   "hormone - phenotype", 
                                   "hormone - protein", 
                                   "phenotype - protein"),
                       selected = c("hormone - hormone",
                                    "phenotype - phenotype",
                                    "hormone - phenotype", 
                                    "hormone - protein", 
                                    "phenotype - protein")),
    hr(),
    checkboxInput(
      inputId = "filter_strong", 
      label = "Show only associations with |estimate| > 0.15", 
      value = TRUE
    ),
    hr(),
    checkboxInput(
      inputId = "filter_two_edges", 
      label = "Show Nodes with Multiple Connections", 
      value = TRUE
    ),
    hr(),
    checkboxInput(
      inputId = "show_diseases",
      label = "Show Associated Diseases (MR)",
      value = FALSE
    ),
    hr(),
    helpText("Click a node to highlight its specific connections. Click empty space to reset.")
  ),
  
  visNetworkOutput("network_plot", height = "800px")
)

server <- function(input, output, session) {
  
  # --- 1. Reactive Data ---
  graph_data <- reactive({
    
    filtered_edges <- edges_with_types %>%
      filter(edge_type_cat %in% input$edge_filters)
    
    if (isTRUE(input$filter_strong)) {
      filtered_edges <- filtered_edges %>%
        filter(strong_assoc == TRUE)
    }
    
    active_node_ids <- unique(c(filtered_edges$from, filtered_edges$to))
    filtered_nodes <- nodes_vis %>% filter(id %in% active_node_ids)
    
    if (isTRUE(input$filter_two_edges)) {
      filtered_nodes <- filtered_nodes %>%
        filter(node_with_two_edges == TRUE)
      
      current_degrees <- table(c(filtered_edges$from, filtered_edges$to))
      nodes_with_2_plus_edges <- names(current_degrees[current_degrees >= 2])
      
      filtered_nodes <- filtered_nodes %>%
        filter(id %in% nodes_with_2_plus_edges)
      
      filtered_edges <- filtered_edges %>%
        filter(from %in% filtered_nodes$id & to %in% filtered_nodes$id)
      
      final_active_ids <- unique(c(filtered_edges$from, filtered_edges$to))
      filtered_nodes <- filtered_nodes %>% filter(id %in% final_active_ids)
    }
    
    if (isTRUE(input$show_diseases)) {
      visible_proteins <- filtered_nodes$id
      disease_edges <- disease_edges_all %>% filter(from %in% visible_proteins)
      
      if (nrow(annot_filtered) > 0) {
        disease_nodes <- disease_nodes_all %>% filter(id %in% disease_edges$to)
        
        filtered_edges <- bind_rows(filtered_edges, disease_edges)
        filtered_nodes <- bind_rows(filtered_nodes, disease_nodes)
      }
    }
    
    # Add Edge IDs for Proxy
    filtered_edges <- filtered_edges %>% mutate(id = paste0("e", row_number()))
    
    filtered_nodes <- filtered_nodes %>% arrange(id)
    
    list(nodes = filtered_nodes, edges = filtered_edges)
  })
  
  # --- 2. Render Network ---
  output$network_plot <- renderVisNetwork({
    dat <- graph_data()
    
    visNetwork(nodes = dat$nodes, edges = dat$edges) %>%
      visNodes(
        shape = "circle",
        widthConstraint = list(minimum = 75, maximum = 75), 
        heightConstraint = list(minimum = 75, valign = "middle"),
        font = list(size = 16, color = "#000000", face = "arial", strokeWidth = 0),
        opacity = 1,
        scaling = list(
          label = list(enabled = TRUE, min = 20, max = 50, maxVisible = 10000)
        )
      ) %>%
      visOptions(
        highlightNearest = FALSE, 
        nodesIdSelection = list(
          enabled = TRUE,
          values = sort(unique(dat$nodes$id))
        )
      ) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -200, 
          centralGravity = 0.005,       
          springLength = 250,           
          springConstant = 0.05         
        ),
        stabilization = list(iterations = 150)
      ) %>%
      visEvents(
        stabilizationIterationsDone = "function() { this.setOptions({physics: false}); }"
      ) %>%
      visInteraction(
        hideNodesOnDrag = FALSE, hideEdgesOnDrag = FALSE
      )
  })
  
  # --- 3. Custom Selection Logic ---
  observeEvent(input$network_plot_selected, {
    
    selected_id <- input$network_plot_selected
    dat <- graph_data()
    current_nodes <- dat$nodes
    current_edges <- dat$edges
    
    proxy <- visNetworkProxy("network_plot")
    
    if (is.null(selected_id) || selected_id == "") {
      
      # --- RESET LOGIC ---
      # We must explicitly set opacity = 1 because the highlight logic reduced it
      nodes_reset <- current_nodes %>% mutate(opacity = 1)
      edges_reset <- current_edges # Original colors
      
      visUpdateNodes(proxy, nodes = nodes_reset)
      visUpdateEdges(proxy, edges = edges_reset)
      
    } else {
      
      # --- HIGHLIGHT LOGIC ---
      
      # 1. Identify neighbors
      incident_edges <- current_edges %>% 
        filter(from == selected_id | to == selected_id)
      
      neighbor_ids <- unique(c(incident_edges$from, incident_edges$to))
      nodes_to_keep_ids <- unique(c(selected_id, neighbor_ids))
      edges_to_keep_ids <- incident_edges$id
      
      # 2. Dim others (Nodes)
      nodes_update <- current_nodes %>%
        mutate(
          color = ifelse(id %in% nodes_to_keep_ids, color, "rgba(200,200,200,0.3)"),
          opacity = ifelse(id %in% nodes_to_keep_ids, 1, 0.3)
        )
      
      # 3. Dim others (Edges)
      edges_update <- current_edges %>%
        mutate(
          color = ifelse(id %in% edges_to_keep_ids, color, "rgba(220,220,220,0.1)"),
          width = ifelse(id %in% edges_to_keep_ids, width, 1)
        )
      
      visUpdateNodes(proxy, nodes = nodes_update)
      visUpdateEdges(proxy, edges = edges_update)
    }
  })
}

shinyApp(ui = ui, server = server)
