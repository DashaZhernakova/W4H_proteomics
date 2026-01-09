library(shiny)
library(visNetwork)
library(dplyr)
library(readr)

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

# Process nodes data WITH CORRECTED FONT SETTINGS
nodes_vis <- nodes_data %>%
  # 1) FILTER: Keep only nodes that appear in edges
  filter(feature %in% c(edges_vis$from, edges_vis$to)) %>%
  mutate(
    id = feature,
    label = feature,  # This makes the name visible on the node
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
    # Font settings will be applied separately (see below)
  )

edges_with_types <- edges_vis %>%
  left_join(nodes_data %>% select(feature, type_from = type), by = c("from" = "feature")) %>%
  left_join(nodes_data %>% select(feature, type_to = type), by = c("to" = "feature")) %>%
  mutate(
    edge_type_cat = paste(pmin(type_from, type_to), "-", pmax(type_from, type_to))
  )


ui <- fluidPage(
  titlePanel("Interactive Network"),
  sidebarLayout(
    sidebarPanel(
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
      helpText("Disease nodes will appear as Orange Boxes connected to Proteins."),
      hr(),
      helpText("The network will automatically hide edges that point to filtered-out nodes.")
    ),
    mainPanel(
      visNetworkOutput("network_plot", height = "800px")
    )
  )
)

server <- function(input, output) {
  
  output$network_plot <- renderVisNetwork({
    
    # --- 1. Filter edges based on Category ---
    filtered_edges <- edges_with_types %>%
      filter(edge_type_cat %in% input$edge_filters)
    
    # --- 2. Filter edges based on Strong Association ---
    if (isTRUE(input$filter_strong)) {
      filtered_edges <- filtered_edges %>%
        filter(strong_assoc == TRUE)
    }
    
    # --- 3. Initial Node Selection ---
    # Select nodes that are part of the currently filtered edges
    active_node_ids <- unique(c(filtered_edges$from, filtered_edges$to))
    filtered_nodes <- nodes_vis %>% filter(id %in% active_node_ids)
    
    # --- 4. Filter Nodes based on 'Two-Edge' Checkbox ---
    if (isTRUE(input$filter_two_edges)) {
      
      # Step A: Filter by the static column
      filtered_nodes <- filtered_nodes %>%
        filter(node_with_two_edges == TRUE)
      
      # Step B: Filter by VISIBLE degree
      # We calculate how many edges each node has in the CURRENT filtered_edges
      # and remove nodes that have fewer than 2 edges visible.
      current_degrees <- table(c(filtered_edges$from, filtered_edges$to))
      nodes_with_2_plus_edges <- names(current_degrees[current_degrees >= 2])
      
      filtered_nodes <- filtered_nodes %>%
        filter(id %in% nodes_with_2_plus_edges)
      
      # Step C: Re-sync Edges
      # Remove edges connected to nodes we just filtered out
      filtered_edges <- filtered_edges %>%
        filter(from %in% filtered_nodes$id & to %in% filtered_nodes$id)
      
      # Step D: Final Cleanup (Optional but recommended)
      # If removing edges in Step C left some NEW nodes floating (0 edges), remove them.
      final_active_ids <- unique(c(filtered_edges$from, filtered_edges$to))
      filtered_nodes <- filtered_nodes %>% filter(id %in% final_active_ids)
    }
    if (isTRUE(input$show_diseases)) {
      
      # 1. Identify currently visible proteins
      # We assume the 'annot' Assay column matches the node IDs
      visible_proteins <- filtered_nodes$id
      
      # 2. Filter 'annot' for these proteins
      annot_filtered <- annot %>% 
        filter(Assay %in% visible_proteins)
      
      if (nrow(annot_filtered) > 0) {
        # 3. Create Disease Edges 
        disease_edges <- data.frame(
          from = annot_filtered$Assay,
          to = annot_filtered$MR_outcomes,
          color = case_when(
            annot_filtered$Beta > 0 ~ "#FF0000",
            annot_filtered$Beta < 0 ~ "#3498DB",
          ),      
          dashes = TRUE,          # Dashed lines to distinguish from causal network
          width = 1.5,
          arrow = 'to',
          smooth = FALSE,
          title = paste("Protein:", annot_filtered$Assay, "<br>Disease:", annot_filtered$MR_outcomes),
          edge_type_cat = "protein - disease" # Helper column
        )
        
        # 4. Create Disease Nodes
        # We need unique diseases from the filtered annotation
        unique_diseases <- unique(annot_filtered$MR_outcomes)
        
        disease_nodes <- data.frame(
          id = unique_diseases,
          feature = unique_diseases,
          label = unique_diseases,
          title = paste("Disease (MR Outcome):", unique_diseases),
          type = "disease",
          color = "#FFA500",      # Orange color
          shape = "box",          # Box shape to distinguish from molecules
          size = 75,
          borderWidth = 2,
          borderWidthSelected = 4,
          widthConstraint = 150,
          node_with_two_edges = FALSE # Logic flag, not strictly used for display
        )
        
        # 5. Merge with existing Network Data
        # bind_rows will fill missing columns with NA, which visNetwork ignores
        filtered_edges <- bind_rows(filtered_edges, disease_edges)
        filtered_nodes <- bind_rows(filtered_nodes, disease_nodes)
      }
    }
    
    # --- 5. Render Network ---
    visNetwork(nodes = filtered_nodes, edges = filtered_edges) %>%
      visNodes(
        shape = "circle",
        widthConstraint = list(minimum = 75, maximum = 75), 
        heightConstraint = list(minimum = 75, valign = "middle"),
        font = list(size = 16, color = "#000000", face = "arial", strokeWidth = 0),
        opacity = 1,
        scaling = list(
          label = list(
            enabled = TRUE,
            min = 20,
            max = 50,
            maxVisible = 10000
          )
        )
      ) %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE, hideColor = "rgba(200,200,200,0.2)"),
        nodesIdSelection = list(
          enabled = TRUE,
          values = sort(unique(filtered_nodes$id)) # Explicitly sort the IDs here
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
}

shinyApp(ui = ui, server = server)
