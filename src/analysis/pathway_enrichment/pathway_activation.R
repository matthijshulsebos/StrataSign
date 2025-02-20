# Load required libraries
library(dplyr)
library(tidyverse)
library(readr)
library(igraph)
library(visNetwork)
library(shiny)
library(DT)
library(KEGGREST)

# Read cluster annotations and create lineage options
cluster_annot <- read_csv("input_tables/annots_list.csv", show_col_types = FALSE) %>%
  mutate(cluster = paste0("cluster_", cluster))

# Create cache directories
dir.create("src_output/kegg_cache/enzymes", recursive = TRUE, showWarnings = FALSE)
dir.create("src_output/kegg_cache/modules", recursive = TRUE, showWarnings = FALSE)

# Function to read enrichment results for a specific cluster
read_enrichment_results <- function(cluster) {
  file_path <- file.path("src_output/cluster_de", cluster, "enrichment", "metabolic_KEGG_pathways.csv")
  print(paste("Checking file:", file_path)) # Debug print
  
  if(file.exists(file_path)) {
    data <- read_csv(file_path, show_col_types = FALSE) %>%
      mutate(
        geneID = as.character(geneID),
        cluster = cluster,
        ID = as.character(ID)
      )
    print(paste("Found", nrow(data), "pathways for cluster", cluster)) # Debug print
    if(nrow(data) > 0) {
      print(paste("First pathway:", data$Description[1])) # Debug print
      return(data)
    }
  } else {
    print(paste("File not found:", file_path)) # Debug print
  }
  return(NULL)
}

# Function to safely query KEGG with retry
safe_kegg_get <- function(id, max_retries = 3) {
  for(i in 1:max_retries) {
    tryCatch({
      Sys.sleep(0.1)  # Rate limiting
      result <- keggGet(id)[[1]]
      return(result)
    }, error = function(e) {
      if(i == max_retries) {
        print(paste("Failed to retrieve KEGG data for:", id, "-", e$message))
        return(NULL)
      }
      Sys.sleep(1)  # Wait longer before retry
    })
  }
  return(NULL)
}

# Function to get or create enzyme cache
get_cached_enzyme_info <- function(gene_id) {
  cache_file <- file.path("src_output/kegg_cache/enzymes", paste0(gene_id, "_enzymes.rds"))
  
  if(file.exists(cache_file)) {
    enzymes <- readRDS(cache_file)
    if(length(enzymes) > 0) {
      print(paste("Found cached enzymes for gene:", gene_id, "->", paste(enzymes, collapse = ", ")))
      return(enzymes)
    }
  }
  
  print(paste("Fetching enzyme info for gene:", gene_id))
  gene_info <- safe_kegg_get(gene_id)
  enzymes <- character(0)
  
  if(!is.null(gene_info)) {
    # Try direct ENZYME field
    if(!is.null(gene_info$ENZYME)) {
      enzymes <- gene_info$ENZYME
      print(paste("Found direct enzymes for gene:", gene_id, "->", paste(enzymes, collapse = ", ")))
    }
    # Try ORTHOLOGY field for enzyme information
    else if(!is.null(gene_info$ORTHOLOGY)) {
      ortho_id <- names(gene_info$ORTHOLOGY)[1]
      ortho_info <- safe_kegg_get(ortho_id)
      if(!is.null(ortho_info) && !is.null(ortho_info$ENZYME)) {
        enzymes <- ortho_info$ENZYME
        print(paste("Found orthology enzymes for gene:", gene_id, "->", paste(enzymes, collapse = ", ")))
      }
    }
  }
  
  saveRDS(enzymes, cache_file)
  return(enzymes)
}

# Function to get or create module cache
get_cached_modules <- function(gene_id) {
  cache_file <- file.path("src_output/kegg_cache/modules", paste0(gene_id, "_modules.rds"))
  
  if(file.exists(cache_file)) {
    modules <- readRDS(cache_file)
    if(length(modules) > 0) {
      print(paste("Found cached modules for gene:", gene_id, "Count:", length(modules)))
      return(modules)
    }
  }
  
  print(paste("Fetching modules for gene:", gene_id))
  modules <- list()
  
  # Get enzyme information first
  enzymes <- get_cached_enzyme_info(gene_id)
  
  if(length(enzymes) > 0) {
    # For each enzyme, get associated modules
    for(enzyme in enzymes) {
      enzyme_id <- paste0("ec:", enzyme)
      print(paste("Querying KEGG for enzyme:", enzyme_id))
      enzyme_info <- safe_kegg_get(enzyme_id)
      
      if(!is.null(enzyme_info)) {
        # Check for direct module links
        if(!is.null(enzyme_info$MODULE)) {
          for(module_id in names(enzyme_info$MODULE)) {
            if(!module_id %in% names(modules)) {
              module_data <- safe_kegg_get(module_id)
              if(!is.null(module_data)) {
                modules[[module_id]] <- list(
                  name = as.character(enzyme_info$MODULE[[module_id]]),
                  reactions = if(!is.null(module_data$REACTION)) 
                               as.character(module_data$REACTION) 
                             else character(0),
                  enzymes = enzyme
                )
                print(paste("Found module:", module_id, "for enzyme:", enzyme))
              }
            } else {
              modules[[module_id]]$enzymes <- c(modules[[module_id]]$enzymes, enzyme)
            }
          }
        }
        
        # Check reaction links
        if(!is.null(enzyme_info$REACTION)) {
          for(reaction_id in enzyme_info$REACTION) {
            print(paste("Checking reaction:", reaction_id))
            reaction_info <- safe_kegg_get(reaction_id)
            if(!is.null(reaction_info) && !is.null(reaction_info$MODULE)) {
              for(module_id in names(reaction_info$MODULE)) {
                if(!module_id %in% names(modules)) {
                  module_data <- safe_kegg_get(module_id)
                  if(!is.null(module_data)) {
                    modules[[module_id]] <- list(
                      name = as.character(reaction_info$MODULE[[module_id]]),
                      reactions = if(!is.null(module_data$REACTION)) 
                                   as.character(module_data$REACTION) 
                                 else character(0),
                      enzymes = enzyme
                    )
                    print(paste("Found module:", module_id, "through reaction for enzyme:", enzyme))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  print(paste("Found", length(modules), "modules for gene", gene_id))
  saveRDS(modules, cache_file)
  return(modules)
}

# Function to get or create global cache
get_global_cache <- function(gene_ids) {
  if(file.exists(GLOBAL_CACHE_FILE)) {
    cache <- readRDS(GLOBAL_CACHE_FILE)
    missing_genes <- setdiff(gene_ids, names(cache$genes))
    
    if(length(missing_genes) == 0) {
      print("Using global cache for all genes")
      return(cache)
    }
    
    print(paste("Found", length(missing_genes), "new genes to process"))
    return(update_global_cache(missing_genes, cache))
  }
  
  print("Creating new global cache")
  return(update_global_cache(gene_ids))
}

# Function to update global cache
update_global_cache <- function(new_genes, cache = list(genes = list(), modules = list())) {
  withProgress(message = 'Updating KEGG information...', value = 0, {
    total_genes <- length(new_genes)
    
    for(i in seq_along(new_genes)) {
      gene_id <- new_genes[i]
      incProgress(1/total_genes, 
                 detail = sprintf("Processing gene %d of %d...", i, total_genes))
      
      # Get enzyme information
      enzymes <- get_cached_enzyme_info(paste0("hsa:", gene_id))
      cache$genes[[gene_id]] <- list(enzymes = enzymes)
      
      if(length(enzymes) > 0) {
        # Process each enzyme
        for(enzyme in enzymes) {
          enzyme_id <- paste0("ec:", enzyme)
          enzyme_info <- safe_kegg_get(enzyme_id)
          
          if(!is.null(enzyme_info) && !is.null(enzyme_info$MODULE)) {
            for(module_id in names(enzyme_info$MODULE)) {
              if(!module_id %in% names(cache$modules)) {
                module_data <- safe_kegg_get(module_id)
                if(!is.null(module_data)) {
                  cache$modules[[module_id]] <- list(
                    name = as.character(enzyme_info$MODULE[[module_id]]),
                    reactions = if(!is.null(module_data$REACTION)) 
                                 as.character(module_data$REACTION) 
                               else character(0),
                    genes = gene_id,
                    enzymes = enzyme
                  )
                }
              } else {
                cache$modules[[module_id]]$genes <- unique(c(cache$modules[[module_id]]$genes, gene_id))
                cache$modules[[module_id]]$enzymes <- unique(c(cache$modules[[module_id]]$enzymes, enzyme))
              }
            }
          }
        }
      }
    }
  })
  
  saveRDS(cache, GLOBAL_CACHE_FILE)
  return(cache)
}

# Function to get module network
get_module_network <- function(enrichment_data) {
  print("Starting module network creation")
  
  # Get all metabolic pathways (including energy, lipid, amino acid metabolism etc.)
  metabolic_pathways <- enrichment_data %>%
    filter(grepl("^hsa00|^hsa01100|^hsa01200|^hsa01210|^hsa01212|^hsa01230", ID)) %>%
    arrange(p.adjust)
  
  print(paste("Found", nrow(metabolic_pathways), "metabolic pathways"))
  print("Sample pathways:")
  print(head(metabolic_pathways$Description))
  
  if(nrow(metabolic_pathways) == 0) {
    return(create_empty_network("No metabolic pathways found"))
  }
  
  # Create nodes for each pathway
  nodes <- metabolic_pathways %>%
    mutate(
      id = ID,
      label = substr(Description, 1, 30),
      group = "pathway",
      title = paste0(
        Description,
        "\n\nAdjusted p-value: ", format(p.adjust, digits = 3),
        "\nGenes: ", Count,
        "\nClusters: ", gsub(";", ", ", clusters)
      ),
      value = Count,
      significance = -log10(p.adjust)
    ) %>%
    select(id, label, group, title, value, significance)
  
  # Color nodes by significance
  nodes$color <- colorRampPalette(c("#FED976", "#FD8D3C", "#E31A1C"))(100)[
    cut(nodes$significance, breaks = 100, labels = FALSE)
  ]
  
  # Create edges between pathways that share genes
  edges <- data.frame()
  if(nrow(nodes) > 1) {
    for(i in 1:(nrow(metabolic_pathways)-1)) {
      for(j in (i+1):nrow(metabolic_pathways)) {
        genes1 <- unlist(strsplit(metabolic_pathways$geneID[i], "/"))
        genes2 <- unlist(strsplit(metabolic_pathways$geneID[j], "/"))
        shared_genes <- intersect(genes1, genes2)
        
        if(length(shared_genes) > 0) {
          edges <- rbind(edges, data.frame(
            from = metabolic_pathways$ID[i],
            to = metabolic_pathways$ID[j],
            arrows = "to",
            title = paste0("Shared genes: ", length(shared_genes)),
            width = length(shared_genes)/5 + 1,  # Scale width by number of shared genes
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  print(paste("Created network with", nrow(nodes), "nodes and", nrow(edges), "edges"))
  
  return(list(
    nodes = nodes %>% distinct(),
    edges = edges %>% distinct()
  ))
}

# Helper function to create empty network
create_empty_network <- function(message) {
  list(
    nodes = data.frame(
      id = "empty",
      label = message,
      group = "message",
      title = message,
      value = 1,
      color = "#cccccc",
      stringsAsFactors = FALSE
    ),
    edges = data.frame()
  )
}

# Helper function to create module edges
create_module_edges <- function(modules) {
  module_ids <- names(modules)
  map_df(1:(length(module_ids)-1), function(i) {
    map_df((i+1):length(module_ids), function(j) {
      shared_reactions <- intersect(
        modules[[module_ids[i]]]$reactions,
        modules[[module_ids[j]]]$reactions
      )
      if(length(shared_reactions) > 0) {
        data.frame(
          from = module_ids[i],
          to = module_ids[j],
          arrows = "to",
          title = paste0("Shared reactions: ", length(shared_reactions)),
          width = length(shared_reactions),
          stringsAsFactors = FALSE
        )
      }
    })
  })
}

# Create Shiny app
ui <- fluidPage(
  titlePanel("Metabolic Module Network"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("lineage", "Select Cell Type:",
                  choices = c("All", sort(unique(na.omit(cluster_annot$lineage))))),
      selectInput("sub_lineage", "Select Sub-type:",
                  choices = c("All", sort(unique(na.omit(cluster_annot$sub_lineage))))),
      hr(),
      verbatimTextOutput("debug_info"),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Module Network", 
                 visNetworkOutput("network_plot", height = "600px")),
        tabPanel("Enriched Pathways", 
                 DTOutput("pathway_table"))
      ),
      width = 9
    )
  )
)

server <- function(input, output, session) {
  # Reactive expression for selected clusters
  selected_clusters <- reactive({
    print("Updating selected clusters...") # Debug print
    clusters <- cluster_annot$cluster
    
    if(input$lineage != "All") {
      clusters <- cluster_annot %>%
        filter(lineage == input$lineage) %>%
        pull(cluster)
    }
    
    if(input$sub_lineage != "All") {
      clusters <- cluster_annot %>%
        filter(sub_lineage == input$sub_lineage) %>%
        pull(cluster)
    }
    
    print(paste("Selected clusters:", paste(clusters, collapse = ", "))) # Debug print
    return(clusters)
  })
  
  # Get enrichment data
  enrichment_data <- reactive({
    # Get selected clusters
    clusters <- selected_clusters()
    print(paste("Selected clusters:", paste(clusters, collapse = ", ")))
    
    # Read and combine enrichment results
    all_data <- map_df(clusters, function(cluster) {
      file_path <- file.path("src_output/cluster_de", cluster, "enrichment", "metabolic_KEGG_pathways.csv")
      print(paste("Reading file:", file_path))
      
      if(file.exists(file_path)) {
        data <- read_csv(file_path, show_col_types = FALSE)
        print(paste("Found", nrow(data), "pathways for cluster", cluster))
        if(nrow(data) > 0) {
          # Ensure geneID is character type
          data$geneID <- as.character(data$geneID)
          return(data %>% mutate(cluster = cluster))
        }
      } else {
        print(paste("File not found:", file_path))
      }
      return(NULL)
    })
    
    print(paste("Total pathways found across all clusters:", nrow(all_data)))
    
    if(nrow(all_data) > 0) {
      # Process and combine the data
      processed_data <- all_data %>%
        mutate(
          ID = gsub("path:", "", ID),
          geneID = as.character(geneID)  # Ensure consistent type
        ) %>%
        group_by(ID, Description) %>%
        summarise(
          p.adjust = min(p.adjust),
          Count = sum(Count),
          geneID = paste(unique(unlist(strsplit(geneID, "/"))), collapse = "/"),
          clusters = paste(unique(cluster), collapse = ";"),
          .groups = 'drop'
        )
      
      print(paste("Processed pathways:", nrow(processed_data)))
      print("Sample of processed data:")
      print(head(processed_data))
      
      return(processed_data)
    }
    
    return(data.frame())
  })
  
  # Create network visualization
  output$network_plot <- renderVisNetwork({
    withProgress(message = 'Building module network...', {
      data <- enrichment_data()
      
      if(nrow(data) == 0) {
        return(visNetwork(
          data.frame(id = 1, label = "No enriched pathways found", group = "message"),
          data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
        ))
      }
      
      network <- get_module_network(data)
      
      # Create the network
      visNetwork(network$nodes, network$edges) %>%
        visNodes(
          shape = "box",
          shadow = TRUE,
          font = list(size = 14)
        ) %>%
        visEdges(
          smooth = TRUE,
          color = list(color = "#666666", highlight = "#E41A1C"),
          scaling = list(min = 1, max = 10)
        ) %>%
        visOptions(
          highlightNearest = list(enabled = TRUE, degree = 1),
          nodesIdSelection = TRUE
        ) %>%
        visLayout(
          randomSeed = 123,
          improvedLayout = TRUE
        ) %>%
        visPhysics(
          solver = "forceAtlas2Based",
          forceAtlas2Based = list(
            gravitationalConstant = -50,
            centralGravity = 0.01,
            springLength = 200,
            springConstant = 0.08
          ),
          stabilization = list(
            enabled = TRUE,
            iterations = 1000
          )
        ) %>%
        visLegend(
          useGroups = FALSE,
          addNodes = data.frame(
            label = c("Most significant", "Moderately significant", "Less significant"),
            shape = "box",
            color = c("#E31A1C", "#FD8D3C", "#FED976"),
            size = 25
          ),
          addEdges = data.frame(
            label = c("Many shared reactions", "Few shared reactions"),
            width = c(8, 2),
            color = c("#666666", "#666666")
          ),
          main = "Module Network Legend",
          width = 0.2,
          position = "right",
          ncol = 1,
          stepX = 100,
          stepY = 100
        )
    })
  })
  
  # Helper function to add legend
  addLegend <- function(vis, title, position = "right", text) {
    footer <- sprintf(
      'document.getElementById("%s").insertAdjacentHTML("afterend", 
      "<div style=\'position:absolute;%s:0;top:0;padding:10px;
      background-color:rgba(255,255,255,0.8);border:1px solid #ddd;
      font-family:sans-serif;font-size:12px\'><b>%s</b><br><pre>%s</pre></div>");',
      vis$elementId,
      position,
      title,
      text
    )
    
    vis %>% onRender(footer)
  }
  
  # Show pathway table
  output$pathway_table <- renderDataTable({
    data <- enrichment_data()
    print(paste("Enrichment table rows:", nrow(data)))
    
    if(nrow(data) > 0) {
      formatted_data <- data %>%
        select(Description, p.adjust, Count, clusters) %>%  # Only select columns we know exist
        mutate(
          p.adjust = format(p.adjust, digits = 3),
          clusters = gsub(";", ", ", clusters)
        ) %>%
        rename(
          Pathway = Description,
          "Adjusted p-value" = p.adjust,
          "Gene count" = Count,
          Clusters = clusters
        ) %>%
        arrange(`Adjusted p-value`)
      
      print("Formatted enrichment table:")
      print(head(formatted_data))
      formatted_data
    } else {
      data.frame(
        Pathway = character(),
        "Adjusted p-value" = character(),
        "Gene count" = integer(),
        Clusters = character(),
        stringsAsFactors = FALSE
      )
    }
  }, options = list(
    pageLength = 15,
    order = list(list(1, 'asc')),  # Sort by p-value
    dom = 'lftip',
    lengthMenu = list(c(10, 15, 25, 50), c('10', '15', '25', '50')),
    rownames = FALSE
  ))
  
  # Debug information
  output$debug_info <- renderPrint({
    data <- enrichment_data()
    cat("Selected clusters:", paste(selected_clusters(), collapse = ", "), "\n\n")
    cat("Number of enriched pathways:", nrow(data), "\n")
    if(nrow(data) > 0) {
      cat("Sample of gene IDs:\n")
      cat(head(data$geneID, 1), "\n\n")
      total_genes <- length(unique(unlist(strsplit(data$geneID, "/"))))
      cat("Total unique genes:", total_genes, "\n\n")
      cat("Top 5 enriched pathways:\n")
      print(head(data %>% 
                  select(Description, p.adjust, geneID, clusters) %>% 
                  arrange(p.adjust), 5))
    }
  })
  
  # Update sub-lineage choices
  observe({
    if(input$lineage != "All") {
      sub_lineages <- cluster_annot %>%
        filter(lineage == input$lineage) %>%
        pull(sub_lineage) %>%
        unique() %>%
        na.omit()
      
      updateSelectInput(session, "sub_lineage",
                       choices = c("All", sort(sub_lineages)))
    } else {
      updateSelectInput(session, "sub_lineage",
                       choices = c("All", sort(unique(na.omit(cluster_annot$sub_lineage)))))
    }
  })
  
  # Update the enrichment table output
  output$enrichment_table <- renderDataTable({
    data <- enrichment_data()
    print(paste("Enrichment table rows:", nrow(data)))
    
    if(nrow(data) > 0) {
      formatted_data <- data %>%
        select(Description, p.adjust, Count, clusters) %>%
        mutate(
          p.adjust = format(p.adjust, digits = 3),
          clusters = gsub(";", ", ", clusters)
        ) %>%
        rename(
          Pathway = Description,
          "Adjusted p-value" = p.adjust,
          "Gene count" = Count,
          Clusters = clusters
        ) %>%
        arrange(`Adjusted p-value`)
      
      print("Formatted enrichment table:")
      print(head(formatted_data))
      formatted_data
    } else {
      data.frame(
        Pathway = character(),
        "Adjusted p-value" = character(),
        "Gene count" = integer(),
        Clusters = character(),
        stringsAsFactors = FALSE
      )
    }
  }, options = list(
    pageLength = 15,
    order = list(list(1, 'asc')),  # Sort by p-value (now column 1)
    dom = 'lftip',
    lengthMenu = list(c(10, 15, 25, 50), c('10', '15', '25', '50')),
    rownames = FALSE
  ))
}

# Run the app
print("Starting Metabolic Module Network Viewer...")
runApp(shinyApp(ui = ui, server = server), launch.browser = TRUE)
