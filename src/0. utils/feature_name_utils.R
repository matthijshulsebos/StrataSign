# Creates a feature identifier gene@cluster_name_part_cluster_id
create_feature_identifier <- function(gene, cluster_name_part, cluster_id_numeric) {
  # Input type validation
  if (!is.character(gene) || !is.character(cluster_name_part) || !is.numeric(cluster_id_numeric)) {
    stop("Inputs must be: character, character, numeric.")
  }
  # Combine cluster name part and ID
  cluster_identifier_string <- paste(cluster_name_part, cluster_id_numeric, sep = "_")
  # Combine gene with the full cluster identifier
  feature_name <- paste(gene, cluster_identifier_string, sep = "@")
  return(feature_name)
}

# Parses feature identifiers and returns gene, cluster_string, cluster_name, cluster_id
parse_feature_identifier <- function(feature_identifiers) {
  results <- lapply(feature_identifiers, function(fi_original) {
    fi <- fi_original

    # Init vars
    gene_val <- NA_character_
    full_cluster_string_val <- NA_character_
    cluster_name_part_val <- NA_character_
    cluster_id_val <- NA_real_

    # Check that feature contains @
    if (is.character(fi) && !is.na(fi) && grepl("@", fi, fixed = TRUE)) {
      parts <- strsplit(fi, "@", fixed = TRUE)[[1]]
      # Check for gene and cluster_string
      if (length(parts) == 2) {
        gene_val <- parts[1]
        full_cluster_string_val <- parts[2]

        # Process the cluster string
        if (nchar(full_cluster_string_val) > 0) {
          # Regex to find name and a ID
          match_data <- regexec("^(.*?)_(\\d+)$", full_cluster_string_val)
          if (match_data[[1]][1] != -1) {
            extracted_parts <- regmatches(full_cluster_string_val, match_data)[[1]]
            cluster_name_part_val <- extracted_parts[2] 
            cluster_id_val <- as.numeric(extracted_parts[3])
          } else {
            # Default to entire cluster string if no ID found
            cluster_name_part_val <- full_cluster_string_val
          }
        } else { 
            # Edge case when @ is trailed by nothing
            cluster_name_part_val <- "" 
        }
      }
    }
    
    # Return a data frame row for each input identifier
    return(data.frame(
      feature_identifier = fi_original,
      gene = gene_val,
      full_cluster_string = full_cluster_string_val,
      cluster_name_part = cluster_name_part_val,
      cluster_id = cluster_id_val,
      stringsAsFactors = FALSE
    ))
  })
  
  # Combine all rows into a data frame
  do.call(rbind, results)
}


# Extracts gene name from feature identifier
get_gene_from_feature <- function(feature_identifier) {
  sapply(feature_identifier, function(fi) {
    # Basic validation and check for @
    if (!is.character(fi) || is.na(fi) || !grepl("@", fi, fixed = TRUE)) return(NA_character_)
    # Split by @ and take the first part (gene)
    strsplit(fi, "@", fixed = TRUE)[[1]][1]
  }, USE.NAMES = FALSE)
}

# Extracts full cluster string from feature identifier
get_cluster_string_from_feature <- function(feature_identifier) {
  sapply(feature_identifier, function(fi) {
    # Check for @
    if (!is.character(fi) || is.na(fi) || !grepl("@", fi, fixed = TRUE)) return(NA_character_)
    parts <- strsplit(fi, "@", fixed = TRUE)[[1]]
    # Return the second part the cluster string
    if (length(parts) == 2) parts[2] else NA_character_
  }, USE.NAMES = FALSE)
}

# Gets simplified (sub)lineage name from a feature identifier
get_simplified_sublineage <- function(identifier, default_value = "None") {
  
  # Helper function to process a single identifier string
  process_single_id <- function(id) {
    # Basic input validation
    if (!is.character(id) || is.na(id)) return(default_value)
    
    cluster_string_part <- id
    # If it contains @ then split
    if (grepl("@", id, fixed = TRUE)) { 
      parts <- strsplit(id, "@", fixed = TRUE)[[1]]
      # Retrieve the second part
      if (length(parts) == 2 && nchar(parts[2]) > 0) {
        cluster_string_part <- parts[2]
      } else { 
        # Nothing trails @ so something is wrong
        return(default_value)
      }
    }
    
    # If the cluster string is empty
    if (nchar(cluster_string_part) == 0) {
        return(default_value)
    }

    # Regex to remove cluster ID
    simplified <- gsub("_[0-9]+(_[a-zA-Z]+)?$", "", cluster_string_part)
    
    # If (sub)lineage name is empty
    if (nchar(simplified) == 0) {
        return(default_value) 
    }
    
    return(simplified)
  }
  
  # Apply to each element in input
  sapply(identifier, process_single_id, USE.NAMES = FALSE)
}
