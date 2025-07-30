
# Combines gene and cluster components into a feature identifier
create_feature_identifier <- function(gene, cluster_name_part, cluster_id_numeric) {
  # Check for data types of components
  if (!is.character(gene) || !is.character(cluster_name_part) || !is.numeric(cluster_id_numeric)) {
    stop("Inputs must be char, char, num.")
  }
  
  # Create the features identifier string by combining
  cluster_identifier_string <- paste(cluster_name_part, cluster_id_numeric, sep = "_")
  feature_name <- paste(gene, cluster_identifier_string, sep = "@")

  return(feature_name)
}


# Deconstructs feature identifiers into their components
parse_feature_identifier <- function(feature_identifiers) {
  # Apply the following to all feature identifiers in the input vector
  results <- lapply(feature_identifiers, function(fi_original) {
    fi <- fi_original

    # Init variables
    gene_val <- NA
    full_cluster_string_val <- NA
    cluster_name_part_val <- NA
    cluster_id_val <- NA

    # Check that feature contains @ for gene and cluster split
    if (is.character(fi) && !is.na(fi) && grepl("@", fi, fixed = TRUE)) {
      parts <- strsplit(fi, "@", fixed = TRUE)[[1]]

      # First part of split is gene and second is cluster string
      if (length(parts) == 2) {
        gene_val <- parts[1]
        full_cluster_string_val <- parts[2]

        if (nchar(full_cluster_string_val) > 0) {
          # Regex to find name and a ID by checking if there is a number after underscore
          match_data <- regexec("^(.*?)_(\\d+)$", full_cluster_string_val)

          # If regex match is successful
          if (match_data[[1]][1] != -1) {
            # Get all part of the match
            extracted_parts <- regmatches(full_cluster_string_val, match_data)[[1]]
            # Part is the cluster name and the second is the cluster ID
            cluster_name_part_val <- extracted_parts[2] 
            cluster_id_val <- as.numeric(extracted_parts[3])
          } else {
            # Fall back on entire cluster string if no cluster ID is found
            cluster_name_part_val <- full_cluster_string_val
          }
        } else { 
            # Edge case when @ is trailed by nothing
            cluster_name_part_val <- "" 
        }
      }
    }
    
    # Return all the components as a dataframe
    return(data.frame(
      feature_identifier = fi_original,
      gene = gene_val,
      full_cluster_string = full_cluster_string_val,
      cluster_name_part = cluster_name_part_val,
      cluster_id = cluster_id_val,
      stringsAsFactors = FALSE
    ))
  })
  
  # Combine all parsed feature identifier rows into a data frame
  do.call(rbind, results)
}


# Extracts gene from feature identifier
get_gene_from_feature <- function(feature_identifier) {
  # Apply this function to each feature identifier in input list
  sapply(feature_identifier, function(fi) {
    # Check if feature identifier contains @
    if (!is.character(fi) || is.na(fi) || !grepl("@", fi, fixed = TRUE)) return(NA)
    # Split by @ and return first part which is the gene
    strsplit(fi, "@", fixed = TRUE)[[1]][1]
  }, USE.NAMES = FALSE)
}

# Extracts full cluster string from feature identifier
get_cluster_string_from_feature <- function(feature_identifier) {
  # Apply this function to each feature identifier in input list
  sapply(feature_identifier, function(fi) {
    # Check if feature identifier contains @
    if (!is.character(fi) || is.na(fi) || !grepl("@", fi, fixed = TRUE)) return(NA)
    # Split by @ and return second part which is the cluster string
    parts <- strsplit(fi, "@", fixed = TRUE)[[1]]
    if (length(parts) == 2) parts[2] else NA
  }, USE.NAMES = FALSE)
}

# Gets simplified (sub)lineage name from a feature identifier
get_simplified_sublineage <- function(identifier, default_value = "None") {
  # Helper function to process a single identifier string
  process_single_id <- function(id) {
    # Check if feature identifier is a character
    if (!is.character(id) || is.na(id)) return(default_value)
    
    # Init cluster variable to parameter
    cluster_string_part <- id
    
    # Check if identifier contains @
    if (grepl("@", id, fixed = TRUE)) { 
      # Get cluster parts by splitting on @
      parts <- strsplit(id, "@", fixed = TRUE)[[1]]
      
      # If the split was succesful it has two parts then focus on the second part
      if (length(parts) == 2 && nchar(parts[2]) > 0) {
        # Second part is the cluster string
        cluster_string_part <- parts[2]
      } else {
        # Default is none or the passed argument
        return(default_value)
      }
    }
    
    # If cluster string is empty return default value
    if (nchar(cluster_string_part) == 0) {
        return(default_value)
    }

    # Regex to remove cluster ID after the cluster name
    simplified <- gsub("_[0-9]+(_[a-zA-Z]+)?$", "", cluster_string_part)
    
    # If simplified is empty return default value
    if (nchar(simplified) == 0) {
        return(default_value) 
    }
    
    return(simplified)
  }
  
  # Same as earlier function apply to all identifiers
  sapply(identifier, process_single_id, USE.NAMES = FALSE)
}

# Get cluster ID from feature identifier
get_cluster_id_from_feature <- function(feature_identifier) {
  # Apply this function to each feature identifier in input list
  sapply(feature_identifier, function(fi) {
    # Check if feature identifier contains @
    if (!is.character(fi) || is.na(fi) || !grepl("@", fi, fixed = TRUE)) return(NA)
    
    # Split by @
    parts <- strsplit(fi, "@", fixed = TRUE)[[1]]

    # Check if split was successful
    if (length(parts) != 2) return(NA)
    
    # Get the second part which is the cluster string
    cluster_string_part <- parts[2]

    # If the cluster string part is empty return NA
    if (nchar(cluster_string_part) == 0) return(NA)

    # Regex to find a numeric ID at the end of the cluster string
    match_data <- regexec(".*?_(\\d+)$", cluster_string_part)
    
    # If regex match is successful
    if (match_data[[1]][1] != -1) {
      # Extract the numeric ID from the match
      extracted_parts <- regmatches(cluster_string_part, match_data)[[1]]

      # Return the final part which is the cluster ID
      return(extracted_parts[2])
    } else {
      # If regex match fails then return NA
      return(NA)
    }
  }, USE.NAMES = FALSE)
}
