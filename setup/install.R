# Function to install packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    print(paste("Installing these packages:", paste(new_packages, collapse=", ")))
    for(package in new_packages) {
      tryCatch({
        install.packages(package, repos="https://cloud.r-project.org")
        print(paste("Successfully installed", package))
      }, error = function(e) {
        print(paste("Error installing", package, ":", e$message))
      })
    }
  } else {
    print("All required packages are already installed:")
    print(packages)
  }
}

# Function to install Bioconductor packages if not already installed
install_bioc_if_missing <- function(packages) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    print("Installed BiocManager")
  }
  
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    print(paste("Installing these Bioconductor packages:", paste(new_packages, collapse=", ")))
    for(package in new_packages) {
      tryCatch({
        BiocManager::install(package, update = FALSE, ask = FALSE)
        print(paste("Successfully installed", package))
      }, error = function(e) {
        print(paste("Error installing", package, ":", e$message))
      })
    }
  } else {
    print("All required Bioconductor packages are already installed:")
    print(packages)
  }
}

# List of required CRAN packages
required_packages <- c(
  "dplyr",
  "tidyverse",
  "readr",
  "igraph",
  "visNetwork",
  "shiny",
  "DT",
  "png",
  "RCurl",
  "XML"
)

# List of required Bioconductor packages
required_bioc_packages <- c(
  "pathview",
  "KEGGREST",
  "graphite",
  "KEGGgraph",
  "org.Hs.eg.db"
)

# Install packages
print("Checking and installing required CRAN packages...")
install_if_missing(required_packages)

print("\nChecking and installing required Bioconductor packages...")
install_bioc_if_missing(required_bioc_packages)

# Verify installations
print("\nVerifying installations...")
all_packages <- c(required_packages, required_bioc_packages)
installed <- all_packages %in% installed.packages()[,"Package"]
for(i in seq_along(all_packages)) {
  status <- if(installed[i]) "installed" else "NOT installed"
  print(paste(all_packages[i], "-", status))
}

print("\nSetup complete. Please restart R session if any new packages were installed.") 