# Setup Instructions

This directory contains scripts for setting up the required R environment for the StrataSign project.

## Installation Script (`install.R`)

The `install.R` script automatically installs all required R packages from both CRAN and Bioconductor repositories.

## Usage

1. Make sure R is installed on your system
2. Open R or RStudio
3. Run the installation script:

```r
source("setup/install.R")
```

The script will:
1. Check for missing packages
2. Install any packages that aren't already present
3. Verify successful installation of all packages
4. Prompt you to restart R if any new packages were installed

## Troubleshooting

If you encounter any installation errors:
1. Check your internet connection
2. Ensure you have write permissions to your R library
3. For Windows users: Make sure Rtools is installed if compiling packages from source
4. Check the error messages in the console for specific package installation issues

## Dependencies

- R version 4.0.0 or higher
- BiocManager (will be installed automatically if missing)