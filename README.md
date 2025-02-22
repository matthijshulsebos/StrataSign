<p align="center">
  <img src="stratasign.png" alt="StrataSign Overview" width="800"/>
</p>

# Extended Analysis of Leader et al.

This repository contains an extended analysis of the Leader et al. study, focusing on machine learning approaches to analyze metabolic gene expression patterns in tumor samples.

## Repository Structure

```
StrataSign/
├── base/                  # Original repository from Leader et al.
│   ├── data/             # Sequencing data
│   ├── input_tables/     # Cluster annotations and metadata
│   ├── intermediates/    # Intermediate files
│   ├── output/           # Paper figures
│   └── scripts/          # Reproduction code
│
├── data/                 # Extended data processing
│   ├── ablation/         # Machine learning ablation study
│   ├── bulk_rna/         # Pseudobulk processed data
│   ├── kegg/            # Metabolic genes and pathways
│   ├── lcam/            # LCAM score calculations
│   ├── temp/            # Temporary files
│   └── validation/      # Data reconstruction validation
│
├── results/             # Analysis outputs
│   ├── ablation/        # Ablation study results
│   │   ├── figures/     # Visualizations
│   │   ├── models/      # Trained models
│   │   └── summary/     # Performance metrics
│   ├── analysis/        # DE and pathway enrichment
│   ├── archive/         # Archived analyses
│   └── validation/      # LCAM distributions etc.
│
├── setup/               # Setup scripts
│   └── install.R        # Package installation
│
├── src/                 # Source code (R)
│   ├── analysis/        # DE and PE analysis
│   ├── data/           # Data preprocessing
│   ├── models/         # Model training/evaluation
│   └── utils/          # Utility functions
│
└── README.md            # This file
```

## Project Components

### Base Repository
Contains the original Leader et al. study materials, including sequencing data, cluster annotations, and reproduction code.

This work extends the analysis from Leader et al. ([PMC8728963](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728963/)).

## Setup

To install required packages:
```R
source("setup/install.R")
```

## Author

Matthijs Hulsebos
