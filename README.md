# Extended Analysis of Leader et al.

This repository contains an extended analysis of the Leader et al. study, focusing on machine learning approaches to analyze metabolic gene expression patterns in tumor samples.

## Repository Structure

```
IMM/
├── original/           # Original repository content from Leader et al.
│   ├── data/          # Original datasets
│   ├── input_tables/  # Supporting tables
│   └── output/        # Original analysis output
├── src/               # Extended analysis code
│   ├── scripts/       # Analysis scripts
│   │   ├── DE/       # Differential expression analysis
│   │   ├── PE/       # Pathway enrichment analysis
│   │   ├── kegg/     # KEGG metabolic gene extraction
│   │   ├── modeling/ # Machine learning models
│   │   ├── preprocessing/ # Data transformation scripts
│   │   └── utils/    # General utility functions
│   └── output/       # Output from extended analysis
└── README.md         # This file
```

## Extended Analysis

The extended analysis in this repository focuses on:
- Preprocessing of bulk and single-cell RNA sequencing data
- Machine learning models to analyze metabolic gene patterns
- Feature importance analysis using SHAP values
- Comparative analysis across different sample groups

## Original Study

This work extends the analysis from Leader et al. (https://pmc.ncbi.nlm.nih.gov/articles/PMC8728963/).

## Author

Matthijs Hulsebos