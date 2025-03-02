# This script filters non-metabolic genes and healthy samples from the pseudo bulk dataset 

library(Matrix)

# Load datasets
if (!exists("lung_ldm")) {
  message("Loading Mt. Sinai scRNAseq data into R")
  load("data/lung_ldm.rd")
}

pseudo_bulk <- read.csv("src_output/pb_sample.csv",r=1,h=1,stringsAsFactors = F, check.names = FALSE)
sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)
imm_genes <- read.csv("src_output/hsa01100_genes.csv",r=1,h=1,stringsAsFactors = F)

# Filter expression data on metabolic genes (rows) and tumor samples (columns)
pseudo_bulk_filtered <- pseudo_bulk[intersect(unlist(imm_genes), rownames(pseudo_bulk)), ]

# Write to file
write.csv(
  pseudo_bulk_filtered,
  file = "src_output/pb_sample_metabolic.csv"
)

# Clear memory
rm(lung_ldm)
