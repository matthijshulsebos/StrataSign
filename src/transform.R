# This script converts the leader_et_al dataset from scRNAseq to bulk

library(Matrix)

if (!exists("lung_ldm")) {
  message("Loading Mt. Sinai scRNAseq data into R")
  load("data/lung_ldm.rd")
}

sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)

group_labels <- lung_ldm$dataset$cell_to_sample
group_labels <- factor(group_labels)

# Create a one-hot encoding model matrix 
model_mat <- model.matrix(~ 0 + group_labels)
colnames(model_mat) <- levels(group_labels)

# Multiply to aggregate counts by sample
pseudo_bulk <- lung_ldm$dataset$umitab %*% model_mat

# Convert to a dense matrix 
pseudo_bulk_matrix <- as.matrix(pseudo_bulk)
pseudo_bulk_df <- as.data.frame(pseudo_bulk_matrix)

# Write to file
write.csv(
  pseudo_bulk_df,
  file = "src_output/pseudo_bulk.csv"
)

rm("lung_ldm")
