# This script converts the leader_et_al dataset from scRNAseq to bulk

library(Matrix)

# Load dataset if it is not in environment
if (!exists("lung_ldm")) {
  message("Loading Mt. Sinai scRNAseq data into R")
  load("data/lung_ldm.rd")
}

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

# Reuse colnames from model_mat to ensure no 'X' prefix
colnames(pseudo_bulk_df) <- as.character(colnames(model_mat))

# Write to file
write.csv(
  pseudo_bulk_df,
  file = "src_output/pb_sample.csv"
)

# Clear memory
rm("lung_ldm")
