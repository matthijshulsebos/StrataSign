# This script converts pseudo bulk data at sample level to patient level

library(Matrix)

# Load the datasets. Use pseudo_bulk.csv if you do not want filtered genes and samples
pb_sample_filtered <- read.csv("src_output/pb_sample_filtered.csv",r=1,h=1,stringsAsFactors = F, check.names = FALSE)
sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)

# Create a group factor linking sample ID and patient ID
group_labels <- setNames(sample_annots$patient_ID, rownames(sample_annots))
group_labels <- factor(group_labels)

# Extract grouping variable (patient IDs in the order of the samples)
grouping <- group_labels[colnames(pb_sample_filtered)]

# Transpose data to [samples x genes] for aggregation
pb_sample_filtered_t <- t(pb_sample_filtered)

# Aggregate data by patient
pb_patient_filtered_t <- rowsum(pb_sample_filtered_t, group = grouping)

# Transpose back
pb_patient_filtered <- t(pb_patient_filtered_t)

# Write to file
write.csv(
  pb_patient_filtered,
  file = "src_output/pb_patient_filtered.csv"
)
