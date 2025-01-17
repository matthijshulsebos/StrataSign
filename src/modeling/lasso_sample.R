# Load necessary libraries
library(glmnet)    # For LASSO model and cross-validation
library(dplyr)     # For data manipulation
library(data.table) # For reading large CSV files with fread
library(ggplot2)   # (Optional) For enhanced plotting

# Load samples filtered on metabolic genes
gene_expression_matrix  <- read.csv("src_output/pb_sample_metabolic.csv",r=1,h=1,stringsAsFactors = F, check.names = FALSE)
gene_expression_matrix <- t(gene_expression_matrix)

# Load target data
target_data <- fread("src_output/lcam_sample.csv")  # Load CSV

# Load sample metadata
sample_annots <- read.csv("input_tables/table_s1_sample_table.csv",r=1,h=1,stringsAsFactors = F)

# Convert rownames of sample_annots into a column named 'sample_ID'
sample_annots$sample_ID <- rownames(sample_annots)

# Ensure that V1 in target_data is a character to match sample_ID
target_data$V1 <- as.character(target_data$V1)
sample_annots$sample_ID <- as.character(sample_annots$sample_ID)

# Merge target_data and sample_annots based on sample ID
merged_data <- merge(target_data, sample_annots, by.x = "V1", by.y = "sample_ID", all.x = TRUE)

# Make sure there are no missing values in the merged data
# merged_data <- na.omit(merged_data)

# Create patient train/test split to ensure all samples from a patient are either in train or test
set.seed(123)  # For reproducibility
patients <- unique(merged_data$patient_ID)
train_patients <- sample(patients, size = floor(0.8 * length(patients)), replace = FALSE)  # 80% for training
test_patients <- setdiff(patients, train_patients)  # Remaining patients for testing

# Separate train and test sets based on patient IDs
train_samples <- merged_data %>% filter(patient_ID %in% train_patients)
test_samples <- merged_data %>% filter(patient_ID %in% test_patients)

# Subset gene expression data for train and test samples based on sample IDs
train_expression <- gene_expression_matrix[rownames(gene_expression_matrix) %in% train_samples$V1, ]
test_expression <- gene_expression_matrix[rownames(gene_expression_matrix) %in% test_samples$V1, ]

# Prepare the response variable (difference) for training
train_response <- train_samples$difference
test_response <- test_samples$difference

# Fit LASSO model on training data
lasso_model <- cv.glmnet(x = as.matrix(train_expression), y = train_response, alpha = 1)

# Plot cross-validation results to select the best lambda
plot(lasso_model)

# Make predictions on the test set
test_predictions <- predict(lasso_model, s = "lambda.min", newx = as.matrix(test_expression))

# Evaluate performance using Mean Squared Error (MSE)
mse <- mean((test_predictions - test_response)^2)
cat("Mean Squared Error: ", mse, "\n")

# Optional: Save the model for later use
save(lasso_model, file = "lasso_model.RData")

# Optionally, you can also store the predictions along with their respective sample IDs
predictions_df <- data.frame(sample_ID = test_samples$V1, predicted_difference = test_predictions)
write.csv(predictions_df, "predictions.csv", row.names = FALSE)

# Predicted vs Actual Plot
plot(test_response, test_predictions, 
     xlab = "Actual LCAM Difference", ylab = "Predicted LCAM Difference",
     main = "Predicted vs Actual LCAM Differences", 
     pch = 16, col = "blue")
abline(a = 0, b = 1, col = "red", lty = 2)  # Add a reference line y = x

# Residuals Plot
residuals <- test_predictions - test_response
plot(test_predictions, residuals, 
     xlab = "Predicted LCAM Difference", ylab = "Residuals",
     main = "Residuals Plot", 
     pch = 16, col = "blue")
abline(h = 0, col = "red", lty = 2)  # Add horizontal line at y = 0

# Histogram of Residuals
hist(residuals, breaks = 20, col = "blue", 
     main = "Histogram of Residuals", 
     xlab = "Residuals", 
     ylab = "Frequency")

# Learning Curve Plot (example for training and validation error)
train_errors <- lasso_model$cvm  # Cross-validation error for training data
validation_errors <- lasso_model$cvsd  # Standard deviation of validation errors
plot(log10(lasso_model$lambda), train_errors, type = "b", col = "blue", 
     xlab = "Log(Lambda)", ylab = "Cross-Validation Error",
     main = "Learning Curve")
lines(log10(lasso_model$lambda), validation_errors, type = "b", col = "red")
legend("topright", legend = c("Train Error", "Validation Error"), col = c("blue", "red"), lty = 1)
