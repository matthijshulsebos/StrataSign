
library(Matrix)
library(data.table)
library(glmnet)
library(caret)

# Load datasets
if (!exists("lung_ldm")) {
  message("Loading Mt. Sinai scRNAseq data into R")
  load("data/lung_ldm.rd")
}

# Load sparse matrix and target data
counts <- lung_ldm$dataset$counts  # Sparse matrix with dimensions 91 x 33660 x 60

# Load target data
target_data <- fread("src_output/lcam_sample.csv")  # Load CSV
y <- target_data$difference  # Target column to predict

# Reshape the counts matrix into 2D (flattening)
dim(counts) <- c(nrow(counts), prod(dim(counts)[-1]))  # 91 x (33660 * 60)

# Ensure that the dimensions match
if (nrow(counts) != length(y)) {
  stop("Mismatch in the number of samples between counts and target data!")
}

# Split data into training and test sets
set.seed(123)  # For reproducibility
train_index <- createDataPartition(y, p = 0.8, list = FALSE)

# Training and test sets
X_train <- counts[train_index, ]
y_train <- y[train_index]
X_test <- counts[-train_index, ]
y_test <- y[-train_index]

# Train a Lasso regression model using glmnet
lasso_model <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian", standardize = FALSE)

# Extract the best lambda value
best_lambda <- lasso_model$lambda.min
cat("Best Lambda:", best_lambda, "\n")

# Evaluate the model on the test set
predictions <- predict(lasso_model, s = "lambda.min", newx = X_test)

# Calculate RMSE (Root Mean Square Error)
rmse <- sqrt(mean((y_test - predictions)^2))
cat("RMSE:", rmse, "\n")

# Optional: Save the model and results for future use
saveRDS(lasso_model, "lasso_model.rds")
write.csv(data.frame(True = y_test, Predicted = as.vector(predictions)), "predictions.csv")
