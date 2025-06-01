# Models

This directory contains the implementation of different machine learning models and the main training pipeline script.

## Training Pipeline Architecture (`train_pipeline.R`)

The `train_pipeline.R` script orchestrates the model training process. Here's how it works:

1.  **Configuration Loading:**
    *   It reads parameters from `config.yaml` located in the same directory. This file defines which normalization types, cell types, gene types, and models to process.
    *   The `extract_active_items` helper function parses the configuration to get lists of active parameters.

2.  **Path Management:**
    *   It utilizes `src/0. utils/path_manager.R` for creating standardized output directories for model results via the `create_model_output_dir` function.

3.  **Model Script Sourcing:**
    *   The `source_all_model_scripts` function dynamically finds and sources all `.R` files located in the `src/2. models/modeling/` subdirectory. Each of these files is expected to define a model training function.

4.  **Iterative Training Loop:**
    *   The script creates a grid of all combinations of active `norm_type`, `cell_type`, and `gene_type` from the configuration.
    *   It then loops through each unique combination:
        *   **Data Loading:** For each combination, it constructs paths to the preprocessed `X_train`, `X_test`, `y_train`, and `y_test` CSV files located under `output/1. data preprocessing/training datasets/<norm_type>/<cell_type>/<gene_type>/`. The `load_dataset_csvs` function reads these files.
        *   **Model Iteration:** Within each data combination loop, it iterates through the list of models specified in `config.yaml`.
            *   If a model is marked with `train: true`, the pipeline attempts to train it.
            *   **Dynamic Function Call:** It constructs the model training function name (e.g., `train_elasticnet_model` for a model named "elasticnet") and calls this function, passing the loaded datasets (`X_train_df`, `X_test_df`, `y_train_df`, `y_test_df`).
            *   **Error Handling:** `tryCatch` blocks are used around data loading and model training to allow the pipeline to continue with other combinations/models if one fails.

5.  **Results Handling:**
    *   The `results_to_output` function takes the list returned by a model's training function.
    *   It saves the following artifacts to the directory created by `create_model_output_dir` (typically `output/2. models/<norm_type>/<cell_type>/<gene_type>/<model_name>/`):
        *   The trained model object (e.g., `elasticnet_model_<cell_type>_<gene_type>.rds`).
        *   A CSV file with predictions (e.g., `predictions_<cell_type>_<gene_type>.csv`).
        *   A CSV file with feature importance scores (e.g., `feature_importance_<cell_type>_<gene_type>.csv`).

## Guidelines for Model Training Scripts

To integrate a new model into the pipeline, its corresponding R script (placed in `src/2. models/modeling/`) must adhere to the following structure and conventions:

1.  **File Location:**
    *   The script must be located in the `src/2. models/modeling/` directory.

2.  **Naming Convention:**
    *   The R script file name (e.g., `newmodel.R`) should ideally match the `name` attribute used for that model in the `config.yaml` (e.g., if `name: "newmodel"` in `config.yaml`, the script could be `newmodel.R`).

3.  **Main Training Function:**
    *   Each script must define a primary training function.
    *   The name of this function **must** follow the pattern: `train_<model_name>_model`, where `<model_name>` is the `name` specified in `config.yaml`.
        *   Example: If `config.yaml` has `name: "elasticnet"`, the function must be `train_elasticnet_model`.

4.  **Function Signature:**
    *   The training function must accept exactly four arguments in the following order:
        1.  `X_train_df`: A data frame of training features.
        2.  `X_test_df`: A data frame of test features.
        3.  `y_train_df`: A single-column data frame of training labels.
        4.  `y_test_df`: A single-column data frame of test labels.

5.  **Return Value (List):**
    *   The training function **must** return a named list containing the following three elements:
        *   `model_object`: The trained model object itself. This object should be directly savable by R's `saveRDS()` function.
        *   `predictions_df`: A data frame containing predictions on the test set. It **must** include:
            *   `y_test`: A column with the true numeric labels (0 or 1) from the test set.
            *   `y_pred`: A column with the predicted numeric class labels (0 or 1) for the test set.
            *   It *may optionally* include `y_pred_prob` for predicted probabilities if the model provides them.
        *   `raw_feature_importance`: A data frame containing feature importance scores. It **must** include:
            *   `Feature`: A column with the names of the features.
            *   `Value`: A column with the corresponding importance score or coefficient for each feature.

6.  **Model-Specific Preprocessing:**
    *   The model script is responsible for any preprocessing steps specific to that model (e.g., scaling, one-hot encoding if necessary, specific data transformations). The input data frames (`X_train_df`, etc.) are provided as loaded from the CSVs for the current data combination.
    *   For example, the `elasticnet.R` script includes its own `preprocess_data` function for tasks like NZV filtering and factor-to-numeric conversion of the target variable.

7.  **Dependencies:**
    *   Any R packages required by the model (e.g., `glmnet`, `keras`, `randomForest`) should be loaded using `library()` at the top of the model script. The pipeline assumes these packages are already installed in the R environment.

8.  **Error Handling:**
    *   While the main pipeline uses `tryCatch` to handle errors during the call to the model training function, model scripts can implement their own internal error checks and use `stop()` for critical issues that prevent successful training or result generation.

By following these guidelines, new models can be seamlessly integrated into the automated training pipeline.