# Models

This directory contains the implementation of different machine learning models and analyses for the StrataSign project.

## Directory Structure

```
models/
├── ablation/              # Ablation study implementations
│   ├── feature_groups/    # Models with different feature subsets
│   └── architectures/     # Different model architectures
└── evaluation/          # Model evaluation scripts
    ├── metrics/        # Performance metric calculations
    └── visualization/  # Performance visualization
```

## Model Implementations

### Ablation Studies
- Feature group importance
  - Metabolic pathways
  - Gene expression signatures
- Architecture comparisons
  - Model complexity analysis
  - Performance trade-offs


## Output

Model outputs are saved to:
- Trained models: `results/models/{model_type}/{date}_model.rds`
- Performance metrics: `results/analysis/{model_type}/metrics_{date}.csv`
- Feature importance: `results/analysis/{model_type}/importance_{date}.csv`
