# Supported models & tuning cheat-sheet

`ukb_ml_supported_models(outcome_type)` returns the live registry. Below is
the static snapshot — verify in-session before generating code.

## Classification / regression (R/ml_model.R)

| `model` | Label | Outcomes | Package | Default tuning params |
|---------|-------|----------|---------|------------------------|
| `logistic` | Logistic regression | binary | `stats` (glm) | *(none)* |
| `linear` | Linear regression | continuous | `stats` (lm) | *(none)* |
| `rf` | Random Forest | all | `ranger` | `num.trees`, `mtry`, `min.node.size` |
| `xgboost` | XGBoost | all | `xgboost` | `nrounds`, `max_depth`, `eta` |
| `glmnet` | Elastic Net | all | `glmnet` | `alpha` |
| `svm` | Support Vector Machine | all | `e1071` | `cost`, `gamma` |
| `nnet` | Neural Network | all | `nnet` | `size`, `decay` |
| `rpart` | Decision Tree | all | `rpart` | `cp`, `maxdepth`, `minsplit` |
| `naive_bayes` | Naive Bayes | binary/multiclass | `e1071` | `laplace` |

## Survival (R/ml_survival.R)

| `model` | Label | Package | Default tuning params |
|---------|-------|---------|------------------------|
| `rsf` | Random Survival Forest | `randomForestSRC` | `ntree`, `mtry`, `nodesize` |
| `gbm_surv` | Gradient Boosting | `gbm` | `n.trees`, `interaction.depth`, `shrinkage` |
| `coxnet` | Penalised Cox | `glmnet` | `alpha`, `lambda` |

## `param_grid` examples

### XGBoost (binary)

```r
tune_params = list(
  search = "grid",
  param_grid = expand.grid(
    nrounds   = c(100, 200, 400),
    max_depth = c(3, 6, 9),
    eta       = c(0.05, 0.1, 0.2)
  ),
  resampling = "validation",
  metric     = "auc", maximize = TRUE
)
```

### glmnet logistic (random search)

```r
tune_params = list(
  search = "random",
  n_iter = 30,
  param_space = list(
    alpha = c(0, 1)              # uniform between 0 and 1
  ),
  resampling = "cv", folds = 5,
  metric = "auc", maximize = TRUE
)
```

### Random Survival Forest

```r
tune_params = list(
  search = "grid",
  param_grid = expand.grid(
    ntree    = c(500, 1000),
    mtry     = c(3, 5, 8),
    nodesize = c(5, 15, 30)
  ),
  metric = "c_index"
)
```

## Metrics

| `metric` | When |
|----------|------|
| `"auc"` | binary classification (default) |
| `"logloss"` | calibrated probability scoring |
| `"accuracy"` | balanced classes |
| `"f1"` | imbalanced, want PR-aware |
| `"brier"` | proper scoring rule, calibration-aware |
| `"rmse"` | continuous regression |
| `"mae"` | continuous regression, robust |
| `"c_index"` | survival (Harrell's C) |

## Choosing the search method

| Need | Search |
|------|--------|
| ≤ 20 grid points | `"grid"` |
| Large grid, time-budgeted | `"random"` with `n_iter` |
| Smooth continuous param space | `"bayes"` (requires `rBayesianOptimization`) |

## Choosing `resampling`

| Setting | When |
|---------|------|
| `"validation"` | Train / validation / test split available — fastest |
| `"cv"` with `folds = 5` | Smaller cohort or want variance estimate |

## Quick decision tree

1. Continuous outcome → `linear` (baseline) + `xgboost` (flexible).
2. Binary outcome, well-balanced, interpretability matters → `logistic` +
   `glmnet`.
3. Binary outcome, imbalanced, performance matters → `xgboost` with
   class-weighting via `param_grid$scale_pos_weight`.
4. Multiclass → `xgboost` (objective `multi:softprob`) or `rf`.
5. Survival, semi-parametric inference → `coxnet`.
6. Survival, non-linear interactions, lots of features → `rsf`.
