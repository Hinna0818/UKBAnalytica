#' Machine Learning Module for UK Biobank Data Analysis
#'
#' @description
#' Unified machine learning interface for UK Biobank data analysis.
#' Supports classification, regression, and provides consistent API
#' across different ML algorithms.
#'
#' @name ml_model
#' @keywords internal
#' @importFrom pROC roc auc
NULL

# Helper Functions

#' Check if ML package is available
#' @keywords internal
.check_ml_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    auto_install <- isTRUE(getOption("UKBAnalytica.auto_install_ml", FALSE))
    if (auto_install) {
      message(sprintf(
        "Package '%s' is required and not installed. Installing because option UKBAnalytica.auto_install_ml = TRUE.",
        pkg
      ))
      utils::install.packages(pkg)
    }
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is required for this function. Install with install.packages('%s'), or set options(UKBAnalytica.auto_install_ml = TRUE) before calling the selected model.",
      pkg, pkg
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Warn about deprecated legacy ML interfaces
#' @keywords internal
#' @noRd
.ukb_ml_deprecated <- function(old, new) {
  warning(sprintf(
    "%s() is deprecated; please use %s instead.",
    old, new
  ), call. = FALSE)
}

#' Get model type label
#' @keywords internal
.get_model_label <- function(model) {
  labels <- c(
    rf = "Random Forest",
    xgboost = "XGBoost",
    glmnet = "Elastic Net",
    svm = "Support Vector Machine",
    nnet = "Neural Network",
    logistic = "Logistic Regression",
    rpart = "Decision Tree",
    naive_bayes = "Naive Bayes"
  )
  labels[model]
}

#' Parse formula to get response and predictors
#' @keywords internal
.parse_formula <- function(formula, data) {
  terms_obj <- terms(formula, data = data)
  response <- all.vars(formula)[1]
  predictors <- attr(terms_obj, "term.labels")
  
  # Handle . notation

if (length(predictors) == 1 && predictors == ".") {
    predictors <- setdiff(names(data), response)
  }
  
  list(response = response, predictors = predictors)
}

#' Prepare model matrix
#' @keywords internal
.prepare_model_data <- function(formula, data, task) {
  parsed <- .parse_formula(formula, data)
  
  # Remove rows with NA
  complete_idx <- complete.cases(data[, c(parsed$response, parsed$predictors), drop = FALSE])
  clean_data <- data[complete_idx, , drop = FALSE]
  
  if (nrow(clean_data) < nrow(data)) {
    message(sprintf("Removed %d rows with missing values", nrow(data) - nrow(clean_data)))
  }
  
  # Prepare X and y
  X <- clean_data[, parsed$predictors, drop = FALSE]
  y <- clean_data[[parsed$response]]
  
  # Convert factors for classification
  if (task == "classification") {
    if (!is.factor(y)) {
      y <- as.factor(y)
    }
  }
  
  list(
    X = X,
    y = y,
    data = clean_data,
    response = parsed$response,
    predictors = parsed$predictors,
    complete_idx = which(complete_idx)
  )
}

#' Split data into train/test
#' @keywords internal
.split_data <- function(data, y, split_ratio, stratify, seed) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(data)
  
  if (stratify && is.factor(y)) {
    # Stratified split
    train_idx <- integer(0)
    for (level in levels(y)) {
      level_idx <- which(y == level)
      n_train <- round(length(level_idx) * split_ratio)
      train_idx <- c(train_idx, sample(level_idx, n_train))
    }
  } else {
    n_train <- round(n * split_ratio)
    train_idx <- sample(n, n_train)
  }
  
  test_idx <- setdiff(seq_len(n), train_idx)
  
  list(train_idx = train_idx, test_idx = test_idx)
}

# Legacy implementation kept for reference; the exported implementation lives
# in R/ml_workflow.R and preserves this calling style.
.ukb_ml_split_data_legacy <- function(df,
                                      split_ratio = 0.8,
                                      stratify_by = NULL,
                                      seed = NULL,
                                      verbose = TRUE) {

  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame or data.table", call. = FALSE)
  }

  n <- nrow(df)
  if (is.null(n) || n < 2) {
    stop("`df` must contain at least 2 rows", call. = FALSE)
  }

  if (!is.numeric(split_ratio) || length(split_ratio) != 1 ||
      is.na(split_ratio) || split_ratio <= 0 || split_ratio >= 1) {
    stop("`split_ratio` must be a single numeric value strictly between 0 and 1", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(stratify_by)) {
    if (!is.character(stratify_by) || length(stratify_by) != 1) {
      stop("`stratify_by` must be NULL or a single column name", call. = FALSE)
    }
    if (!stratify_by %in% names(df)) {
      stop(sprintf("Column not found in `df`: %s", stratify_by), call. = FALSE)
    }

    strat_vec <- as.character(df[[stratify_by]])
    strat_vec[is.na(strat_vec)] <- "<NA>"

    strata <- split(seq_len(n), strat_vec)
    train_idx <- integer(0)

    for (group_idx in strata) {
      n_group <- length(group_idx)

      n_group_train <- as.integer(round(n_group * split_ratio))
      if (n_group > 1) {
        n_group_train <- max(1L, min(n_group - 1L, n_group_train))
      } else {
        n_group_train <- if (split_ratio >= 0.5) 1L else 0L
      }

      if (n_group_train > 0L) {
        train_idx <- c(train_idx, sample(group_idx, n_group_train))
      }
    }
  } else {
    n_train <- as.integer(round(n * split_ratio))
    n_train <- max(1L, min(n - 1L, n_train))
    train_idx <- sample(seq_len(n), n_train)
  }

  train_idx <- sort(unique(train_idx))
  internal_validation_idx <- setdiff(seq_len(n), train_idx)

  train_df <- df[train_idx, , drop = FALSE]
  internal_validation_df <- df[internal_validation_idx, , drop = FALSE]

  if (isTRUE(verbose)) {
    if (is.null(stratify_by)) {
      message(sprintf(
        "Split complete: train=%d, internal_validation=%d (random)",
        nrow(train_df), nrow(internal_validation_df)
      ))
    } else {
      message(sprintf(
        "Split complete: train=%d, internal_validation=%d (stratified by '%s')",
        nrow(train_df), nrow(internal_validation_df), stratify_by
      ))
    }
  }

  list(
    train = train_df,
    internal_validation = internal_validation_df
  )
}

# Main Model Training Function

#' Train a Machine Learning Model
#'
#' @description
#' Unified interface for training machine learning models on UK Biobank data.
#' Supports random forest, XGBoost, elastic net, SVM, and neural networks.
#'
#' @param formula Model formula (e.g., outcome ~ var1 + var2)
#' @param data Data frame containing variables
#' @param model Model type: "rf", "xgboost", "glmnet", "svm", "nnet", "logistic"
#' @param task Task type: "classification" or "regression"
#' @param split_ratio Train/test split ratio (default 0.8)
#' @param stratify Logical; use stratified sampling for classification (default TRUE)
#' @param seed Random seed for reproducibility
#' @param sample_n Optional; subsample data for large datasets
#' @param params List of model-specific parameters
#' @param cv Logical; perform cross-validation (default FALSE)
#' @param cv_folds Number of CV folds (default 5)
#' @param verbose Logical; print progress messages
#' @param ... Additional arguments passed to model function
#'
#' @return An object of class "ukb_ml" containing:
#' \itemize{
#'   \item model: The fitted model object
#'   \item model_type: Type of model used
#'   \item task: Task type (classification/regression)
#'   \item predictors: Names of predictor variables
#'   \item outcome: Name of outcome variable
#'   \item train_data: Training data
#'   \item test_data: Test data
#'   \item metrics: Model performance metrics
#' }
#'
#' @export
ukb_ml_model <- function(formula,
                         data,
                         model = c("rf", "xgboost", "glmnet", "svm", "nnet", "logistic"),
                         task = c("classification", "regression"),
                         split_ratio = 0.8,
                         stratify = TRUE,
                         seed = NULL,
                         sample_n = NULL,
                         params = list(),
                         cv = FALSE,
                         cv_folds = 5,
                         verbose = TRUE,
                         ...) {
  .ukb_ml_deprecated("ukb_ml_model", "ukb_ml_workflow()")
  
  model <- match.arg(model)
  task <- match.arg(task)
  
  # Sample data if requested
  if (!is.null(sample_n) && sample_n < nrow(data)) {
    if (!is.null(seed)) set.seed(seed)
    data <- data[sample(nrow(data), sample_n), ]
    if (verbose) message(sprintf("Sampled %d observations from data", sample_n))
  }
  
  # Prepare data
  prep <- .prepare_model_data(formula, data, task)
  
  # Split data
  split <- .split_data(prep$data, prep$y, split_ratio, stratify, seed)
  
  train_data <- prep$data[split$train_idx, ]
  test_data <- prep$data[split$test_idx, ]
  
  X_train <- prep$X[split$train_idx, , drop = FALSE]
  y_train <- prep$y[split$train_idx]
  X_test <- prep$X[split$test_idx, , drop = FALSE]
  y_test <- prep$y[split$test_idx]
  
  if (verbose) {
    message(sprintf("Training %s for %s", .get_model_label(model), task))
    message(sprintf("Train: %d, Test: %d observations", nrow(train_data), nrow(test_data)))
  }
  
  # Train model
  fitted_model <- switch(model,
    rf = .fit_rf(X_train, y_train, task, params, verbose, ...),
    xgboost = .fit_xgboost(X_train, y_train, task, params, verbose, ...),
    glmnet = .fit_glmnet(X_train, y_train, task, params, verbose, ...),
    svm = .fit_svm(X_train, y_train, task, params, verbose, ...),
    nnet = .fit_nnet(X_train, y_train, task, params, verbose, ...),
    logistic = .fit_logistic(formula, train_data, task, params, verbose, ...)
  )
  
  # Create result object
  result <- list(
    model = fitted_model,
    model_type = model,
    task = task,
    formula = formula,
    predictors = prep$predictors,
    outcome = prep$response,
    train_data = train_data,
    test_data = test_data,
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test,
    train_idx = split$train_idx,
    test_idx = split$test_idx,
    seed = seed,
    call = match.call()
  )
  
  class(result) <- "ukb_ml"
  
  # Calculate metrics on test set
  result$metrics <- ukb_ml_metrics(result, verbose = FALSE)
  
  if (verbose) {
    if (task == "classification") {
      message(sprintf("Test AUC: %.3f", result$metrics["auc"]))
    } else {
      message(sprintf("Test RMSE: %.3f", result$metrics["rmse"]))
    }
  }
  
  # Cross-validation if requested
  if (cv) {
    result$cv_results <- ukb_ml_cv(
      formula, prep$data, model = model, task = task,
      folds = cv_folds, seed = seed, params = params, verbose = verbose
    )
  }
  
  result
}

# Model Fitting Functions

#' Fit Random Forest
#' @keywords internal
.fit_rf <- function(X, y, task, params, verbose, ...) {
  .check_ml_package("ranger")
  
  # Default parameters
  default_params <- list(
    num.trees = 500,
    mtry = NULL,
    min.node.size = if (task == "classification") 1 else 5,
    importance = "permutation"
  )
  
  # Merge with user params
  model_params <- modifyList(default_params, params)
  
  # Prepare data
  train_df <- cbind(y = y, X)
  
  # Set probability for classification
  if (task == "classification") {
    model_params$probability <- TRUE
  }
  
  # Fit model
  do.call(ranger::ranger, c(
    list(formula = y ~ ., data = train_df),
    model_params
  ))
}

#' Fit XGBoost
#' @keywords internal
.fit_xgboost <- function(X, y, task, params, verbose, ...) {
  .check_ml_package("xgboost")
  
  # Convert to matrix
  X_mat <- as.matrix(X)
  
  # Handle factors
  for (i in seq_len(ncol(X_mat))) {
    if (is.factor(X[, i])) {
      X_mat[, i] <- as.numeric(X[, i]) - 1
    }
  }
  mode(X_mat) <- "numeric"
  
  # Prepare labels
  if (task == "classification") {
    y_numeric <- as.numeric(y) - 1
    objective <- if (length(unique(y)) == 2) "binary:logistic" else "multi:softprob"
    eval_metric <- "auc"
  } else {
    y_numeric <- as.numeric(y)
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
  }
  
  # Default parameters
  default_params <- list(
    nrounds = 100,
    max_depth = 6,
    eta = 0.3,
    objective = objective,
    eval_metric = eval_metric,
    verbose = if (verbose) 1 else 0
  )
  
  # Merge with user params
  model_params <- modifyList(default_params, params)
  
  # Extract nrounds
  nrounds <- model_params$nrounds
  model_params$nrounds <- NULL
  
  # Create DMatrix
  dtrain <- xgboost::xgb.DMatrix(data = X_mat, label = y_numeric)
  
  # Add num_class for multiclass
  if (task == "classification" && length(unique(y)) > 2) {
    model_params$num_class <- length(unique(y))
  }
  
  # Fit model
  xgboost::xgb.train(
    params = model_params,
    data = dtrain,
    nrounds = nrounds
  )
}

#' Fit GLMNet (Elastic Net)
#' @keywords internal
.fit_glmnet <- function(X, y, task, params, verbose, ...) {
  .check_ml_package("glmnet")
  
  # Convert to matrix
  X_mat <- model.matrix(~ . - 1, data = X)
  
  # Default parameters
  default_params <- list(
    alpha = 1,  # LASSO
    family = if (task == "classification") "binomial" else "gaussian"
  )
  
  # Merge with user params
  model_params <- modifyList(default_params, params)
  
  # Fit model with CV to select lambda
  do.call(glmnet::cv.glmnet, c(
    list(x = X_mat, y = y),
    model_params
  ))
}

#' Fit SVM
#' @keywords internal
.fit_svm <- function(X, y, task, params, verbose, ...) {
  .check_ml_package("e1071")
  
  # Default parameters
  default_params <- list(
    kernel = "radial",
    probability = TRUE,
    scale = TRUE
  )
  
  # Merge with user params
  model_params <- modifyList(default_params, params)
  
  # Prepare data
  train_df <- cbind(y = y, X)
  
  # Fit model
  do.call(e1071::svm, c(
    list(formula = y ~ ., data = train_df),
    model_params
  ))
}

#' Fit Neural Network
#' @keywords internal
.fit_nnet <- function(X, y, task, params, verbose, ...) {
  .check_ml_package("nnet")
  
  # Default parameters
  default_params <- list(
    size = 10,
    decay = 0.01,
    maxit = 200,
    trace = verbose
  )
  
  # Merge with user params
  model_params <- modifyList(default_params, params)
  
  # Prepare data
  train_df <- cbind(y = y, X)
  
  # Fit model
  if (task == "classification") {
    do.call(nnet::nnet, c(
      list(formula = y ~ ., data = train_df),
      model_params
    ))
  } else {
    model_params$linout <- TRUE
    do.call(nnet::nnet, c(
      list(formula = y ~ ., data = train_df),
      model_params
    ))
  }
}

#' Fit Logistic/Linear Regression
#' @keywords internal
.fit_logistic <- function(formula, data, task, params, verbose, ...) {
  if (task == "classification") {
    glm(formula, data = data, family = binomial())
  } else {
    lm(formula, data = data)
  }
}

# Prediction Function

#' Predict from ML Model
#'
#' @description
#' Generate predictions from a trained ukb_ml model.
#'
#' @param object A ukb_ml object from ukb_ml_model()
#' @param newdata Optional new data for prediction. If NULL, uses test data.
#' @param type Prediction type: "response", "prob", "class", "link"
#' @param ... Additional arguments
#'
#' @return Predictions as vector or matrix
#'
#' @export
ukb_ml_predict <- function(object, newdata = NULL, type = c("response", "prob", "class", "link"), ...) {
  .ukb_ml_deprecated("ukb_ml_predict", "predict(ml$final_model, newdata = ..., type = ...)")
  UseMethod("ukb_ml_predict")
}

#' @export
ukb_ml_predict.ukb_ml <- function(object, newdata = NULL, type = c("response", "prob", "class", "link"), ...) {
  type <- match.arg(type)
  
  if (is.null(newdata)) {
    newdata <- object$test_data
    X_new <- object$X_test
  } else {
    X_new <- newdata[, object$predictors, drop = FALSE]
  }
  
  model <- object$model
  model_type <- object$model_type
  task <- object$task
  
  # Get predictions based on model type
  pred <- switch(model_type,
    rf = .predict_rf(model, X_new, type, task),
    xgboost = .predict_xgboost(model, X_new, type, task, object),
    glmnet = .predict_glmnet(model, X_new, type, task),
    svm = .predict_svm(model, X_new, type, task),
    nnet = .predict_nnet(model, X_new, type, task),
    logistic = .predict_logistic(model, newdata, type, task)
  )
  
  pred
}

#' @keywords internal
.predict_rf <- function(model, X, type, task) {
  pred <- predict(model, data = X)
  
  if (task == "classification") {
    if (type %in% c("prob", "response")) {
      return(pred$predictions)
    } else if (type == "class") {
      return(colnames(pred$predictions)[apply(pred$predictions, 1, which.max)])
    }
  } else {
    return(pred$predictions)
  }
}

#' @keywords internal
.predict_xgboost <- function(model, X, type, task, object) {
  .check_ml_package("xgboost")
  
  X_mat <- as.matrix(X)
  for (i in seq_len(ncol(X_mat))) {
    if (is.factor(X[, i])) {
      X_mat[, i] <- as.numeric(X[, i]) - 1
    }
  }
  mode(X_mat) <- "numeric"
  
  pred <- predict(model, X_mat)
  
  if (task == "classification") {
    if (type == "class") {
      return(ifelse(pred > 0.5, levels(object$y_train)[2], levels(object$y_train)[1]))
    } else {
      # Return probability matrix
      return(cbind(1 - pred, pred))
    }
  } else {
    return(pred)
  }
}

#' @keywords internal
.predict_glmnet <- function(model, X, type, task) {
  .check_ml_package("glmnet")
  
  X_mat <- model.matrix(~ . - 1, data = X)
  
  if (task == "classification") {
    if (type == "class") {
      predict(model, newx = X_mat, s = "lambda.min", type = "class")
    } else {
      predict(model, newx = X_mat, s = "lambda.min", type = "response")
    }
  } else {
    predict(model, newx = X_mat, s = "lambda.min")
  }
}

#' @keywords internal
.predict_svm <- function(model, X, type, task) {
  if (task == "classification") {
    if (type == "class") {
      predict(model, newdata = X)
    } else {
      attr(predict(model, newdata = X, probability = TRUE), "probabilities")
    }
  } else {
    predict(model, newdata = X)
  }
}

#' @keywords internal
.predict_nnet <- function(model, X, type, task) {
  if (task == "classification") {
    pred <- predict(model, newdata = X, type = "raw")
    if (type == "class") {
      colnames(pred)[apply(pred, 1, which.max)]
    } else {
      pred
    }
  } else {
    predict(model, newdata = X)
  }
}

#' @keywords internal
.predict_logistic <- function(model, newdata, type, task) {
  if (task == "classification") {
    if (type == "class") {
      pred <- predict(model, newdata = newdata, type = "response")
      ifelse(pred > 0.5, 1, 0)
    } else {
      pred <- predict(model, newdata = newdata, type = "response")
      cbind(1 - pred, pred)
    }
  } else {
    predict(model, newdata = newdata)
  }
}

# Variable Importance

#' Get Variable Importance
#'
#' @description
#' Extract variable importance from a trained ML model.
#'
#' @param object A ukb_ml object
#' @param type Importance type (model-specific)
#' @param ... Additional arguments
#'
#' @return Data frame with variable importance scores
#'
#' @export
ukb_ml_importance <- function(object, type = NULL, ...) {
  .ukb_ml_deprecated("ukb_ml_importance", "ukb_shap() or ml$feature_result$selected_features")
  
  model <- object$model
  model_type <- object$model_type
  
  importance <- switch(model_type,
    rf = {
      imp <- ranger::importance(model)
      data.frame(
        variable = names(imp),
        importance = as.numeric(imp),
        stringsAsFactors = FALSE
      )
    },
    xgboost = {
      .check_ml_package("xgboost")
      imp <- xgboost::xgb.importance(model = model)
      data.frame(
        variable = imp$Feature,
        importance = imp$Gain,
        stringsAsFactors = FALSE
      )
    },
    glmnet = {
      coefs <- as.matrix(coef(model, s = "lambda.min"))
      data.frame(
        variable = rownames(coefs)[-1],
        importance = abs(coefs[-1, 1]),
        stringsAsFactors = FALSE
      )
    },
    logistic = {
      coefs <- coef(model)
      data.frame(
        variable = names(coefs)[-1],
        importance = abs(coefs[-1]),
        stringsAsFactors = FALSE
      )
    },
    {
      warning("Variable importance not available for this model type")
      return(NULL)
    }
  )
  
  # Sort by importance
  importance <- importance[order(importance$importance, decreasing = TRUE), ]
  rownames(importance) <- NULL
  
  importance
}

# Cross-Validation

#' Cross-Validation for ML Models
#'
#' @description
#' Perform k-fold cross-validation for ML models.
#'
#' @param formula Model formula
#' @param data Data frame
#' @param model Model type
#' @param task Task type
#' @param folds Number of folds (default 5)
#' @param repeats Number of repeats (default 1)
#' @param stratify Use stratified folds for classification
#' @param metrics Metrics to compute
#' @param params Model parameters
#' @param seed Random seed
#' @param verbose Print progress
#' @param ... Additional arguments
#'
#' @return ukb_ml_cv object with cross-validation results
#'
#' @export
ukb_ml_cv <- function(formula,
                      data,
                      model = "rf",
                      task = "classification",
                      folds = 5,
                      repeats = 1,
                      stratify = TRUE,
                      metrics = NULL,
                      params = list(),
                      seed = NULL,
                      verbose = TRUE,
                      ...) {
  .ukb_ml_deprecated("ukb_ml_cv", "ukb_ml_tune(..., resampling = 'cv') inside ukb_ml_workflow()")
  
  if (!is.null(seed)) set.seed(seed)
  
  # Prepare data
  prep <- .prepare_model_data(formula, data, task)
  
  # Default metrics
  if (is.null(metrics)) {
    metrics <- if (task == "classification") {
      c("auc", "accuracy", "sensitivity", "specificity")
    } else {
      c("rmse", "mae", "rsquared")
    }
  }
  
  all_results <- list()
  
  for (rep in seq_len(repeats)) {
    # Create folds
    if (stratify && task == "classification") {
      fold_idx <- .create_stratified_folds(prep$y, folds)
    } else {
      fold_idx <- .create_folds(nrow(prep$data), folds)
    }
    
    fold_results <- list()
    
    for (fold in seq_len(folds)) {
      if (verbose) {
        message(sprintf("Repeat %d, Fold %d/%d", rep, fold, folds))
      }
      
      # Split data
      test_idx <- fold_idx[[fold]]
      train_idx <- unlist(fold_idx[-fold])
      
      train_data <- prep$data[train_idx, ]
      
      X_train <- prep$X[train_idx, , drop = FALSE]
      y_train <- prep$y[train_idx]
      X_test <- prep$X[test_idx, , drop = FALSE]
      y_test <- prep$y[test_idx]
      
      # Fit model
      fitted <- switch(model,
        rf = .fit_rf(X_train, y_train, task, params, FALSE, ...),
        xgboost = .fit_xgboost(X_train, y_train, task, params, FALSE, ...),
        glmnet = .fit_glmnet(X_train, y_train, task, params, FALSE, ...),
        svm = .fit_svm(X_train, y_train, task, params, FALSE, ...),
        nnet = .fit_nnet(X_train, y_train, task, params, FALSE, ...),
        logistic = .fit_logistic(formula, train_data, task, params, FALSE, ...)
      )
      
      # Create temp object for prediction
      temp_obj <- list(
        model = fitted,
        model_type = model,
        task = task,
        predictors = prep$predictors,
        X_test = X_test,
        y_train = y_train
      )
      class(temp_obj) <- "ukb_ml"
      
      # Get predictions
      pred <- ukb_ml_predict(temp_obj, type = "prob")
      
      # Calculate metrics
      fold_metrics <- .calculate_metrics(y_test, pred, task, metrics)
      fold_results[[fold]] <- fold_metrics
    }
    
    all_results[[rep]] <- do.call(rbind, fold_results)
  }
  
  # Combine results
  combined <- do.call(rbind, all_results)
  
  result <- list(
    cv_metrics = combined,
    mean_metrics = colMeans(combined),
    sd_metrics = apply(combined, 2, sd),
    folds = folds,
    repeats = repeats,
    model = model,
    task = task
  )
  
  class(result) <- "ukb_ml_cv"
  
  if (verbose) {
    message("\nCross-validation results:")
    for (m in names(result$mean_metrics)) {
      message(sprintf("  %s: %.3f (SD: %.3f)", m, result$mean_metrics[m], result$sd_metrics[m]))
    }
  }
  
  result
}

#' @keywords internal
.create_folds <- function(n, k) {
  idx <- sample(n)
  split(idx, cut(seq_along(idx), k, labels = FALSE))
}

#' @keywords internal
.create_stratified_folds <- function(y, k) {
  folds <- vector("list", k)
  
  for (level in levels(y)) {
    level_idx <- which(y == level)
    level_folds <- .create_folds(length(level_idx), k)
    
    for (i in seq_len(k)) {
      folds[[i]] <- c(folds[[i]], level_idx[level_folds[[i]]])
    }
  }
  
  folds
}

#' @keywords internal
.calculate_metrics <- function(y_true, y_pred, task, metrics) {
  result <- numeric()
  
  if (task == "classification") {
    .check_ml_package("pROC")
    
    # Get probability for positive class
    if (is.matrix(y_pred)) {
      prob <- y_pred[, 2]
    } else {
      prob <- y_pred
    }
    
    # Binary predictions
    pred_class <- ifelse(prob > 0.5, levels(y_true)[2], levels(y_true)[1])
    pred_class <- factor(pred_class, levels = levels(y_true))
    
    # Confusion matrix elements
    tp <- sum(pred_class == levels(y_true)[2] & y_true == levels(y_true)[2])
    tn <- sum(pred_class == levels(y_true)[1] & y_true == levels(y_true)[1])
    fp <- sum(pred_class == levels(y_true)[2] & y_true == levels(y_true)[1])
    fn <- sum(pred_class == levels(y_true)[1] & y_true == levels(y_true)[2])
    
    for (m in metrics) {
      val <- switch(m,
        auc = as.numeric(auc(roc(y_true, prob, quiet = TRUE))),
        accuracy = (tp + tn) / (tp + tn + fp + fn),
        sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA,
        recall = if (tp + fn > 0) tp / (tp + fn) else NA,
        specificity = if (tn + fp > 0) tn / (tn + fp) else NA,
        ppv = if (tp + fp > 0) tp / (tp + fp) else NA,
        precision = if (tp + fp > 0) tp / (tp + fp) else NA,
        npv = if (tn + fn > 0) tn / (tn + fn) else NA,
        f1 = if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA,
        NA
      )
      result[m] <- val
    }
  } else {
    y_true <- as.numeric(y_true)
    y_pred <- as.numeric(y_pred)
    
    for (m in metrics) {
      val <- switch(m,
        rmse = sqrt(mean((y_true - y_pred)^2)),
        mae = mean(abs(y_true - y_pred)),
        rsquared = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2),
        NA
      )
      result[m] <- val
    }
  }
  
  result
}

# Model Comparison

#' Compare Multiple ML Models
#'
#' @description
#' Compare performance of multiple trained ML models.
#'
#' @param ... ukb_ml objects to compare
#' @param models Alternative: list of ukb_ml objects
#' @param metrics Metrics to compare
#' @param test_data Optional common test data
#' @param plot Whether to create comparison plot
#'
#' @return ukb_ml_compare object with comparison results
#'
#' @export
ukb_ml_compare <- function(...,
                           models = list(),
                           metrics = NULL,
                           test_data = NULL,
                           plot = TRUE) {
  .ukb_ml_deprecated("ukb_ml_compare", "multiple ukb_ml_workflow() calls compared by final_test_metrics")
  
  # Collect models
  dots <- list(...)
  if (length(dots) > 0) {
    models <- c(dots, models)
  }
  
  if (length(models) < 2) {
    stop("At least 2 models required for comparison")
  }
  
  # Get model names
  model_names <- sapply(models, function(m) {
    paste0(.get_model_label(m$model_type))
  })
  
  # Make unique names if duplicates
  if (any(duplicated(model_names))) {
    model_names <- make.unique(model_names)
  }
  
  task <- models[[1]]$task
  
  # Default metrics
  if (is.null(metrics)) {
    metrics <- if (task == "classification") {
      c("auc", "accuracy", "sensitivity", "specificity")
    } else {
      c("rmse", "mae", "rsquared")
    }
  }
  
  # Calculate metrics for each model
  comparison <- data.frame(model = model_names, stringsAsFactors = FALSE)
  
  for (i in seq_along(models)) {
    m <- models[[i]]
    
    if (!is.null(test_data)) {
      # Use common test data
      pred <- ukb_ml_predict(m, newdata = test_data, type = "prob")
      y_true <- test_data[[m$outcome]]
      if (task == "classification") y_true <- as.factor(y_true)
      met <- .calculate_metrics(y_true, pred, task, metrics)
    } else {
      met <- m$metrics[metrics]
    }
    
    for (metric in names(met)) {
      comparison[i, metric] <- met[metric]
    }
  }
  
  result <- list(
    comparison = comparison,
    models = model_names,
    task = task,
    metrics = metrics
  )
  
  class(result) <- "ukb_ml_compare"
  
  if (plot) {
    result$plot <- plot_ml_compare(result)
  }
  
  result
}

# S3 Methods

#' @export
print.ukb_ml <- function(x, ...) {
  cat("\n")
  cat("UKB Machine Learning Model\n")
  cat("\n")
  cat(sprintf("Model: %s\n", .get_model_label(x$model_type)))
  cat(sprintf("Task: %s\n", x$task))
  cat(sprintf("Outcome: %s\n", x$outcome))
  cat(sprintf("Predictors: %d variables\n", length(x$predictors)))
  cat(sprintf("Train size: %d\n", nrow(x$train_data)))
  cat(sprintf("Test size: %d\n", nrow(x$test_data)))
  cat("\n")
  
  if (!is.null(x$metrics)) {
    cat("Test Metrics:\n")
    for (m in names(x$metrics)) {
      cat(sprintf("  %s: %.3f\n", m, x$metrics[m]))
    }
  }
  
  invisible(x)
}

#' @export
summary.ukb_ml <- function(object, ...) {
  print(object)
  
  cat("\nVariable Importance (Top 10):\n")
  imp <- ukb_ml_importance(object)
  if (!is.null(imp)) {
    print(head(imp, 10))
  }
  
  invisible(object)
}

#' @export
predict.ukb_ml <- function(object, newdata = NULL, type = "response", ...) {
  ukb_ml_predict(object, newdata = newdata, type = type, ...)
}

#' @export
print.ukb_ml_cv <- function(x, ...) {
  cat("\n")
  cat("UKB ML Cross-Validation Results\n")
  cat("\n")
  cat(sprintf("Model: %s\n", .get_model_label(x$model)))
  cat(sprintf("Folds: %d, Repeats: %d\n", x$folds, x$repeats))
  cat("\nMetrics (Mean +/- SD):\n")
  for (m in names(x$mean_metrics)) {
    cat(sprintf("  %s: %.3f +/- %.3f\n", m, x$mean_metrics[m], x$sd_metrics[m]))
  }
  invisible(x)
}

#' @export
print.ukb_ml_compare <- function(x, ...) {
  cat("\n")
  cat("UKB ML Model Comparison\n")
  cat("\n")
  print(x$comparison, row.names = FALSE)
  invisible(x)
}
