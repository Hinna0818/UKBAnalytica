#' SHAP Explanations for Machine Learning Models
#'
#' @description
#' Compute and visualize SHAP (SHapley Additive exPlanations) values
#' for interpreting ML model predictions.
#'
#' @name ml_shap
#' @keywords internal
NULL

# Main SHAP Computation

#' Compute SHAP Values
#'
#' @description
#' Calculate SHAP values for model interpretation. SHAP values explain
#' each feature's contribution to individual predictions.
#'
#' @param object A \code{ukb_ml_workflow}, \code{ukb_ml_final}, or legacy
#'   \code{ukb_ml} object.
#' @param data Data for SHAP computation. If \code{object} is a
#'   \code{ukb_ml_workflow} and \code{data = NULL}, the frozen test set is used.
#'   If \code{object} is a \code{ukb_ml_final}, \code{data} is required.
#' @param nsim Number of Monte Carlo samples for SHAP estimation (default 100).
#'   Ignored when \code{method = "xgboost"}.
#' @param sample_n Optional; subsample observations for large datasets
#' @param seed Random seed
#' @param verbose Print progress
#' @param class_level Optional class to explain for multiclass
#'   \code{ukb_ml_workflow}/\code{ukb_ml_final} objects.
#' @param method SHAP backend. \code{"auto"} uses the native XGBoost
#'   contribution backend for XGBoost models and \code{fastshap} otherwise.
#' @param ... Additional arguments
#'
#' @return A ukb_shap object containing:
#' \itemize{
#'   \item shap_values: Matrix of SHAP values (n x p)
#'   \item baseline: Model baseline (expected) value
#'   \item feature_names: Names of features
#'   \item feature_values: Original feature values
#' }
#'
#' @examples
#' \dontrun{
#' # Recommended workflow object
#' ml <- ukb_ml_workflow(diabetes ~ age + bmi + sbp, data, model = "rf")
#'
#' # Compute SHAP values on the frozen test set
#' shap <- ukb_shap(ml, sample_n = 500)
#'
#' # Summary plot
#' ukb_shap_summary(shap)
#' }
#'
#' @export
ukb_shap <- function(object,
                     data = NULL,
                     nsim = 100,
                     sample_n = NULL,
                     seed = NULL,
                     verbose = TRUE,
                     class_level = NULL,
                     method = c("auto", "fastshap", "xgboost"),
                     ...) {
  method <- match.arg(method)

  if (inherits(object, "ukb_ml_workflow")) {
    if (is.null(object$final_model)) {
      stop("The workflow does not contain a final model. Run ukb_ml_workflow(..., fit_final = TRUE).", call. = FALSE)
    }
    if (is.null(data)) {
      data <- object$split$test
    }
    return(.ukb_shap_final_model(
      object = object$final_model,
      data = data,
      nsim = nsim,
      sample_n = sample_n,
      seed = seed,
      verbose = verbose,
      class_level = class_level,
      method = method,
      ...
    ))
  }

  if (inherits(object, "ukb_ml_final")) {
    if (is.null(data)) {
      stop("`data` is required when calling ukb_shap() on a ukb_ml_final object. You can also call ukb_shap() on the full ukb_ml_workflow object to use its frozen test set.", call. = FALSE)
    }
    return(.ukb_shap_final_model(
      object = object,
      data = data,
      nsim = nsim,
      sample_n = sample_n,
      seed = seed,
      verbose = verbose,
      class_level = class_level,
      method = method,
      ...
    ))
  }

  if (method == "xgboost") {
    stop("method = 'xgboost' is only available for ukb_ml_workflow or ukb_ml_final objects fitted with model = 'xgboost'.", call. = FALSE)
  }

  .check_ml_package("fastshap")
  
  if (!is.null(seed)) set.seed(seed)
  
  # Get data
  if (is.null(data)) {
    X <- object$X_test
    feature_values <- object$test_data[, object$predictors, drop = FALSE]
  } else {
    X <- data[, object$predictors, drop = FALSE]
    feature_values <- X
  }
  
  # Subsample if requested
  if (!is.null(sample_n) && sample_n < nrow(X)) {
    idx <- sample(nrow(X), sample_n)
    X <- X[idx, , drop = FALSE]
    feature_values <- feature_values[idx, , drop = FALSE]
    if (verbose) message(sprintf("Using %d sampled observations for SHAP", sample_n))
  }
  
  if (verbose) {
    message(sprintf("Computing SHAP values for %d observations...", nrow(X)))
  }
  
  # Create prediction wrapper
  pred_wrapper <- .create_shap_predict_wrapper(object)
  
  # Compute SHAP values
  shap_values <- fastshap::explain(
    object = object$model,
    feature_names = object$predictors,
    X = X,
    pred_wrapper = pred_wrapper,
    nsim = nsim,
    ...
  )
  
  # Calculate baseline (expected value)
  baseline <- .calculate_baseline(object)
  
  result <- list(
    shap_values = as.matrix(shap_values),
    baseline = baseline,
    feature_names = object$predictors,
    feature_values = feature_values,
    model_type = object$model_type,
    task = object$task
  )
  
  class(result) <- "ukb_shap"
  
  if (verbose) message("SHAP computation complete")
  
  result
}

#' Compute SHAP for the new ukb_ml_workflow final model object
#' @keywords internal
#' @noRd
.ukb_shap_final_model <- function(object,
                                  data,
                                  nsim = 100,
                                  sample_n = NULL,
                                  seed = NULL,
                                  verbose = TRUE,
                                  class_level = NULL,
                                  method = c("auto", "fastshap", "xgboost"),
                                  ...) {
  if (!inherits(object, "ukb_ml_final")) {
    stop("`object` must be a ukb_ml_final object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)

  if ((method == "auto" && identical(object$model, "xgboost")) || method == "xgboost") {
    if (!identical(object$model, "xgboost")) {
      stop("method = 'xgboost' requires a ukb_ml_final object fitted with model = 'xgboost'.", call. = FALSE)
    }
    return(.ukb_shap_xgboost_final_model(
      object = object,
      data = data,
      sample_n = sample_n,
      seed = seed,
      verbose = verbose
    ))
  }

  .check_ml_package("fastshap")

  features <- object$selected_features %||% object$predictors
  data <- as.data.frame(data)
  missing_features <- setdiff(features, names(data))
  if (length(missing_features) > 0L) {
    stop("SHAP data is missing feature column(s): ", paste(missing_features, collapse = ", "), call. = FALSE)
  }

  X <- data[, features, drop = FALSE]
  feature_values <- X

  if (!is.null(sample_n) && sample_n < nrow(X)) {
    idx <- sample(nrow(X), sample_n)
    X <- X[idx, , drop = FALSE]
    feature_values <- feature_values[idx, , drop = FALSE]
    if (isTRUE(verbose)) message(sprintf("Using %d sampled observations for SHAP", sample_n))
  }

  if (object$outcome_type == "multiclass") {
    if (is.null(class_level)) {
      stop("`class_level` is required for multiclass SHAP explanations.", call. = FALSE)
    }
    if (!class_level %in% object$classes) {
      stop("`class_level` must be one of: ", paste(object$classes, collapse = ", "), call. = FALSE)
    }
  }

  pred_wrapper <- function(model, newdata) {
    if (model$outcome_type == "continuous") {
      return(as.numeric(.ukb_ml_predict_core(model, newdata, type = "response")))
    }
    if (model$outcome_type == "binary") {
      return(as.numeric(.ukb_ml_predict_core(model, newdata, type = "prob")))
    }
    prob <- .ukb_ml_predict_core(model, newdata, type = "prob")
    as.numeric(prob[, class_level])
  }

  if (isTRUE(verbose)) {
    target <- if (object$outcome_type == "multiclass") paste0("class=", class_level) else object$outcome_type
    message(sprintf("Computing SHAP values for %d observations (%s)...", nrow(X), target))
  }

  shap_values <- fastshap::explain(
    object = object,
    feature_names = features,
    X = X,
    pred_wrapper = pred_wrapper,
    nsim = nsim,
    ...
  )

  baseline <- mean(pred_wrapper(object, X), na.rm = TRUE)

  result <- list(
    shap_values = as.matrix(shap_values),
    baseline = baseline,
    feature_names = features,
    feature_values = feature_values,
    model_type = object$model,
    task = if (object$outcome_type == "continuous") "regression" else "classification",
    outcome_type = object$outcome_type,
    class_level = class_level
  )

  class(result) <- "ukb_shap"
  if (isTRUE(verbose)) message("SHAP computation complete")
  result
}

#' Compute native XGBoost SHAP values for a final model object
#' @keywords internal
#' @noRd
.ukb_shap_xgboost_final_model <- function(object,
                                          data,
                                          sample_n = NULL,
                                          seed = NULL,
                                          verbose = TRUE) {
  .check_ml_package("xgboost")

  if (!inherits(object, "ukb_ml_final")) {
    stop("`object` must be a ukb_ml_final object.", call. = FALSE)
  }
  if (!identical(object$model, "xgboost")) {
    stop("Native XGBoost SHAP requires model = 'xgboost'.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  data <- as.data.frame(data)
  features <- object$selected_features %||% object$predictors
  missing_features <- setdiff(features, names(data))
  if (length(missing_features) > 0L) {
    stop("SHAP data is missing feature column(s): ", paste(missing_features, collapse = ", "), call. = FALSE)
  }

  if (!is.null(sample_n) && sample_n < nrow(data)) {
    idx <- sample(nrow(data), sample_n)
    data <- data[idx, , drop = FALSE]
    if (isTRUE(verbose)) message(sprintf("Using %d sampled observations for SHAP", sample_n))
  }

  mf <- stats::model.frame(
    object$x_terms,
    data = data,
    na.action = stats::na.pass,
    xlev = object$xlevels
  )
  X <- stats::model.matrix(object$x_terms, data = mf, contrasts.arg = object$contrasts)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  expected_cols <- object$feature_names
  missing_matrix_cols <- setdiff(expected_cols, colnames(X))
  if (length(missing_matrix_cols) > 0L) {
    add <- matrix(0, nrow = nrow(X), ncol = length(missing_matrix_cols))
    colnames(add) <- missing_matrix_cols
    X <- cbind(X, add)
  }
  extra_matrix_cols <- setdiff(colnames(X), expected_cols)
  if (length(extra_matrix_cols) > 0L) {
    X <- X[, setdiff(colnames(X), extra_matrix_cols), drop = FALSE]
  }
  X <- X[, expected_cols, drop = FALSE]
  storage.mode(X) <- "numeric"

  if (isTRUE(verbose)) {
    message(sprintf("Computing native XGBoost SHAP values for %d observations...", nrow(X)))
  }
  contrib <- predict(
    object$fitted_model,
    xgboost::xgb.DMatrix(X),
    predcontrib = TRUE
  )
  contrib <- as.matrix(contrib)
  bias_col <- ncol(contrib)
  shap_values <- contrib[, -bias_col, drop = FALSE]
  colnames(shap_values) <- expected_cols

  result <- list(
    shap_values = shap_values,
    baseline = mean(contrib[, bias_col], na.rm = TRUE),
    feature_names = expected_cols,
    feature_values = as.data.frame(X, check.names = FALSE),
    model_type = object$model,
    task = if (object$outcome_type == "continuous") "regression" else "classification",
    outcome_type = object$outcome_type,
    class_level = if (object$outcome_type == "binary") object$positive_class else NULL,
    method = "xgboost"
  )

  class(result) <- "ukb_shap"
  if (isTRUE(verbose)) message("SHAP computation complete")
  result
}

#' Create prediction wrapper for SHAP
#' @keywords internal
.create_shap_predict_wrapper <- function(object) {
  model_type <- object$model_type
  task <- object$task
  
  switch(model_type,
    rf = function(model, newdata) {
      pred <- predict(model, data = newdata)
      if (task == "classification") {
        pred$predictions[, 2]  # Probability of positive class
      } else {
        pred$predictions
      }
    },
    xgboost = function(model, newdata) {
      .check_ml_package("xgboost")
      X_mat <- as.matrix(newdata)
      for (i in seq_len(ncol(X_mat))) {
        if (is.factor(newdata[, i])) {
          X_mat[, i] <- as.numeric(newdata[, i]) - 1
        }
      }
      mode(X_mat) <- "numeric"
      predict(model, X_mat)
    },
    glmnet = function(model, newdata) {
      .check_ml_package("glmnet")
      X_mat <- model.matrix(~ . - 1, data = newdata)
      as.numeric(predict(model, newx = X_mat, s = "lambda.min", type = "response"))
    },
    logistic = function(model, newdata) {
      predict(model, newdata = newdata, type = "response")
    },
    {
      # Default wrapper
      function(model, newdata) {
        pred <- predict(model, newdata = newdata)
        if (is.matrix(pred)) pred[, 1] else pred
      }
    }
  )
}

#' Calculate model baseline
#' @keywords internal
.calculate_baseline <- function(object) {
  if (object$task == "classification") {
    # Baseline is the average predicted probability
    pred <- ukb_ml_predict(object, newdata = object$train_data, type = "prob")
    if (is.matrix(pred)) mean(pred[, 2]) else mean(pred)
  } else {
    mean(object$y_train)
  }
}

# SHAP Summary

#' SHAP Summary Statistics
#'
#' @description
#' Calculate summary statistics from SHAP values.
#'
#' @param object A ukb_shap object
#' @param n Number of top features to show (default 20)
#' @param ... Additional arguments
#'
#' @return Data frame with feature importance based on SHAP
#'
#' @export
ukb_shap_summary <- function(object, n = 20, ...) {
  
  shap <- object$shap_values
  
  # Calculate mean absolute SHAP
  importance <- data.frame(
    feature = object$feature_names,
    mean_abs_shap = colMeans(abs(shap)),
    mean_shap = colMeans(shap),
    sd_shap = apply(shap, 2, sd),
    stringsAsFactors = FALSE
  )
  
  importance <- importance[order(importance$mean_abs_shap, decreasing = TRUE), ]
  rownames(importance) <- NULL
  
  if (!is.null(n) && n < nrow(importance)) {
    importance <- importance[seq_len(n), ]
  }
  
  importance
}

# SHAP Dependence

#' SHAP Dependence Values
#'
#' @description
#' Get SHAP dependence data for a specific feature.
#'
#' @param object A ukb_shap object
#' @param feature Feature name to analyze
#' @param color_feature Optional feature for coloring (interaction analysis)
#' @param ... Additional arguments
#'
#' @return Data frame with feature values and SHAP values
#'
#' @export
ukb_shap_dependence <- function(object, feature, color_feature = NULL, ...) {
  
  if (!feature %in% object$feature_names) {
    stop(sprintf("Feature '%s' not found", feature))
  }
  
  feature_idx <- which(object$feature_names == feature)
  
  dep_data <- data.frame(
    feature_value = object$feature_values[[feature]],
    shap_value = object$shap_values[, feature_idx]
  )
  
  if (!is.null(color_feature)) {
    if (!color_feature %in% object$feature_names) {
      stop(sprintf("Color feature '%s' not found", color_feature))
    }
    dep_data$color_value <- object$feature_values[[color_feature]]
  }
  
  attr(dep_data, "feature") <- feature
  attr(dep_data, "color_feature") <- color_feature
  
  dep_data
}

# SHAP Force Plot Data
#' SHAP Force Plot Data
#'
#' @description
#' Get SHAP contribution data for a single observation (force plot).
#'
#' @param object A ukb_shap object
#' @param row_id Row index to explain
#' @param max_features Maximum features to show
#' @param ... Additional arguments
#'
#' @return Data frame with feature contributions for the observation
#'
#' @export
ukb_shap_force <- function(object, row_id = 1, max_features = 10, ...) {
  
  if (row_id > nrow(object$shap_values)) {
    stop("row_id exceeds number of observations")
  }
  
  shap_row <- object$shap_values[row_id, ]
  feature_row <- object$feature_values[row_id, ]
  
  force_data <- data.frame(
    feature = object$feature_names,
    value = as.character(unlist(feature_row)),
    shap = shap_row,
    stringsAsFactors = FALSE
  )
  
  # Sort by absolute contribution
  force_data <- force_data[order(abs(force_data$shap), decreasing = TRUE), ]
  
  if (!is.null(max_features) && max_features < nrow(force_data)) {
    force_data <- force_data[seq_len(max_features), ]
  }
  
  attr(force_data, "baseline") <- object$baseline
  attr(force_data, "prediction") <- object$baseline + sum(shap_row)
  
  force_data
}

# S3 Methods

#' @export
print.ukb_shap <- function(x, ...) {
  cat("\n")
  cat("SHAP Explanation Object\n")
  cat("\n")
  cat(sprintf("Observations: %d\n", nrow(x$shap_values)))
  cat(sprintf("Features: %d\n", ncol(x$shap_values)))
  cat(sprintf("Baseline: %.4f\n", x$baseline))
  cat("\nTop features by mean |SHAP|:\n")
  
  summary <- ukb_shap_summary(x, n = 10)
  print(summary[, c("feature", "mean_abs_shap")], row.names = FALSE)
  
  invisible(x)
}

#' @export
summary.ukb_shap <- function(object, n = 20, ...) {
  ukb_shap_summary(object, n = n, ...)
}
