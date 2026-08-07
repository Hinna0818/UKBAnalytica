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
#' Calculate SHAP values through native TreeSHAP, KernelSHAP, or FastSHAP while
#' preserving one output structure across supported machine-learning models.
#'
#' @param object A \code{ukb_ml_workflow}, \code{ukb_ml_final}, or legacy
#'   \code{ukb_ml} object.
#' @param data Data to explain. A workflow uses its frozen test set when this is
#'   \code{NULL}. A final-model object requires explicit data.
#' @param nsim Number of Monte Carlo samples for FastSHAP. Ignored by TreeSHAP
#'   and KernelSHAP.
#' @param sample_n Optional number of observations to explain.
#' @param seed Optional random seed.
#' @param verbose Print progress messages.
#' @param class_level Class probability to explain for multiclass models.
#' @param method SHAP backend. \code{"auto"} uses TreeSHAP for XGBoost and an
#'   installed model-agnostic backend otherwise. Legacy values
#'   \code{"xgboost"} and \code{"permutation"} are aliases for
#'   \code{"treeshap"} and \code{"fastshap"}. \code{"kernalshap"} is accepted
#'   as an alias for \code{"kernelshap"}.
#' @param background_data Optional reference data for KernelSHAP or FastSHAP.
#'   A workflow uses its training split by default. A final model uses
#'   \code{data} when no background is supplied.
#' @param background_n Maximum number of background observations. Use
#'   \code{NULL} to retain all supplied background observations.
#' @param ... Additional arguments passed to the selected backend.
#'
#' @return A \code{ukb_shap} object containing the SHAP matrix, baseline,
#'   feature names and values, backend name, output scale, background sample
#'   size, and a local-accuracy diagnostic where available.
#'
#' @export
ukb_shap <- function(object,
                     data = NULL,
                     nsim = 100,
                     sample_n = NULL,
                     seed = NULL,
                     verbose = TRUE,
                     class_level = NULL,
                     method = c(
                       "auto", "treeshap", "kernelshap", "fastshap",
                       "xgboost", "permutation", "kernalshap"
                     ),
                     background_data = NULL,
                     background_n = 200,
                     ...) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  if (inherits(object, "ukb_ml_workflow")) {
    if (is.null(object$final_model)) {
      stop("The workflow does not contain a final model. Run ukb_ml_workflow(..., fit_final = TRUE).", call. = FALSE)
    }
    if (is.null(data)) {
      data <- object$split$test
    }
    if (is.null(background_data)) {
      background_data <- object$split$train
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
      background_data = background_data,
      background_n = background_n,
      ...
    ))
  }

  if (inherits(object, "ukb_ml_final")) {
    if (is.null(data)) {
      stop("`data` is required when calling ukb_shap() on a ukb_ml_final object.", call. = FALSE)
    }
    if (is.null(background_data)) {
      background_data <- data
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
      background_data = background_data,
      background_n = background_n,
      ...
    ))
  }

  if (!inherits(object, "ukb_ml")) {
    stop("`object` must be a ukb_ml_workflow, ukb_ml_final, or legacy ukb_ml object.", call. = FALSE)
  }

  if (is.null(data)) {
    X <- object$X_test
    feature_values <- object$test_data[, object$predictors, drop = FALSE]
  } else {
    X <- data[, object$predictors, drop = FALSE]
    feature_values <- X
  }
  if (!is.null(sample_n) && sample_n < nrow(X)) {
    idx <- sample(nrow(X), sample_n)
    X <- X[idx, , drop = FALSE]
    feature_values <- feature_values[idx, , drop = FALSE]
    if (isTRUE(verbose)) message(sprintf("Using %d sampled observations for SHAP", sample_n))
  }

  resolved_method <- .ukb_shap_resolve_method(method, object$model_type)
  if (identical(resolved_method, "treeshap")) {
    return(.ukb_shap_xgboost_legacy(
      object = object,
      X = X,
      feature_values = feature_values,
      class_level = class_level,
      verbose = verbose
    ))
  }

  if (is.null(background_data)) {
    background_data <- object$train_data
  }
  background_X <- background_data[, object$predictors, drop = FALSE]
  background_X <- .ukb_shap_sample_background(background_X, background_n)
  pred_wrapper <- .create_shap_predict_wrapper(object)
  backend <- .ukb_shap_model_agnostic(
    method = resolved_method,
    object = object$model,
    X = X,
    background_X = background_X,
    pred_wrapper = pred_wrapper,
    nsim = nsim,
    verbose = verbose,
    ...
  )

  result <- list(
    shap_values = backend$shap_values,
    baseline = backend$baseline,
    feature_names = object$predictors,
    feature_values = feature_values,
    model_type = object$model_type,
    task = object$task,
    method = resolved_method,
    output_scale = if (object$task == "classification") "probability" else "response",
    background_n = nrow(background_X),
    local_accuracy_error = backend$local_accuracy_error
  )
  class(result) <- "ukb_shap"
  if (isTRUE(verbose)) message("SHAP computation complete")
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
                                  method = c(
                                    "auto", "treeshap", "kernelshap", "fastshap",
                                    "xgboost", "permutation", "kernalshap"
                                  ),
                                  background_data = NULL,
                                  background_n = 200,
                                  ...) {
  if (!inherits(object, "ukb_ml_final")) {
    stop("`object` must be a ukb_ml_final object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)
  resolved_method <- .ukb_shap_resolve_method(method, object$model)

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

  if (identical(resolved_method, "treeshap")) {
    return(.ukb_shap_xgboost_final_model(
      object = object,
      data = data,
      sample_n = sample_n,
      seed = seed,
      verbose = verbose,
      class_level = class_level
    ))
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

  if (is.null(background_data)) {
    background_data <- data
  }
  background_data <- as.data.frame(background_data)
  missing_background <- setdiff(features, names(background_data))
  if (length(missing_background) > 0L) {
    stop(
      "SHAP background data is missing feature column(s): ",
      paste(missing_background, collapse = ", "),
      call. = FALSE
    )
  }
  background_X <- background_data[, features, drop = FALSE]
  background_X <- .ukb_shap_sample_background(background_X, background_n)

  backend <- .ukb_shap_model_agnostic(
    method = resolved_method,
    object = object,
    X = X,
    background_X = background_X,
    pred_wrapper = pred_wrapper,
    nsim = nsim,
    verbose = verbose,
    ...
  )

  result <- list(
    shap_values = backend$shap_values,
    baseline = backend$baseline,
    feature_names = features,
    feature_values = feature_values,
    model_type = object$model,
    task = if (object$outcome_type == "continuous") "regression" else "classification",
    outcome_type = object$outcome_type,
    class_level = class_level,
    method = resolved_method,
    output_scale = if (object$outcome_type == "continuous") "response" else "probability",
    background_n = nrow(background_X),
    local_accuracy_error = backend$local_accuracy_error
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
                                          verbose = TRUE,
                                          class_level = NULL) {
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

  mf <- model.frame(
    object$x_terms,
    data = data,
    na.action = na.pass,
    xlev = object$xlevels
  )
  X <- model.matrix(object$x_terms, data = mf, contrasts.arg = object$contrasts)
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
    message(sprintf("Computing native TreeSHAP values for %d observations...", nrow(X)))
  }

  class_index <- 1L
  explained_class <- NULL
  if (object$outcome_type == "binary") {
    explained_class <- object$positive_class
  } else if (object$outcome_type == "multiclass") {
    if (is.null(class_level) || !class_level %in% object$classes) {
      stop("`class_level` must identify one class for multiclass TreeSHAP.", call. = FALSE)
    }
    class_index <- match(class_level, object$classes)
    explained_class <- class_level
  }

  contrib <- .ukb_xgboost_contribution_matrix(
    model = object$fitted_model,
    X = X,
    feature_names = expected_cols,
    class_index = class_index
  )
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
    class_level = explained_class,
    method = "treeshap",
    output_scale = if (object$outcome_type == "continuous") "response" else "margin",
    background_n = object$train_rows,
    local_accuracy_error = attr(contrib, "local_accuracy_error")
  )

  class(result) <- "ukb_shap"
  if (isTRUE(verbose)) message("SHAP computation complete")
  result
}

#' Compute native TreeSHAP for a legacy ukb_ml object
#' @keywords internal
#' @noRd
.ukb_shap_xgboost_legacy <- function(object,
                                     X,
                                     feature_values,
                                     class_level = NULL,
                                     verbose = TRUE) {
  if (!identical(object$model_type, "xgboost")) {
    stop("TreeSHAP requires an XGBoost model.", call. = FALSE)
  }
  X_mat <- .ukb_shap_numeric_matrix(X)
  classes <- if (is.factor(object$y_train)) levels(object$y_train) else unique(object$y_train)
  class_index <- 1L
  explained_class <- NULL
  if (identical(object$task, "classification") && length(classes) > 2L) {
    if (is.null(class_level) || !class_level %in% classes) {
      stop("`class_level` must identify one class for multiclass TreeSHAP.", call. = FALSE)
    }
    class_index <- match(class_level, classes)
    explained_class <- class_level
  } else if (identical(object$task, "classification")) {
    explained_class <- as.character(classes[[2]])
  }
  contrib <- .ukb_xgboost_contribution_matrix(
    model = object$model,
    X = X_mat,
    feature_names = object$predictors,
    class_index = class_index
  )
  bias_col <- ncol(contrib)
  result <- list(
    shap_values = contrib[, -bias_col, drop = FALSE],
    baseline = mean(contrib[, bias_col], na.rm = TRUE),
    feature_names = object$predictors,
    feature_values = feature_values,
    model_type = object$model_type,
    task = object$task,
    outcome_type = if (object$task == "regression") "continuous" else if (length(classes) > 2L) "multiclass" else "binary",
    class_level = explained_class,
    method = "treeshap",
    output_scale = if (object$task == "regression") "response" else "margin",
    background_n = nrow(object$X_train),
    local_accuracy_error = attr(contrib, "local_accuracy_error")
  )
  class(result) <- "ukb_shap"
  if (isTRUE(verbose)) message("SHAP computation complete")
  result
}

#' Resolve SHAP backend names and legacy aliases
#' @keywords internal
#' @noRd
.ukb_shap_resolve_method <- function(method, model_type) {
  aliases <- c(
    xgboost = "treeshap",
    permutation = "fastshap",
    kernalshap = "kernelshap"
  )
  if (method %in% names(aliases)) {
    method <- unname(aliases[[method]])
  }
  if (identical(method, "auto")) {
    if (identical(model_type, "xgboost")) {
      return("treeshap")
    }
    if (requireNamespace("kernelshap", quietly = TRUE)) {
      return("kernelshap")
    }
    if (requireNamespace("fastshap", quietly = TRUE)) {
      return("fastshap")
    }
    stop(
      "Model-agnostic SHAP requires package 'kernelshap' or 'fastshap'. ",
      "Install one of these packages or use an XGBoost model with TreeSHAP.",
      call. = FALSE
    )
  }
  if (identical(method, "treeshap") && !identical(model_type, "xgboost")) {
    stop("TreeSHAP is only available for XGBoost models.", call. = FALSE)
  }
  method
}

#' Sample a reproducible SHAP background set
#' @keywords internal
#' @noRd
.ukb_shap_sample_background <- function(background_X, background_n = 200) {
  background_X <- as.data.frame(background_X)
  if (nrow(background_X) == 0L) {
    stop("SHAP background data must contain at least one row.", call. = FALSE)
  }
  if (is.null(background_n)) {
    return(background_X)
  }
  if (!is.numeric(background_n) || length(background_n) != 1L ||
      is.na(background_n) || background_n < 1) {
    stop("`background_n` must be NULL or a positive integer.", call. = FALSE)
  }
  background_n <- as.integer(background_n)
  if (nrow(background_X) > background_n) {
    background_X <- background_X[sample.int(nrow(background_X), background_n), , drop = FALSE]
  }
  background_X
}

#' Run a model-agnostic SHAP backend
#' @keywords internal
#' @noRd
.ukb_shap_model_agnostic <- function(method,
                                     object,
                                     X,
                                     background_X,
                                     pred_wrapper,
                                     nsim = 100,
                                     verbose = TRUE,
                                     ...) {
  X <- as.data.frame(X)
  background_X <- as.data.frame(background_X)
  if (!identical(names(X), names(background_X))) {
    stop("SHAP explanation and background data must have identical feature columns.", call. = FALSE)
  }

  extra_args <- list(...)
  if (length(extra_args) > 0L && (is.null(names(extra_args)) || any(!nzchar(names(extra_args))))) {
    stop("Additional SHAP backend arguments in `...` must be named.", call. = FALSE)
  }

  if (identical(method, "kernelshap")) {
    if (!requireNamespace("kernelshap", quietly = TRUE)) {
      stop("Package 'kernelshap' is required for method = 'kernelshap'.", call. = FALSE)
    }
    kernelshap_fun <- getExportedValue("kernelshap", "kernelshap")
    args <- modifyList(
      list(
        object = object,
        X = X,
        bg_X = background_X,
        pred_fun = pred_wrapper,
        verbose = verbose
      ),
      extra_args
    )
    fit <- do.call(kernelshap_fun, args)
    shap_values <- as.matrix(fit$S)
    baseline <- as.numeric(fit$baseline)[1]
  } else if (identical(method, "fastshap")) {
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("Package 'fastshap' is required for method = 'fastshap'.", call. = FALSE)
    }
    nsim <- as.integer(nsim)
    if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
      stop("`nsim` must be a positive integer for FastSHAP.", call. = FALSE)
    }
    fastshap_fun <- getExportedValue("fastshap", "explain")
    args <- modifyList(
      list(
        object = object,
        feature_names = names(X),
        X = background_X,
        newdata = X,
        pred_wrapper = pred_wrapper,
        nsim = nsim,
        adjust = nsim > 1L,
        shap_only = FALSE,
        parallel = FALSE
      ),
      extra_args
    )
    fit <- do.call(fastshap_fun, args)
    shap_values <- as.matrix(fit$shapley_values)
    baseline <- as.numeric(fit$baseline)[1]
  } else {
    stop("Unsupported model-agnostic SHAP backend: ", method, call. = FALSE)
  }

  if (!identical(dim(shap_values), c(nrow(X), ncol(X)))) {
    stop("The SHAP backend returned an unexpected matrix shape.", call. = FALSE)
  }
  colnames(shap_values) <- names(X)
  pred <- as.numeric(pred_wrapper(object, X))
  local_error <- if (length(pred) == nrow(X)) {
    max(abs(baseline + rowSums(shap_values) - pred), na.rm = TRUE)
  } else {
    NA_real_
  }
  if (!is.finite(local_error)) local_error <- NA_real_

  list(
    shap_values = shap_values,
    baseline = baseline,
    local_accuracy_error = local_error
  )
}

#' Convert legacy predictors to the numeric matrix used by XGBoost
#' @keywords internal
#' @noRd
.ukb_shap_numeric_matrix <- function(X) {
  X <- as.data.frame(X)
  if (any(vapply(X, is.character, logical(1)))) {
    stop("Character predictors are not supported by legacy XGBoost SHAP.", call. = FALSE)
  }
  X_mat <- as.matrix(X)
  for (i in seq_along(X)) {
    if (is.factor(X[[i]])) {
      X_mat[, i] <- as.numeric(X[[i]]) - 1
    }
  }
  storage.mode(X_mat) <- "numeric"
  X_mat
}

#' Return one class-specific XGBoost contribution matrix
#' @keywords internal
#' @noRd
.ukb_xgboost_contribution_matrix <- function(model,
                                              X,
                                              feature_names,
                                              class_index = 1L) {
  .check_ml_package("xgboost")
  X <- as.matrix(X)
  dmat <- xgboost::xgb.DMatrix(X)
  contrib_raw <- tryCatch(
    predict(model, dmat, predcontrib = TRUE, strict_shape = TRUE),
    error = function(e) predict(model, dmat, predcontrib = TRUE)
  )
  dims <- dim(contrib_raw)
  expected_columns <- length(feature_names) + 1L

  if (length(dims) == 3L) {
    if (class_index < 1L || class_index > dims[[2]]) {
      stop("Requested TreeSHAP class index is outside the model output.", call. = FALSE)
    }
    if (dims[[1]] == nrow(X) && dims[[3]] == expected_columns) {
      selected <- contrib_raw[, class_index, , drop = FALSE]
      contrib <- matrix(selected, nrow = nrow(X), ncol = expected_columns)
    } else if (dims[[1]] == expected_columns && dims[[3]] == nrow(X)) {
      selected <- contrib_raw[, class_index, , drop = FALSE]
      contrib <- t(matrix(selected, nrow = expected_columns, ncol = nrow(X)))
    } else {
      stop("Unexpected multidimensional TreeSHAP output shape.", call. = FALSE)
    }
  } else {
    contrib <- as.matrix(contrib_raw)
    if (nrow(contrib) != nrow(X) || ncol(contrib) != expected_columns) {
      stop("Unexpected TreeSHAP output shape.", call. = FALSE)
    }
  }
  colnames(contrib) <- c(feature_names, "BIAS")

  margin <- predict(model, dmat, outputmargin = TRUE)
  if (is.matrix(margin)) {
    margin <- margin[, class_index]
  } else if (length(margin) != nrow(X) && length(margin) %% nrow(X) == 0L) {
    margin <- matrix(margin, nrow = nrow(X), byrow = TRUE)[, class_index]
  }
  local_error <- if (length(margin) == nrow(X)) {
    max(abs(rowSums(contrib) - as.numeric(margin)), na.rm = TRUE)
  } else {
    NA_real_
  }
  attr(contrib, "local_accuracy_error") <- local_error
  contrib
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
