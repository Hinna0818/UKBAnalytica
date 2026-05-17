#' Machine Learning Model Evaluation
#'
#' @description
#' Functions for evaluating and visualizing ML model performance.
#'
#' @name ml_evaluate
#' @keywords internal
#' @importFrom pROC roc auc ci.auc
NULL

#' Create ROC Curve Data for Binary ML Predictions
#'
#' @description
#' Converts binary outcome predictions into a tidy ROC curve table with AUC and
#' optional 95% confidence interval. This helper is useful for plotting one or
#' more model ROC curves without re-running model evaluation.
#'
#' @param truth True binary outcome values.
#' @param prob Predicted probability for the positive class.
#' @param model_id Optional model identifier.
#' @param model_label Optional model label used in plots.
#' @param positive_class Optional positive class label. Defaults to the second
#'   level after converting \code{truth} to a factor.
#' @param ci Logical. Calculate AUC 95% confidence interval.
#' @param ci_method Method passed to \code{ci.auc()}; usually
#'   \code{"delong"} or \code{"bootstrap"}.
#' @param quiet Logical passed to \code{roc()}.
#'
#' @return A data.frame with specificity, sensitivity, false-positive rate,
#'   threshold, AUC, and optional confidence interval columns.
#'
#' @examples
#' \dontrun{
#' roc_df <- ukb_ml_roc_data(
#'   truth = test$case,
#'   prob = pred_prob,
#'   model_label = "XGBoost"
#' )
#' plot_ml_roc_compare(roc_df)
#' }
#'
#' @export
ukb_ml_roc_data <- function(truth,
                            prob,
                            model_id = NULL,
                            model_label = NULL,
                            positive_class = NULL,
                            ci = TRUE,
                            ci_method = c("delong", "bootstrap"),
                            quiet = TRUE) {
  .check_ml_package("pROC")

  ci_method <- match.arg(ci_method)
  truth <- droplevels(as.factor(truth))
  if (nlevels(truth) != 2L) {
    stop("ukb_ml_roc_data() only supports binary outcomes.", call. = FALSE)
  }

  if (is.null(positive_class)) {
    positive_class <- levels(truth)[2]
  }
  if (!positive_class %in% levels(truth)) {
    stop("`positive_class` must be one of: ", paste(levels(truth), collapse = ", "), call. = FALSE)
  }
  negative_class <- setdiff(levels(truth), positive_class)[1]

  prob <- as.numeric(prob)
  keep <- !is.na(truth) & !is.na(prob)
  truth <- droplevels(truth[keep])
  prob <- prob[keep]
  if (length(unique(truth)) != 2L) {
    stop("Both outcome classes must be present after removing missing values.", call. = FALSE)
  }

  roc_obj <- roc(
    response = truth,
    predictor = prob,
    levels = c(negative_class, positive_class),
    quiet = quiet
  )
  auc_value <- as.numeric(auc(roc_obj))

  if (isTRUE(ci)) {
    auc_ci <- as.numeric(ci.auc(roc_obj, method = ci_method))
    lower95 <- auc_ci[[1]]
    upper95 <- auc_ci[[3]]
  } else {
    lower95 <- NA_real_
    upper95 <- NA_real_
  }

  if (is.null(model_id)) {
    model_id <- "model"
  }
  if (is.null(model_label)) {
    model_label <- model_id
  }

  out <- data.frame(
    model_id = model_id,
    model_label = model_label,
    specificity = as.numeric(roc_obj$specificities),
    sensitivity = as.numeric(roc_obj$sensitivities),
    fpr = 1 - as.numeric(roc_obj$specificities),
    threshold = as.numeric(roc_obj$thresholds),
    auc = auc_value,
    lower95 = lower95,
    upper95 = upper95,
    ci_method = if (isTRUE(ci)) ci_method else NA_character_,
    positive_class = positive_class,
    n = length(prob),
    n_event = sum(truth == positive_class),
    stringsAsFactors = FALSE
  )
  class(out) <- c("ukb_ml_roc_data", class(out))
  out
}

# Metrics Function

#' Calculate Model Performance Metrics
#'
#' @description
#' Compute performance metrics for a trained ML model.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data for evaluation
#' @param metrics Specific metrics to compute (NULL for defaults)
#' @param ci Logical; compute confidence intervals (default FALSE)
#' @param ci_method Method for CI: "bootstrap" or "delong" (for AUC)
#' @param n_boot Number of bootstrap samples
#' @param verbose Print results
#' @param ... Additional arguments
#'
#' @return Named vector or list with metrics and optional CIs
#'
#' @export
ukb_ml_metrics <- function(object,
                           newdata = NULL,
                           metrics = NULL,
                           ci = FALSE,
                           ci_method = c("bootstrap", "delong"),
                           n_boot = 1000,
                           verbose = TRUE,
                           ...) {
  .ukb_ml_deprecated("ukb_ml_metrics", "ukb_ml_evaluate_test() or ml$final_test_metrics")
  
  ci_method <- match.arg(ci_method)
  task <- object$task
  
  # Get predictions and true values
  if (is.null(newdata)) {
    y_true <- object$y_test
    pred <- ukb_ml_predict(object, type = "prob")
  } else {
    y_true <- newdata[[object$outcome]]
    if (task == "classification") y_true <- as.factor(y_true)
    pred <- ukb_ml_predict(object, newdata = newdata, type = "prob")
  }
  
  # Default metrics
  if (is.null(metrics)) {
    metrics <- if (task == "classification") {
      c("auc", "accuracy", "sensitivity", "specificity", "ppv", "npv", "f1", "brier")
    } else {
      c("rmse", "mae", "rsquared")
    }
  }
  
  # Calculate metrics
  result <- .calculate_all_metrics(y_true, pred, task, metrics)
  
  # Add confidence intervals if requested
  if (ci) {
    if (task == "classification" && ci_method == "delong" && "auc" %in% metrics) {
      .check_ml_package("pROC")
      prob <- if (is.matrix(pred)) pred[, 2] else pred
      roc_obj <- roc(y_true, prob, quiet = TRUE)
      auc_ci <- as.numeric(ci.auc(roc_obj))
      result <- list(
        metrics = result,
        ci = list(auc = c(lower = auc_ci[1], estimate = auc_ci[2], upper = auc_ci[3]))
      )
    } else {
      # Bootstrap CI
      boot_results <- .bootstrap_metrics(y_true, pred, task, metrics, n_boot)
      result <- list(
        metrics = result,
        ci = boot_results
      )
    }
  }
  
  if (verbose && !ci) {
    cat("\nModel Performance Metrics:\n")
    for (m in names(result)) {
      cat(sprintf("  %s: %.4f\n", m, result[m]))
    }
  }
  
  result
}

#' @keywords internal
.calculate_all_metrics <- function(y_true, y_pred, task, metrics) {
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
        recall      = if (tp + fn > 0) tp / (tp + fn) else NA,
        specificity = if (tn + fp > 0) tn / (tn + fp) else NA,
        ppv       = if (tp + fp > 0) tp / (tp + fp) else NA,
        precision = if (tp + fp > 0) tp / (tp + fp) else NA,
        npv = if (tn + fn > 0) tn / (tn + fn) else NA,
        f1 = if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA,
        brier = mean((prob - as.numeric(y_true == levels(y_true)[2]))^2),
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
        mape = mean(abs((y_true - y_pred) / y_true)) * 100,
        NA
      )
      result[m] <- val
    }
  }
  
  result
}

#' @keywords internal
.bootstrap_metrics <- function(y_true, y_pred, task, metrics, n_boot) {
  n <- length(y_true)
  boot_mat <- matrix(NA, nrow = n_boot, ncol = length(metrics))
  colnames(boot_mat) <- metrics
  
  for (b in seq_len(n_boot)) {
    idx <- sample(n, replace = TRUE)
    y_b <- y_true[idx]
    pred_b <- if (is.matrix(y_pred)) y_pred[idx, ] else y_pred[idx]
    boot_mat[b, ] <- .calculate_all_metrics(y_b, pred_b, task, metrics)
  }
  
  # Calculate CIs
  ci_list <- list()
  for (m in metrics) {
    ci_list[[m]] <- c(
      lower = quantile(boot_mat[, m], 0.025, na.rm = TRUE),
      estimate = mean(boot_mat[, m], na.rm = TRUE),
      upper = quantile(boot_mat[, m], 0.975, na.rm = TRUE)
    )
  }
  
  ci_list
}

#' @keywords internal
.get_binary_truth_prob <- function(object, newdata = NULL, caller = "Evaluation") {
  if (object$task != "classification") {
    stop(sprintf("%s is only for classification models", caller))
  }

  if (is.null(newdata)) {
    y_true <- object$y_test
    pred <- ukb_ml_predict(object, type = "prob")
  } else {
    y_true <- newdata[[object$outcome]]
    pred <- ukb_ml_predict(object, newdata = newdata, type = "prob")
  }

  y_true <- droplevels(as.factor(y_true))
  prob <- if (is.matrix(pred)) pred[, 2] else pred

  keep <- !is.na(prob) & !is.na(y_true)
  if (!all(keep)) {
    y_true <- droplevels(y_true[keep])
    prob <- prob[keep]
  }

  if (length(prob) == 0) {
    stop(sprintf("%s: no valid observations after removing missing values", caller))
  }

  if (nlevels(y_true) != 2) {
    stop(sprintf("%s requires a binary outcome with exactly 2 classes", caller))
  }

  if (any(!is.finite(prob))) {
    stop(sprintf("%s: predicted probabilities contain non-finite values", caller))
  }

  if (any(prob < 0 | prob > 1)) {
    stop(sprintf("%s: predicted probabilities must be within [0, 1]", caller))
  }

  y_binary <- as.integer(y_true == levels(y_true)[2])
  list(y_true = y_true, prob = prob, y_binary = y_binary)
}

# ROC Curve

#' ROC Curve Analysis
#'
#' @description
#' Generate ROC curve and calculate AUC with optional confidence intervals.
#'
#' @param object A ukb_ml object or list of objects
#' @param newdata Optional new data
#' @param plot Whether to create ROC plot (default TRUE)
#' @param ci Compute confidence interval for AUC
#' @param ci_method Method: "delong" (default) or "bootstrap"
#' @param ... Additional arguments
#'
#' @return ukb_ml_roc object with ROC curve data
#'
#' @export
ukb_ml_roc <- function(object,
                       newdata = NULL,
                       plot = TRUE,
                       ci = TRUE,
                       ci_method = c("delong", "bootstrap"),
                       ...) {
  .ukb_ml_deprecated("ukb_ml_roc", "ukb_ml_evaluate_test() or ml$final_test_predictions")
  
  .check_ml_package("pROC")
  ci_method <- match.arg(ci_method)
  
  # Handle single model or list
  if (inherits(object, "ukb_ml")) {
    objects <- list(object)
    names(objects) <- .get_model_label(object$model_type)
  } else if (is.list(object)) {
    objects <- object
    if (is.null(names(objects))) {
      names(objects) <- sapply(objects, function(o) .get_model_label(o$model_type))
    }
  }
  
  roc_list <- list()
  auc_df <- data.frame(
    model = character(),
    auc = numeric(),
    lower = numeric(),
    upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(objects)) {
    obj <- objects[[name]]
    
    if (is.null(newdata)) {
      y_true <- obj$y_test
      pred <- ukb_ml_predict(obj, type = "prob")
    } else {
      y_true <- as.factor(newdata[[obj$outcome]])
      pred <- ukb_ml_predict(obj, newdata = newdata, type = "prob")
    }
    
    prob <- if (is.matrix(pred)) pred[, 2] else pred
    
    roc_obj <- roc(y_true, prob, quiet = TRUE)
    roc_list[[name]] <- roc_obj
    
    # AUC with CI
    auc_val <- as.numeric(auc(roc_obj))
    
    if (ci) {
      if (ci_method == "delong") {
        auc_ci <- as.numeric(ci.auc(roc_obj))
      } else {
        auc_ci <- as.numeric(ci.auc(roc_obj, method = "bootstrap"))
      }
      auc_df <- rbind(auc_df, data.frame(
        model = name,
        auc = auc_val,
        lower = auc_ci[1],
        upper = auc_ci[3],
        stringsAsFactors = FALSE
      ))
    } else {
      auc_df <- rbind(auc_df, data.frame(
        model = name,
        auc = auc_val,
        lower = NA,
        upper = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  result <- list(
    roc = roc_list,
    auc = auc_df
  )
  
  class(result) <- "ukb_ml_roc"
  
  if (plot) {
    result$plot <- plot_ml_roc(result)
  }
  
  result
}

#' @export
print.ukb_ml_roc <- function(x, ...) {
  cat("\nROC Analysis Results\n")
  cat("\n")
  print(x$auc, row.names = FALSE)
  invisible(x)
}

# Calibration Curve

#' Calibration Curve Analysis
#'
#' @description
#' Generate calibration curve to assess prediction reliability.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data
#' @param n_bins Number of bins for calibration (default 10)
#' @param method Smoothing method: "loess", "isotonic", or "none"
#' @param plot Whether to create calibration plot
#' @param ... Additional arguments
#'
#' @return ukb_ml_calibration object
#'
#' @export
ukb_ml_calibration <- function(object,
                               newdata = NULL,
                               n_bins = 10,
                               method = c("none", "loess", "isotonic"),
                               plot = TRUE,
                               ...) {
  .ukb_ml_deprecated("ukb_ml_calibration", "ukb_ml_evaluate_test() and ml$final_test_predictions")
  
  method <- match.arg(method)
  
  if (object$task != "classification") {
    stop("Calibration curves are only for classification models")
  }
  
  # Get predictions
  if (is.null(newdata)) {
    y_true <- object$y_test
    pred <- ukb_ml_predict(object, type = "prob")
  } else {
    y_true <- as.factor(newdata[[object$outcome]])
    pred <- ukb_ml_predict(object, newdata = newdata, type = "prob")
  }
  
  prob <- if (is.matrix(pred)) pred[, 2] else pred
  y_binary <- as.numeric(y_true == levels(y_true)[2])
  
  # Create bins
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bins <- cut(prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # Calculate observed vs predicted by bin
  cal_data <- data.frame(
    bin = seq_len(n_bins),
    predicted = numeric(n_bins),
    observed = numeric(n_bins),
    count = integer(n_bins),
    lower = numeric(n_bins),
    upper = numeric(n_bins)
  )
  
  for (i in seq_len(n_bins)) {
    idx <- which(bins == i)
    if (length(idx) > 0) {
      cal_data$predicted[i] <- mean(prob[idx])
      cal_data$observed[i] <- mean(y_binary[idx])
      cal_data$count[i] <- length(idx)
      
      # Wilson CI for proportion
      if (length(idx) >= 5) {
        prop <- mean(y_binary[idx])
        n <- length(idx)
        z <- 1.96
        denom <- 1 + z^2 / n
        center <- (prop + z^2 / (2 * n)) / denom
        margin <- z * sqrt((prop * (1 - prop) + z^2 / (4 * n)) / n) / denom
        cal_data$lower[i] <- max(0, center - margin)
        cal_data$upper[i] <- min(1, center + margin)
      } else {
        cal_data$lower[i] <- NA
        cal_data$upper[i] <- NA
      }
    } else {
      cal_data$predicted[i] <- (breaks[i] + breaks[i + 1]) / 2
      cal_data$observed[i] <- NA
      cal_data$count[i] <- 0
      cal_data$lower[i] <- NA
      cal_data$upper[i] <- NA
    }
  }
  
  # Remove empty bins
  cal_data <- cal_data[cal_data$count > 0, ]
  
  # Calculate calibration metrics
  # Brier score
  brier <- mean((prob - y_binary)^2)
  
  # Expected Calibration Error (ECE)
  ece <- sum(abs(cal_data$observed - cal_data$predicted) * cal_data$count) / sum(cal_data$count)
  
  result <- list(
    calibration = cal_data,
    brier_score = brier,
    ece = ece,
    method = method
  )
  
  class(result) <- "ukb_ml_calibration"
  
  if (plot) {
    result$plot <- plot_ml_calibration(result)
  }
  
  result
}

#' @export
print.ukb_ml_calibration <- function(x, ...) {
  cat("\nCalibration Analysis\n")
  cat("\n")
  cat(sprintf("Brier Score: %.4f\n", x$brier_score))
  cat(sprintf("Expected Calibration Error: %.4f\n", x$ece))
  cat("\nCalibration by bin:\n")
  print(x$calibration[, c("bin", "predicted", "observed", "count")], row.names = FALSE)
  invisible(x)
}

# Confusion Matrix

#' Confusion Matrix
#'
#' @description
#' Generate confusion matrix for classification model.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data
#' @param threshold Classification threshold (default 0.5)
#' @param plot Whether to create confusion matrix plot
#' @param ... Additional arguments
#'
#' @return ukb_ml_confusion object
#'
#' @export
ukb_ml_confusion <- function(object,
                             newdata = NULL,
                             threshold = 0.5,
                             plot = TRUE,
                             ...) {
  .ukb_ml_deprecated("ukb_ml_confusion", "ukb_ml_evaluate_test(threshold = ...)")
  
  if (object$task != "classification") {
    stop("Confusion matrix is only for classification models")
  }
  
  # Get predictions
  if (is.null(newdata)) {
    y_true <- object$y_test
    pred <- ukb_ml_predict(object, type = "prob")
  } else {
    y_true <- as.factor(newdata[[object$outcome]])
    pred <- ukb_ml_predict(object, newdata = newdata, type = "prob")
  }
  
  prob <- if (is.matrix(pred)) pred[, 2] else pred
  pred_class <- ifelse(prob > threshold, levels(y_true)[2], levels(y_true)[1])
  pred_class <- factor(pred_class, levels = levels(y_true))
  
  # Create confusion matrix
  cm <- table(Predicted = pred_class, Actual = y_true)
  
  # Calculate metrics
  tp <- cm[2, 2]
  tn <- cm[1, 1]
  fp <- cm[2, 1]
  fn <- cm[1, 2]
  
  metrics <- c(
    accuracy = (tp + tn) / sum(cm),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv = tp / (tp + fp),
    npv = tn / (tn + fn),
    f1 = 2 * tp / (2 * tp + fp + fn)
  )
  
  result <- list(
    confusion_matrix = cm,
    metrics = metrics,
    threshold = threshold
  )
  
  class(result) <- "ukb_ml_confusion"
  
  if (plot) {
    result$plot <- plot_ml_confusion(result)
  }
  
  result
}

#' @export
print.ukb_ml_confusion <- function(x, ...) {
  cat("\nConfusion Matrix (threshold =", x$threshold, ")\n")
  cat("\n")
  print(x$confusion_matrix)
  cat("\nMetrics:\n")
  for (m in names(x$metrics)) {
    cat(sprintf("  %s: %.4f\n", m, x$metrics[m]))
  }
  invisible(x)
}

# KS Curve

#' KS Curve Analysis
#'
#' @description
#' Compute Kolmogorov-Smirnov curve (TPR - FPR vs threshold) for a
#' binary classification model.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data for evaluation
#' @param plot Whether to create the KS plot (default TRUE)
#' @param n_thresholds Number of threshold points (default 200)
#' @param ... Additional arguments
#'
#' @return A ukb_ml_ks object with fields: data (threshold/tpr/fpr/ks),
#'   ks_stat (max KS), ks_threshold (threshold at max KS)
#'
#' @export
ukb_ml_ks <- function(object,
                      newdata = NULL,
                      plot = TRUE,
                      n_thresholds = 200,
                      ...) {
  .ukb_ml_deprecated("ukb_ml_ks", "ukb_ml_evaluate_test() and ml$final_test_predictions")

  n_thresholds <- as.integer(n_thresholds)
  if (length(n_thresholds) != 1 || is.na(n_thresholds) || n_thresholds < 2) {
    stop("n_thresholds must be a single integer >= 2")
  }

  eval_data <- .get_binary_truth_prob(object, newdata, caller = "ukb_ml_ks")
  prob <- eval_data$prob
  y_binary <- eval_data$y_binary
  n_pos <- sum(y_binary == 1L)
  n_neg <- sum(y_binary == 0L)

  thresholds <- seq(0, 1, length.out = n_thresholds)

  tpr_vec <- vapply(thresholds, function(t) {
    if (n_pos == 0) return(0)
    sum(prob >= t & y_binary == 1L) / n_pos
  }, numeric(1))

  fpr_vec <- vapply(thresholds, function(t) {
    if (n_neg == 0) return(0)
    sum(prob >= t & y_binary == 0L) / n_neg
  }, numeric(1))

  ks_vec  <- tpr_vec - fpr_vec
  max_idx <- which.max(ks_vec)

  result <- list(
    data          = data.frame(threshold = thresholds, tpr = tpr_vec,
                               fpr = fpr_vec, ks = ks_vec),
    ks_stat       = ks_vec[max_idx],
    ks_threshold  = thresholds[max_idx]
  )
  class(result) <- "ukb_ml_ks"

  if (plot) result$plot <- plot_ml_ks(result)
  result
}

#' @export
print.ukb_ml_ks <- function(x, ...) {
  cat("\nKS Curve Analysis\n")
  cat(sprintf("  KS Statistic : %.4f\n", x$ks_stat))
  cat(sprintf("  At Threshold : %.4f\n", x$ks_threshold))
  invisible(x)
}

# PR Curve

#' Precision-Recall Curve Analysis
#'
#' @description
#' Compute Precision-Recall curve and area under PR curve (AUPRC)
#' for a binary classification model.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data
#' @param plot Whether to create the PR plot (default TRUE)
#' @param n_thresholds Number of threshold points (default 200)
#' @param ... Additional arguments
#'
#' @return A ukb_ml_pr object with fields: data (threshold/precision/recall),
#'   auprc, prevalence
#'
#' @export
ukb_ml_pr <- function(object,
                      newdata = NULL,
                      plot = TRUE,
                      n_thresholds = 200,
                      ...) {
  .ukb_ml_deprecated("ukb_ml_pr", "ukb_ml_evaluate_test() and ml$final_test_predictions")

  n_thresholds <- as.integer(n_thresholds)
  if (length(n_thresholds) != 1 || is.na(n_thresholds) || n_thresholds < 2) {
    stop("n_thresholds must be a single integer >= 2")
  }

  eval_data <- .get_binary_truth_prob(object, newdata, caller = "ukb_ml_pr")
  prob <- eval_data$prob
  y_binary <- eval_data$y_binary
  prevalence <- mean(y_binary)

  thresholds <- seq(0, 1, length.out = n_thresholds)

  precision_vec <- vapply(thresholds, function(t) {
    tp <- sum(prob >= t & y_binary == 1L)
    fp <- sum(prob >= t & y_binary == 0L)
    if (tp + fp == 0) return(NA_real_)
    tp / (tp + fp)
  }, numeric(1))

  recall_vec <- vapply(thresholds, function(t) {
    tp <- sum(prob >= t & y_binary == 1L)
    fn <- sum(prob < t  & y_binary == 1L)
    if (tp + fn == 0) return(0)
    tp / (tp + fn)
  }, numeric(1))

  pr_data <- data.frame(threshold = thresholds,
                        precision = precision_vec,
                        recall    = recall_vec)

  # AUPRC via trapezoidal integration on valid rows, sorted by recall
  valid     <- !is.na(pr_data$precision)
  pr_sorted <- pr_data[valid, ][order(pr_data$recall[valid]), ]
  auprc <- if (nrow(pr_sorted) >= 2) {
    sum(diff(pr_sorted$recall) *
          (head(pr_sorted$precision, -1) + utils::tail(pr_sorted$precision, -1)) / 2)
  } else {
    NA_real_
  }

  result <- list(
    data       = pr_data,
    auprc      = auprc,
    prevalence = prevalence
  )
  class(result) <- "ukb_ml_pr"

  if (plot) result$plot <- plot_ml_pr(result)
  result
}

#' @export
print.ukb_ml_pr <- function(x, ...) {
  cat("\nPR Curve Analysis\n")
  cat(sprintf("  AUPRC      : %.4f\n", x$auprc))
  cat(sprintf("  Prevalence : %.4f\n", x$prevalence))
  invisible(x)
}

# Gain / Lift Curves

#' Gain and Lift Curve Analysis
#'
#' @description
#' Compute Gain and Lift curves for a binary classification model by
#' ranking predictions into decile bins.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data
#' @param plot Whether to create gain and lift plots (default TRUE)
#' @param n_bins Number of bins / deciles (default 10)
#' @param ... Additional arguments
#'
#' @return A ukb_ml_gain_lift object with field data containing:
#'   decile, population_pct, positive_capture_pct, gain, lift
#'
#' @export
ukb_ml_gain_lift <- function(object,
                             newdata = NULL,
                             plot = TRUE,
                             n_bins = 10,
                             ...) {
  .ukb_ml_deprecated("ukb_ml_gain_lift", "ukb_ml_evaluate_test() and ml$final_test_predictions")

  n_bins <- as.integer(n_bins)
  if (length(n_bins) != 1 || is.na(n_bins) || n_bins < 1) {
    stop("n_bins must be a single integer >= 1")
  }

  eval_data <- .get_binary_truth_prob(object, newdata, caller = "ukb_ml_gain_lift")
  prob <- eval_data$prob
  y_binary <- eval_data$y_binary
  n <- length(y_binary)
  total_pos <- sum(y_binary)

  if (n_bins > n) {
    warning(sprintf("n_bins (%d) is larger than sample size (%d); using n_bins = %d", n_bins, n, n),
      call. = FALSE
    )
    n_bins <- n
  }

  # Sort descending by predicted probability
  ord      <- order(prob, decreasing = TRUE)
  y_sorted <- y_binary[ord]

  # Assign to bins (ceiling guarantees each observation falls in [1, n_bins])
  bins <- pmin(ceiling(seq_along(y_sorted) * n_bins / n), n_bins)

  gl_data <- data.frame(
    decile               = seq_len(n_bins),
    population_pct       = numeric(n_bins),
    positive_capture_pct = numeric(n_bins),
    gain                 = numeric(n_bins),
    lift                 = numeric(n_bins)
  )

  for (k in seq_len(n_bins)) {
    in_top_k               <- bins <= k
    pop_pct                <- sum(in_top_k) / n
    pos_pct                <- if (total_pos > 0) sum(y_sorted[in_top_k]) / total_pos else 0
    gl_data$population_pct[k]       <- pop_pct
    gl_data$positive_capture_pct[k] <- pos_pct
    gl_data$gain[k]                  <- pos_pct
    gl_data$lift[k]                  <- if (pop_pct > 0) pos_pct / pop_pct else 0
  }

  result <- list(data = gl_data, n_bins = n_bins)
  class(result) <- "ukb_ml_gain_lift"

  if (plot) {
    result$gain_plot <- plot_ml_gain(result)
    result$lift_plot <- plot_ml_lift(result)
    result$plot      <- result$gain_plot
  }
  result
}

#' @export
print.ukb_ml_gain_lift <- function(x, ...) {
  cat("\nGain / Lift Table\n")
  print(x$data, row.names = FALSE)
  invisible(x)
}

# DCA

#' Decision Curve Analysis
#'
#' @description
#' Compute Decision Curve Analysis (DCA) net benefit across a range of
#' threshold probabilities for a binary classification model.
#'
#' @param object A ukb_ml object
#' @param newdata Optional new data
#' @param plot Whether to create the DCA plot (default TRUE)
#' @param thresholds Numeric vector of threshold probabilities
#'   (default seq(0.01, 0.99, by = 0.01))
#' @param harm Additional harm parameter subtracted from net benefit (default 0)
#' @param ... Additional arguments
#'
#' @return A ukb_ml_dca object with field data containing:
#'   threshold, net_benefit_model, net_benefit_all, net_benefit_none
#'
#' @export
ukb_ml_dca <- function(object,
                       newdata    = NULL,
                       plot       = TRUE,
                       thresholds = seq(0.01, 0.99, by = 0.01),
                       harm       = 0,
                       ...) {
  .ukb_ml_deprecated("ukb_ml_dca", "ukb_ml_evaluate_test() and ml$final_test_predictions")

  if (!is.numeric(thresholds) || length(thresholds) == 0) {
    stop("thresholds must be a non-empty numeric vector")
  }

  thresholds <- sort(unique(as.numeric(thresholds)))
  if (any(!is.finite(thresholds))) {
    stop("thresholds must contain only finite numeric values")
  }

  if (any(thresholds <= 0 | thresholds >= 1)) {
    stop("thresholds must be strictly between 0 and 1")
  }

  if (!is.numeric(harm) || length(harm) != 1 || is.na(harm) || !is.finite(harm)) {
    stop("harm must be a single finite numeric value")
  }

  eval_data <- .get_binary_truth_prob(object, newdata, caller = "ukb_ml_dca")
  prob <- eval_data$prob
  y_binary <- eval_data$y_binary
  n <- length(y_binary)
  prevalence <- mean(y_binary)

  nb_model <- vapply(thresholds, function(t) {
    tp <- sum(prob >= t & y_binary == 1L)
    fp <- sum(prob >= t & y_binary == 0L)
    tp / n - fp / n * (t / (1 - t)) - harm
  }, numeric(1))

  nb_all <- vapply(thresholds, function(t) {
    prevalence - (1 - prevalence) * (t / (1 - t))
  }, numeric(1))

  result <- list(
    data = data.frame(
      threshold          = thresholds,
      net_benefit_model  = nb_model,
      net_benefit_all    = nb_all,
      net_benefit_none   = 0
    ),
    prevalence = prevalence,
    harm       = harm
  )
  class(result) <- "ukb_ml_dca"

  if (plot) result$plot <- plot_ml_dca(result)
  result
}

#' @export
print.ukb_ml_dca <- function(x, ...) {
  cat("\nDecision Curve Analysis\n")
  cat(sprintf("  Prevalence : %.4f\n", x$prevalence))
  cat(sprintf("  Harm       : %.4f\n", x$harm))
  cat("\nNet Benefit (first 5 thresholds):\n")
  print(head(x$data, 5), row.names = FALSE)
  invisible(x)
}
