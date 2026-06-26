#' Multiple Imputation Result Pooling
#'
#' @description
#' Functions for combining results from analyses performed on multiply imputed datasets
#' using Rubin's Rules. Supports various regression models including linear, logistic,
#' Poisson, Cox, and negative binomial regression.
#'
#' @name mi_pool
#' @keywords internal
#' @importFrom mitools MIcombine MIextract imputationList
#' @importFrom survival coxph
#' @importFrom MASS glm.nb
NULL


#' Pool Results from Multiple Imputation Models
#'
#' @description
#' Combines results of regression analyses performed on multiply imputed datasets
#' using Rubin's Rules via the mitools package.
#'
#' @param models A list of fitted model objects (one per imputed dataset).
#'   If NULL, models will be fitted using \code{datasets} and \code{formula}.
#' @param datasets A list of data.frames or an \code{imputationList} object.
#'   Required if \code{models} is NULL.
#' @param formula A formula specifying the model. Required if \code{datasets} is provided.
#' @param model_type Character string specifying the model type:
#'   \describe{
#'     \item{"lm"}{Linear regression}
#'     \item{"logistic"}{Logistic regression (GLM with binomial family)}
#'     \item{"poisson"}{Poisson regression}
#'     \item{"cox"}{Cox proportional hazards model}
#'     \item{"negbin"}{Negative binomial regression}
#'   }
#' @param family A \code{family} object for GLM. If NULL, inferred from \code{model_type}.
#' @param df.complete Complete-data degrees of freedom for small-sample correction.
#'   Default is Inf (large sample approximation).
#' @param conf.level Confidence level for intervals. Default 0.95.
#' @param exponentiate Logical; whether to exponentiate coefficients (for OR/HR/RR).
#'   If NULL, automatically determined based on model type.
#'
#' @return An object of class \code{mi_pooled_result} containing:
#'   \describe{
#'     \item{pooled}{Data frame with pooled estimates, standard errors, CIs, p-values, and FMI}
#'     \item{mi_result}{The raw \code{MIresult} object from mitools}
#'     \item{n_imputations}{Number of imputed datasets}
#'     \item{model_type}{The model type used}
#'     \item{formula}{The model formula}
#'     \item{exponentiated}{Whether estimates are exponentiated}
#'     \item{call}{The function call}
#'   }
#'
#' @export
pool_mi_models <- function(models = NULL,
                           datasets = NULL,
                           formula = NULL,
                           model_type = c("lm", "logistic", "poisson", "cox", "negbin"),
                           family = NULL,
                           df.complete = Inf,
                           conf.level = 0.95,
                           exponentiate = NULL) {

  model_type <- match.arg(model_type)
  call <- match.call()

  # Validate inputs
if (is.null(models) && is.null(datasets)) {
    stop("Either 'models' or 'datasets' must be provided.")
  }

  if (!is.null(datasets) && is.null(formula)) {
    stop("'formula' is required when 'datasets' is provided.")
  }

  # Fit models if not provided
  if (is.null(models)) {
    models <- fit_mi_models(
      datasets = datasets,
      formula = formula,
      model_type = model_type,
      family = family
    )
  }

  # Validate models list
  if (!is.list(models) || length(models) < 2) {
    stop("'models' must be a list with at least 2 model objects.")
  }

  n_imputations <- length(models)

  # Extract formula if not provided
  if (is.null(formula)) {
    formula <- tryCatch(
      stats::formula(models[[1]]),
      error = function(e) NULL
    )
  }

  # Combine using mitools
  mi_result <- tryCatch({
    MIcombine(models, df.complete = df.complete)
  }, error = function(e) {
    # Fallback: extract coefficients and variances manually
    betas <- MIextract(models, fun = stats::coef)
    vars <- MIextract(models, fun = stats::vcov)
    MIcombine(betas, vars, df.complete = df.complete)
  })

  # Build pooled results table
  pooled <- .build_pooled_table(mi_result, conf.level)

  # Determine exponentiation
  if (is.null(exponentiate)) {
    exponentiate <- model_type %in% c("logistic", "poisson", "cox", "negbin")
  }

  # Apply exponentiation if requested
  exponentiated <- FALSE
  if (exponentiate && model_type != "lm") {
    pooled$estimate <- exp(pooled$estimate)
    pooled$conf.low <- exp(pooled$conf.low)
    pooled$conf.high <- exp(pooled$conf.high)
    exponentiated <- TRUE
  }

  # Create result object
  result <- list(
    pooled = pooled,
    mi_result = mi_result,
    n_imputations = n_imputations,
    model_type = model_type,
    formula = formula,
    exponentiated = exponentiated,
    conf.level = conf.level,
    call = call
  )

  class(result) <- "mi_pooled_result"
  return(result)
}


#' Fit Regression Models on Multiple Imputed Datasets
#'
#' @description
#' Fits the specified regression model on each imputed dataset.
#'
#' @param datasets A list of data.frames or an \code{imputationList} object.
#' @param formula A formula specifying the model.
#' @param model_type Character string specifying the model type.
#' @param family A \code{family} object for GLM (optional).
#' @param ... Additional arguments passed to the model fitting function.
#'
#' @return A list of fitted model objects.
#'
#' @export
fit_mi_models <- function(datasets,
                          formula,
                          model_type = c("lm", "logistic", "poisson", "cox", "negbin"),
                          family = NULL,
                          ...) {

  model_type <- match.arg(model_type)

  # Convert imputationList to regular list if needed
  if (inherits(datasets, "imputationList")) {
    datasets <- datasets$imputations
  }

  if (!is.list(datasets)) {
    stop("'datasets' must be a list of data.frames or an imputationList.")
  }

  # Determine fitting function
  fit_fn <- switch(model_type,
    "lm" = function(d) stats::lm(formula, data = d, ...),
    "logistic" = function(d) {
      fam <- if (is.null(family)) stats::binomial() else family
      stats::glm(formula, data = d, family = fam, ...)
    },
    "poisson" = function(d) {
      fam <- if (is.null(family)) stats::poisson() else family
      stats::glm(formula, data = d, family = fam, ...)
    },
    "cox" = function(d) coxph(formula, data = d, ...),
    "negbin" = function(d) glm.nb(formula, data = d, ...)
  )

  # Fit models
  models <- lapply(seq_along(datasets), function(i) {
    tryCatch({
      fit_fn(datasets[[i]])
    }, error = function(e) {
      warning(sprintf("Model fitting failed for imputation %d: %s", i, e$message))
      NULL
    })
  })

  # Remove failed fits
  failed <- sapply(models, is.null)
  if (any(failed)) {
    warning(sprintf("%d out of %d model fits failed.", sum(failed), length(models)))
    models <- models[!failed]
  }

  if (length(models) < 2) {
    stop("At least 2 successful model fits are required.")
  }

  return(models)
}


#' Create an imputationList Object
#'
#' @description
#' Converts a list of data.frames to a \code{imputationList} object for use
#' with mitools functions.
#'
#' @param datasets A list of data.frames (imputed datasets).
#' @param validate Logical; whether to validate that all datasets have the
#'   same structure. Default TRUE.
#'
#' @return An \code{imputationList} object.
#'
#' @export
create_imputation_list <- function(datasets, validate = TRUE) {
  if (!is.list(datasets) || length(datasets) < 2) {
    stop("'datasets' must be a list of at least 2 data.frames.")
  }

  # Validate structure consistency
  if (validate) {
    ref_names <- names(datasets[[1]])
    ref_nrow <- nrow(datasets[[1]])

    for (i in 2:length(datasets)) {
      if (!identical(names(datasets[[i]]), ref_names)) {
        warning(sprintf(
          "Column names differ between imputation 1 and %d. Results may be unreliable.", i
        ))
      }
      if (nrow(datasets[[i]]) != ref_nrow) {
        warning(sprintf(
          "Row counts differ: imputation 1 has %d rows, imputation %d has %d rows.",
          ref_nrow, i, nrow(datasets[[i]])
        ))
      }
    }
  }

  imputationList(datasets)
}


#' Pool Custom Estimates from Multiple Imputations
#'
#' @description
#' Combines custom parameter estimates (not limited to regression coefficients)
#' from multiply imputed datasets using Rubin's Rules.
#'
#' @param estimates A list of numeric vectors containing point estimates from
#'   each imputed dataset. All vectors must have the same length.
#' @param variances A list of variance-covariance matrices (or single variances
#'   as 1x1 matrices) corresponding to the estimates.
#' @param df.complete Complete-data degrees of freedom. Default Inf.
#' @param conf.level Confidence level for intervals. Default 0.95.
#' @param labels Character vector of labels for the estimates. If NULL,
#'   names are taken from the first estimate vector or generated as "est1", "est2", etc.
#'
#' @return An object of class \code{mi_pooled_result}.
#'
#' @export
pool_custom_estimates <- function(estimates,
                                   variances,
                                   df.complete = Inf,
                                   conf.level = 0.95,
                                   labels = NULL) {
  if (!is.list(estimates) || !is.list(variances)) {
    stop("'estimates' and 'variances' must be lists.")
  }

  if (length(estimates) != length(variances)) {
    stop("'estimates' and 'variances' must have the same length.")
  }

  if (length(estimates) < 2) {
    stop("At least 2 sets of estimates are required.")
  }

  # Ensure estimates are numeric vectors
  estimates <- lapply(estimates, as.numeric)

  # Ensure variances are matrices
  variances <- lapply(variances, function(v) {
    if (is.matrix(v)) v else matrix(v, nrow = length(estimates[[1]]))
  })

  # Combine using mitools
  mi_result <- MIcombine(estimates, variances, df.complete = df.complete)

  # Build pooled table
  pooled <- .build_pooled_table(mi_result, conf.level)

  # Add labels
  if (!is.null(labels)) {
    if (length(labels) == nrow(pooled)) {
      pooled$term <- labels
    } else if (length(labels) == 1) {
      pooled$term <- paste0(labels, seq_len(nrow(pooled)))
    }
  } else if (!is.null(names(estimates[[1]]))) {
    pooled$term <- names(estimates[[1]])
  }

  result <- list(
    pooled = pooled,
    mi_result = mi_result,
    n_imputations = length(estimates),
    model_type = "custom",
    formula = NULL,
    exponentiated = FALSE,
    conf.level = conf.level,
    call = match.call()
  )

  class(result) <- "mi_pooled_result"
  return(result)
}


# -----------------------------------------------------------------------------
# Internal helper functions
# -----------------------------------------------------------------------------

#' Build Pooled Results Table from MIresult
#' @noRd
.build_pooled_table <- function(mi_result, conf.level = 0.95) {

  # Extract components
  estimates <- mi_result$coefficients
  se <- sqrt(diag(mi_result$variance))
  df <- mi_result$df

  # Get FMI (fraction of missing information)
  # FMI = (B + B/m) / T where B is between-imputation variance, T is total variance
  fmi <- mi_result$missinfo
  if (is.null(fmi)) {
    fmi <- rep(NA_real_, length(estimates))
  }

  # Calculate test statistics and p-values
  z <- estimates / se
  alpha <- 1 - conf.level

  # Use t-distribution if df is finite
  if (any(is.finite(df))) {
    p.value <- 2 * stats::pt(-abs(z), df = df)
    t_crit <- stats::qt(1 - alpha / 2, df = df)
  } else {
    p.value <- 2 * stats::pnorm(-abs(z))
    t_crit <- stats::qnorm(1 - alpha / 2)
  }

  # Confidence intervals
  conf.low <- estimates - t_crit * se
  conf.high <- estimates + t_crit * se

  # Term names
  term <- names(estimates)
  if (is.null(term)) {
    term <- paste0("V", seq_along(estimates))
  }

  data.frame(
    term = term,
    estimate = estimates,
    std.error = se,
    statistic = z,
    df = df,
    p.value = p.value,
    conf.low = conf.low,
    conf.high = conf.high,
    fmi = as.numeric(fmi),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


# -----------------------------------------------------------------------------
# S3 Methods for mi_pooled_result
# -----------------------------------------------------------------------------

#' @export
print.mi_pooled_result <- function(x, digits = 3, ...) {
  cat("Multiple Imputation Pooled Results\n")
  cat("-----------------------------------\n")
  cat(sprintf("Model type: %s\n", x$model_type))
  cat(sprintf("Imputations: %d\n", x$n_imputations))

  if (!is.null(x$formula)) {
    cat(sprintf("Formula: %s\n", deparse(x$formula)))
  }

  if (x$exponentiated) {
    est_label <- switch(x$model_type,
      "logistic" = "OR",
      "cox" = "HR",
      "poisson" = "RR",
      "negbin" = "RR",
      "exp(est)"
    )
    cat(sprintf("Estimates are exponentiated (%s)\n", est_label))
  }

  cat("\n")

  # Format output
  out <- x$pooled
  out$estimate <- round(out$estimate, digits)
  out$std.error <- round(out$std.error, digits)
  out$statistic <- round(out$statistic, digits)
  out$p.value <- format.pval(out$p.value, digits = digits)
  out$conf.low <- round(out$conf.low, digits)
  out$conf.high <- round(out$conf.high, digits)
  out$fmi <- round(out$fmi, digits)

  print(out, row.names = FALSE)
  invisible(x)
}


#' @export
summary.mi_pooled_result <- function(object, exponentiate = NULL, ...) {
  cat("Multiple Imputation Results (Rubin's Rules)\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat(sprintf("Number of imputations: %d\n", object$n_imputations))
  cat(sprintf("Model type: %s\n", object$model_type))
  cat(sprintf("Confidence level: %.0f%%\n\n", object$conf.level * 100))

  # Determine display format
  if (is.null(exponentiate)) {
    exponentiate <- object$exponentiated
  }

  out <- object$pooled

  # Re-exponentiate or un-exponentiate if needed
  if (exponentiate && !object$exponentiated && object$model_type != "lm") {
    out$estimate <- exp(out$estimate)
    out$conf.low <- exp(out$conf.low)
    out$conf.high <- exp(out$conf.high)
  } else if (!exponentiate && object$exponentiated) {
    out$estimate <- log(out$estimate)
    out$conf.low <- log(out$conf.low)
    out$conf.high <- log(out$conf.high)
  }

  # Add effect type label
  est_label <- if (exponentiate && object$model_type != "lm") {
    switch(object$model_type,
      "logistic" = "OR",
      "cox" = "HR",
      "poisson" = "RR",
      "negbin" = "RR",
      "exp(est)"
    )
  } else {
    "Estimate"
  }

  # Format for display
  cat("Pooled Coefficients:\n")
  cat(sprintf("\n  %-20s %10s %10s %12s %10s %8s\n",
              "Term", est_label, "Std.Err", "95% CI", "p-value", "FMI"))
  cat(paste(rep("-", 78), collapse = ""), "\n")

  for (i in seq_len(nrow(out))) {
    ci_str <- sprintf("[%.3f, %.3f]", out$conf.low[i], out$conf.high[i])
    p_str <- format.pval(out$p.value[i], digits = 3)
    fmi_str <- if (is.na(out$fmi[i])) "NA" else sprintf("%.1f%%", out$fmi[i] * 100)

    cat(sprintf("  %-20s %10.3f %10.3f %12s %10s %8s\n",
                substr(out$term[i], 1, 20),
                out$estimate[i],
                out$std.error[i],
                ci_str,
                p_str,
                fmi_str))
  }

  cat("\n")
  cat("FMI = Fraction of Missing Information\n")

  invisible(object)
}


#' @export
coef.mi_pooled_result <- function(object, ...) {
  est <- object$pooled$estimate
  names(est) <- object$pooled$term
  return(est)
}


#' @export
confint.mi_pooled_result <- function(object, parm = NULL, level = 0.95, ...) {
  # If different level requested, recalculate
  if (level != object$conf.level) {
    mi_res <- object$mi_result
    estimates <- mi_res$coefficients
    se <- sqrt(diag(mi_res$variance))
    df <- mi_res$df
    alpha <- 1 - level

    if (any(is.finite(df))) {
      t_crit <- stats::qt(1 - alpha / 2, df = df)
    } else {
      t_crit <- stats::qnorm(1 - alpha / 2)
    }

    ci <- cbind(
      estimates - t_crit * se,
      estimates + t_crit * se
    )
  } else {
    ci <- cbind(object$pooled$conf.low, object$pooled$conf.high)
  }

  rownames(ci) <- object$pooled$term
  colnames(ci) <- paste0(c(100 * (1 - level) / 2, 100 * (1 - (1 - level) / 2)), "%")

  # Subset if parm specified
  if (!is.null(parm)) {
    if (is.numeric(parm)) {
      ci <- ci[parm, , drop = FALSE]
    } else {
      ci <- ci[parm, , drop = FALSE]
    }
  }

  return(ci)
}


#' @export
vcov.mi_pooled_result <- function(object, ...) {
  v <- object$mi_result$variance
  rownames(v) <- colnames(v) <- object$pooled$term
  return(v)
}


#' Tidy Method for mi_pooled_result
#'
#' @description
#' Returns a tidy data frame of pooled estimates, compatible with broom package style.
#'
#' @param x An \code{mi_pooled_result} object.
#' @param conf.int Logical; include confidence intervals? Default TRUE.
#' @param conf.level Confidence level. Default 0.95.
#' @param exponentiate Logical; exponentiate estimates? Default FALSE.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns: term, estimate, std.error, statistic, p.value,
#'   and optionally conf.low, conf.high, fmi.
#'
#' @export
tidy.mi_pooled_result <- function(x,
                                   conf.int = TRUE,
                                   conf.level = 0.95,
                                   exponentiate = FALSE,
                                   ...) {
  out <- x$pooled

  # Handle exponentiation
  if (exponentiate && !x$exponentiated && x$model_type != "lm") {
    out$estimate <- exp(out$estimate)
    out$conf.low <- exp(out$conf.low)
    out$conf.high <- exp(out$conf.high)
  } else if (!exponentiate && x$exponentiated) {
    out$estimate <- log(out$estimate)
    out$conf.low <- log(out$conf.low)
    out$conf.high <- log(out$conf.high)
  }

  # Recalculate CI if different level
  if (conf.level != x$conf.level) {
    mi_res <- x$mi_result
    estimates_raw <- if (x$exponentiated && !exponentiate) {
      log(x$pooled$estimate)
    } else if (!x$exponentiated && exponentiate) {
      log(out$estimate)
    } else {
      mi_res$coefficients
    }

    se <- sqrt(diag(mi_res$variance))
    df <- mi_res$df
    alpha <- 1 - conf.level

    if (any(is.finite(df))) {
      t_crit <- stats::qt(1 - alpha / 2, df = df)
    } else {
      t_crit <- stats::qnorm(1 - alpha / 2)
    }

    out$conf.low <- estimates_raw - t_crit * se
    out$conf.high <- estimates_raw + t_crit * se

    if (exponentiate && x$model_type != "lm") {
      out$conf.low <- exp(out$conf.low)
      out$conf.high <- exp(out$conf.high)
    }
  }

  # Select columns
  if (conf.int) {
    out <- out[, c("term", "estimate", "std.error", "statistic", "p.value",
                   "conf.low", "conf.high", "fmi")]
  } else {
    out <- out[, c("term", "estimate", "std.error", "statistic", "p.value", "fmi")]
  }

  return(out)
}
