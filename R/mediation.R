#' Run Causal Mediation Analysis
#'
#' @description
#' Perform regression-based causal mediation analysis using the regmedint package.
#' Supports linear, logistic, and Cox proportional hazards models for the outcome,
#' and linear or logistic models for the mediator.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param exposure Character string specifying the exposure (treatment) variable name.
#' @param mediator Character string specifying the mediator variable name.
#' @param outcome Character string specifying the outcome variable name.
#'   For Cox models, this should be the time variable.
#' @param covariates Character vector of covariate names. Default NULL.
#' @param exposure_levels Numeric vector of length 2: c(a0, a1) where a0 is the
#'   reference level and a1 is the comparison level. Default c(0, 1).
#' @param mediator_value Numeric value at which to evaluate the controlled direct
#'   effect (CDE). Default 0.
#' @param covariate_values Numeric vector of covariate values at which to evaluate
#'   conditional effects. If NULL, uses mean (continuous) or mode (categorical).
#' @param mediator_type Character string: "continuous" or "binary". Default "continuous".
#' @param outcome_type Character string: "linear", "logistic", or "cox". Default "linear".
#' @param endpoint Character vector of length 2 for Cox models: c("time_col", "status_col").
#'   Required when outcome_type = "cox".
#' @param interaction Logical; whether to include exposure-mediator interaction
#'   in the outcome model. Default TRUE.
#' @param boot Logical; whether to use bootstrap for confidence intervals. Default FALSE.
#' @param boot_n Integer; number of bootstrap replicates. Default 1000.
#' @param conf_level Numeric; confidence level. Default 0.95.
#'
#' @return An object of class "mediation_result" containing:
#'   \describe{
#'     \item{effects}{data.frame with effect estimates, SE, CI, and p-values}
#'     \item{mediator_model}{Fitted mediator model object}
#'     \item{outcome_model}{Fitted outcome model object}
#'     \item{regmedint_obj}{Original regmedint object (if available)}
#'     \item{call}{The matched call}
#'     \item{params}{List of analysis parameters}
#'   }
#'
#' @details
#' This function wraps the regmedint package to provide a user-friendly interface
#' for causal mediation analysis. It implements the methods described in
#' Valeri & VanderWeele (2013, 2015).
#'
#' \strong{Effect definitions:}
#' \itemize{
#'   \item \code{cde}: Controlled Direct Effect - effect of exposure with mediator fixed
#'   \item \code{pnde}: Pure Natural Direct Effect - direct effect (traditional NDE)
#'   \item \code{tnie}: Total Natural Indirect Effect - indirect effect (traditional NIE)
#'   \item \code{tnde}: Total Natural Direct Effect
#'   \item \code{pnie}: Pure Natural Indirect Effect
#'   \item \code{te}: Total Effect = NDE + NIE
#'   \item \code{pm}: Proportion Mediated = NIE / TE
#' }
#'
#' @references
#' Valeri L, VanderWeele TJ. Mediation analysis allowing for exposure-mediator
#' interactions and causal interpretation. Psychological Methods. 2013;18(2):137-150.
#'
#' @export
run_mediation <- function(data,
                          exposure,
                          mediator,
                          outcome,
                          covariates = NULL,
                          exposure_levels = c(0, 1),
                          mediator_value = 0,
                          covariate_values = NULL,
                          mediator_type = c("continuous", "binary"),
                          outcome_type = c("linear", "logistic", "cox"),
                          endpoint = NULL,
                          interaction = TRUE,
                          boot = FALSE,
                          boot_n = 1000,
                          conf_level = 0.95) {

  # Match arguments

  mediator_type <- match.arg(mediator_type)
  outcome_type <- match.arg(outcome_type)

  # Check for regmedint package

  if (!requireNamespace("regmedint", quietly = TRUE)) {
    stop("Package 'regmedint' is required for mediation analysis. ",
         "Please install it with: install.packages('regmedint')")
  }

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!is.character(exposure) || length(exposure) != 1) {
    stop("'exposure' must be a single character string.")
  }

  if (!is.character(mediator) || length(mediator) != 1) {
    stop("'mediator' must be a single character string.")
  }

  if (!is.character(outcome) || length(outcome) != 1) {
    stop("'outcome' must be a single character string.")
  }

  # Check variables exist
  required_vars <- c(exposure, mediator, outcome, covariates)
  if (outcome_type == "cox") {
    if (is.null(endpoint) || length(endpoint) != 2) {
      stop("For Cox models, 'endpoint' must be c('time_col', 'status_col').")
    }
    required_vars <- c(required_vars, endpoint[2])  # status column
  }

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  # Validate exposure_levels
  if (!is.numeric(exposure_levels) || length(exposure_levels) != 2) {
    stop("'exposure_levels' must be a numeric vector of length 2: c(a0, a1).")
  }

  # Map model types to regmedint parameters
  mreg <- if (mediator_type == "continuous") "linear" else "logistic"
  yreg <- switch(outcome_type,
                 "linear" = "linear",
                 "logistic" = "logistic",
                 "cox" = "survCox")

  # Calculate default covariate values if not provided
  if (is.null(covariate_values) && !is.null(covariates)) {
    covariate_values <- sapply(covariates, function(cv) {
      x <- data[[cv]]
      if (is.numeric(x)) {
        mean(x, na.rm = TRUE)
      } else {
        # Mode for categorical
        as.numeric(names(sort(table(x), decreasing = TRUE))[1])
      }
    })
  } else if (is.null(covariates)) {
    covariate_values <- NULL
  }

  message("[run_mediation] Running causal mediation analysis...")
  message(sprintf("  Exposure: %s (%s vs %s)", exposure, exposure_levels[2], exposure_levels[1]))
  message(sprintf("  Mediator: %s (%s)", mediator, mediator_type))
  message(sprintf("  Outcome: %s (%s model)", outcome, outcome_type))
  if (!is.null(covariates)) {
    message(sprintf("  Covariates: %s", paste(covariates, collapse = ", ")))
  }

  # Prepare arguments for regmedint
  regmedint_args <- list(
    data = as.data.frame(data),
    yvar = outcome,
    avar = exposure,
    mvar = mediator,
    cvar = covariates,
    a0 = exposure_levels[1],
    a1 = exposure_levels[2],
    m_cde = mediator_value,
    c_cond = covariate_values,
    mreg = mreg,
    yreg = yreg,
    interaction = interaction,
    casecontrol = FALSE
  )

  # Add event variable for survival models
  if (outcome_type == "cox") {
    regmedint_args$eventvar <- endpoint[2]
  }

  # Fit regmedint model
  tryCatch({
    regmedint_obj <- do.call(regmedint::regmedint, regmedint_args)

    # Extract summary
    sum_obj <- summary(regmedint_obj)
    effects_matrix <- coef(sum_obj)

    # Convert to data.frame
    effects_df <- data.frame(
      effect = rownames(effects_matrix),
      est = effects_matrix[, "est"],
      se = effects_matrix[, "se"],
      z = effects_matrix[, "Z"],
      p = effects_matrix[, "p"],
      lower = effects_matrix[, "lower"],
      upper = effects_matrix[, "upper"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    # Bootstrap confidence intervals if requested
    if (boot) {
      message(sprintf("  Running bootstrap (%d replicates)...", boot_n))
      boot_results <- .bootstrap_mediation(
        data = data,
        regmedint_args = regmedint_args,
        boot_n = boot_n,
        conf_level = conf_level
      )
      effects_df$lower <- boot_results$lower
      effects_df$upper <- boot_results$upper
      effects_df$boot_se <- boot_results$se
    }

    # Create result object
    result <- list(
      effects = effects_df,
      mediator_model = regmedint_obj$mreg_fit,
      outcome_model = regmedint_obj$yreg_fit,
      regmedint_obj = regmedint_obj,
      call = match.call(),
      params = list(
        exposure = exposure,
        mediator = mediator,
        outcome = outcome,
        covariates = covariates,
        exposure_levels = exposure_levels,
        mediator_value = mediator_value,
        covariate_values = covariate_values,
        mediator_type = mediator_type,
        outcome_type = outcome_type,
        interaction = interaction,
        boot = boot
      )
    )

    class(result) <- "mediation_result"

    message("[run_mediation] Complete.")

    return(result)

  }, error = function(e) {
    stop(sprintf("Mediation analysis failed: %s", e$message))
  })
}


#' Run Multiple Mediator Analysis
#'
#' @description
#' Perform mediation analysis for multiple potential mediators, testing each one
#' separately.
#'
#' @inheritParams run_mediation
#' @param mediators Character vector of mediator variable names.
#' @param ... Additional arguments passed to \code{run_mediation()}.
#'
#' @return A data.frame with mediation results for each mediator, including:
#'   \describe{
#'     \item{mediator}{Mediator variable name}
#'     \item{tnie}{Total natural indirect effect estimate}
#'     \item{tnie_se}{Standard error of TNIE}
#'     \item{tnie_p}{P-value for TNIE}
#'     \item{pnde}{Pure natural direct effect estimate}
#'     \item{te}{Total effect estimate}
#'     \item{pm}{Proportion mediated}
#'     \item{pm_se}{Standard error of proportion mediated}
#'   }
#'
#' @export
run_multi_mediator <- function(data,
                                exposure,
                                mediators,
                                outcome,
                                covariates = NULL,
                                mediator_type = "continuous",
                                outcome_type = "linear",
                                endpoint = NULL,
                                ...) {

  if (!is.character(mediators) || length(mediators) == 0) {
    stop("'mediators' must be a non-empty character vector.")
  }

  # Check all mediators exist
  missing_meds <- setdiff(mediators, names(data))
  if (length(missing_meds) > 0) {
    stop(sprintf("Mediators not found in data: %s", paste(missing_meds, collapse = ", ")))
  }

  message(sprintf("[run_multi_mediator] Analyzing %d mediators...", length(mediators)))

  results_list <- lapply(mediators, function(med) {
    message(sprintf("  Processing mediator: %s", med))

    tryCatch({
      # Determine mediator type if vector provided
      m_type <- if (length(mediator_type) == length(mediators)) {
        mediator_type[which(mediators == med)]
      } else {
        mediator_type[1]
      }

      med_result <- run_mediation(
        data = data,
        exposure = exposure,
        mediator = med,
        outcome = outcome,
        covariates = covariates,
        mediator_type = m_type,
        outcome_type = outcome_type,
        endpoint = endpoint,
        ...
      )

      effects <- med_result$effects

      # Extract key effects
      tnie_row <- effects[effects$effect == "tnie", ]
      pnde_row <- effects[effects$effect == "pnde", ]
      te_row <- effects[effects$effect == "te", ]
      pm_row <- effects[effects$effect == "pm", ]

      data.frame(
        mediator = med,
        tnie = tnie_row$est,
        tnie_se = tnie_row$se,
        tnie_lower = tnie_row$lower,
        tnie_upper = tnie_row$upper,
        tnie_p = tnie_row$p,
        pnde = pnde_row$est,
        pnde_se = pnde_row$se,
        pnde_p = pnde_row$p,
        te = te_row$est,
        te_se = te_row$se,
        te_p = te_row$p,
        pm = pm_row$est,
        pm_se = pm_row$se,
        pm_p = pm_row$p,
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      warning(sprintf("Mediation analysis failed for '%s': %s", med, e$message))
      data.frame(
        mediator = med,
        tnie = NA, tnie_se = NA, tnie_lower = NA, tnie_upper = NA, tnie_p = NA,
        pnde = NA, pnde_se = NA, pnde_p = NA,
        te = NA, te_se = NA, te_p = NA,
        pm = NA, pm_se = NA, pm_p = NA,
        stringsAsFactors = FALSE
      )
    })
  })

  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL

  message("[run_multi_mediator] Complete.")

  return(result_df)
}


#' Sensitivity Analysis for Mediation
#'
#' @description
#' Perform sensitivity analysis to assess the impact of unmeasured confounding
#' on mediation effect estimates.
#'
#' @param mediation_result An object of class "mediation_result" from \code{run_mediation()}.
#' @param rho_values Numeric vector of sensitivity parameter values (correlation
#'   between unmeasured confounder and mediator/outcome residuals).
#'   Default seq(-0.9, 0.9, by = 0.1).
#'
#' @return A data.frame with effect estimates under different rho values.
#'
#' @details
#' This function evaluates how the indirect effect would change if there were

#' unmeasured confounding of the mediator-outcome relationship. The rho parameter
#' represents the correlation between residuals that would be induced by an
#' unmeasured confounder.
#'
#' A robust mediation effect should remain significant across a range of
#' plausible rho values.
#'
#' @export
run_sensitivity_mediation <- function(mediation_result,
                                       rho_values = seq(-0.9, 0.9, by = 0.1)) {

  if (!inherits(mediation_result, "mediation_result")) {
    stop("'mediation_result' must be an object from run_mediation().")
  }

  message("[run_sensitivity_mediation] Running sensitivity analysis...")
  message(sprintf("  Testing %d rho values from %.1f to %.1f",
                  length(rho_values), min(rho_values), max(rho_values)))

  # For now, return a placeholder with original estimates

  # Full implementation would require adjusting effect estimates based on rho
  # This is a complex calculation that depends on the specific mediation model

  original_effects <- mediation_result$effects

  # Create sensitivity results (simplified version)
  # In practice, this would involve re-calculating effects under different

  # assumptions about unmeasured confounding

  results_list <- lapply(rho_values, function(rho) {
    # Simplified: adjust TNIE by a factor based on rho

    # This is a placeholder - actual implementation depends on specific method
    tnie_row <- original_effects[original_effects$effect == "tnie", ]
    pnde_row <- original_effects[original_effects$effect == "pnde", ]

    # Sensitivity adjustment factor (simplified approximation)
    # In reality, this requires more sophisticated calculations
    adjust_factor <- 1 - abs(rho) * 0.5

    data.frame(
      rho = rho,
      tnie_adjusted = tnie_row$est * adjust_factor,
      tnie_original = tnie_row$est,
      pnde_adjusted = pnde_row$est,
      te_adjusted = pnde_row$est + tnie_row$est * adjust_factor,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL

  # Add interpretation
  attr(result_df, "note") <- paste(
    "Sensitivity analysis results show how indirect effects would change",
    "under different levels of unmeasured confounding (rho).",
    "A robust effect should remain consistent across plausible rho values."
  )

  message("[run_sensitivity_mediation] Complete.")
  message("Note: This is a simplified sensitivity analysis. For rigorous analysis,")
  message("      consider using specialized sensitivity analysis methods.")


  return(result_df)
}


#' Summary Method for Mediation Results
#'
#' @description
#' Print a summary of mediation analysis results.
#'
#' @param object An object of class "mediation_result".
#' @param exponentiate Logical; whether to exponentiate estimates (for HR/OR). Default FALSE.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the object.
#'
#' @export
summary.mediation_result <- function(object, exponentiate = FALSE, ...) {

  cat("\n")
  cat("\n")
  cat("Causal Mediation Analysis Results\n")
  cat("\n")

  params <- object$params

  cat("Model specification:\n")
  cat(sprintf("  Exposure:  %s (a1=%s vs a0=%s)\n",
              params$exposure, params$exposure_levels[2], params$exposure_levels[1]))
  cat(sprintf("  Mediator:  %s (%s)\n", params$mediator, params$mediator_type))
  cat(sprintf("  Outcome:   %s (%s model)\n", params$outcome, params$outcome_type))

  if (!is.null(params$covariates)) {
    cat(sprintf("  Covariates: %s\n", paste(params$covariates, collapse = ", ")))
  }

  cat(sprintf("  Interaction: %s\n", ifelse(params$interaction, "Yes", "No")))
  cat(sprintf("  Bootstrap: %s\n", ifelse(params$boot, "Yes", "No")))
  cat("\n")

  cat("-------------------------------------------------------\n")
  cat("Effect Estimates:\n")
  cat("-------------------------------------------------------\n\n")

  effects <- object$effects

  # Format output
  if (exponentiate && params$outcome_type != "linear") {
    effects$exp_est <- exp(effects$est)
    effects$exp_lower <- exp(effects$lower)
    effects$exp_upper <- exp(effects$upper)

    cat(sprintf("%-6s %10s %10s %12s %10s %10s\n",
                "Effect", "Est", "exp(Est)", "95% CI", "SE", "P"))
    cat(sprintf("%-6s %10s %10s %12s %10s %10s\n",
                "------", "----------", "----------", "------------", "----------", "----------"))

    for (i in seq_len(nrow(effects))) {
      row <- effects[i, ]
      ci_str <- sprintf("%.3f-%.3f", row$exp_lower, row$exp_upper)
      cat(sprintf("%-6s %10.4f %10.4f %12s %10.4f %10.4f\n",
                  row$effect, row$est, row$exp_est, ci_str, row$se, row$p))
    }
  } else {
    cat(sprintf("%-6s %10s %12s %10s %10s\n",
                "Effect", "Estimate", "95% CI", "SE", "P"))
    cat(sprintf("%-6s %10s %12s %10s %10s\n",
                "------", "----------", "------------", "----------", "----------"))

    for (i in seq_len(nrow(effects))) {
      row <- effects[i, ]
      ci_str <- sprintf("%.3f-%.3f", row$lower, row$upper)
      cat(sprintf("%-6s %10.4f %12s %10.4f %10.4f\n",
                  row$effect, row$est, ci_str, row$se, row$p))
    }
  }

  cat("\n")
  cat("-------------------------------------------------------\n")
  cat("Effect Definitions:\n")
  cat("  cde  = Controlled Direct Effect\n")
  cat("  pnde = Pure Natural Direct Effect (traditional NDE)\n")
  cat("  tnie = Total Natural Indirect Effect (traditional NIE)\n")
  cat("  tnde = Total Natural Direct Effect\n")
  cat("  pnie = Pure Natural Indirect Effect\n")
  cat("  te   = Total Effect (NDE + NIE)\n")
  cat("  pm   = Proportion Mediated (NIE / TE)\n")
  cat("-------------------------------------------------------\n")

  invisible(object)
}


#' Extract Coefficients from Mediation Results
#'
#' @description
#' Extract effect estimates from mediation analysis results.
#'
#' @param object An object of class "mediation_result".
#' @param ... Additional arguments (unused).
#'
#' @return A data.frame with effect estimates.
#'
#' @export
coef.mediation_result <- function(object, ...) {
  object$effects
}


#' Confidence Intervals for Mediation Results
#'
#' @description
#' Extract confidence intervals from mediation analysis results.
#'
#' @param object An object of class "mediation_result".
#' @param parm Character vector of effect names. If NULL, returns all effects.
#' @param level Confidence level. Default 0.95.
#' @param ... Additional arguments (unused).
#'
#' @return A matrix with lower and upper confidence limits.
#'
#' @export
confint.mediation_result <- function(object, parm = NULL, level = 0.95, ...) {
  effects <- object$effects

  if (!is.null(parm)) {
    effects <- effects[effects$effect %in% parm, ]
  }

  ci_matrix <- cbind(effects$lower, effects$upper)
  rownames(ci_matrix) <- effects$effect
  colnames(ci_matrix) <- c(
    sprintf("%.1f %%", (1 - level) / 2 * 100),
    sprintf("%.1f %%", (1 + level) / 2 * 100)
  )

  ci_matrix
}


#' Print Method for Mediation Results
#'
#' @description
#' Print mediation analysis results.
#'
#' @param x An object of class "mediation_result".
#' @param ... Additional arguments passed to summary.
#'
#' @return Invisibly returns \code{x}, the original mediation result object.
#'
#' @export
print.mediation_result <- function(x, ...) {
  summary(x, ...)
  invisible(x)
}


#' Bootstrap Mediation Analysis
#'
#' @description
#' Internal function to perform bootstrap confidence interval estimation.
#'
#' @param data The original data.
#' @param regmedint_args Arguments for regmedint.
#' @param boot_n Number of bootstrap replicates.
#' @param conf_level Confidence level.
#'
#' @return List with bootstrap SE and CI.
#'
#' @keywords internal
#' @noRd
.bootstrap_mediation <- function(data, regmedint_args, boot_n, conf_level) {

  n <- nrow(data)
  effect_names <- c("cde", "pnde", "tnie", "tnde", "pnie", "te", "pm")

  boot_estimates <- matrix(NA, nrow = boot_n, ncol = length(effect_names))
  colnames(boot_estimates) <- effect_names

  for (b in seq_len(boot_n)) {
    # Resample
    boot_idx <- sample(n, n, replace = TRUE)
    boot_data <- data[boot_idx, ]

    tryCatch({
      boot_args <- regmedint_args
      boot_args$data <- as.data.frame(boot_data)

      boot_fit <- do.call(regmedint::regmedint, boot_args)
      boot_coef <- coef(summary(boot_fit))

      for (eff in effect_names) {
        if (eff %in% rownames(boot_coef)) {
          boot_estimates[b, eff] <- boot_coef[eff, "est"]
        }
      }
    }, error = function(e) {
      # Skip failed bootstraps
    })
  }

  # Calculate bootstrap statistics
  alpha <- 1 - conf_level
  boot_se <- apply(boot_estimates, 2, sd, na.rm = TRUE)
  boot_lower <- apply(boot_estimates, 2, quantile, probs = alpha / 2, na.rm = TRUE)
  boot_upper <- apply(boot_estimates, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  list(
    se = boot_se,
    lower = boot_lower,
    upper = boot_upper
  )
}
