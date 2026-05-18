#' @importFrom stats relevel anova pchisq
#' @importFrom survival Surv coxph
NULL

#' Run Subgroup Analysis
#'
#' @description
#' Perform subgroup analysis by fitting regression models within each level of a
#' subgroup variable and calculating interaction p-values.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param exposure Character string specifying the exposure variable name.
#' @param outcome Character string specifying the outcome variable name.
#'   For Cox models, this can be NULL if endpoint is specified.
#' @param subgroup_var Character string specifying the subgroup variable name.
#' @param covariates Character vector of covariate names to adjust for. Default NULL.
#' @param model_type Character string specifying model type: \code{"cox"},
#'   \code{"logistic"}, \code{"linear"}, \code{"glm"}, or \code{"negbin"}.
#' @param family For \code{model_type = "glm"} only: the GLM family.  Accepts
#'   a character string, function, or family object (see
#'   \code{\link{runmulti_glm}}).  Default \code{"poisson"}.
#' @param endpoint Character vector of length 2 for Cox models: c("time", "status").
#'   Required when model_type = "cox".
#' @param ref_level Character string specifying the reference level for the subgroup variable.
#'   If NULL, the first level is used as reference.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{subgroup_var}{Name of the subgroup variable}
#'     \item{subgroup}{Subgroup level}
#'     \item{n}{Sample size in subgroup}
#'     \item{n_event}{Number of events (for Cox/logistic models)}
#'     \item{estimate}{Effect estimate (HR for Cox, OR for logistic, Beta for linear)}
#'     \item{lower95}{Lower 95\\% CI}
#'     \item{upper95}{Upper 95\\% CI}
#'     \item{pvalue}{P-value for the exposure effect}
#'     \item{p_interaction}{P-value for interaction between exposure and subgroup}
#'   }
#'
#' @examples
#' \dontrun{
#' library(survival)
#' # Cox model subgroup analysis
#' result <- run_subgroup_analysis(
#'   data = lung,
#'   exposure = "ph.ecog",
#'   subgroup_var = "sex",
#'   model_type = "cox",
#'   endpoint = c("time", "status")
#' )
#'
#' # Logistic regression subgroup analysis
#' mtcars$am_binary <- ifelse(mtcars$am == 1, 1, 0)
#' result <- run_subgroup_analysis(
#'   data = mtcars,
#'   exposure = "hp",
#'   outcome = "am_binary",
#'   subgroup_var = "cyl",
#'   model_type = "logistic"
#' )
#' }
#'
#' @importFrom stats as.formula glm binomial lm coef confint qnorm
#' @export
run_subgroup_analysis <- function(data,
                                   exposure,
                                   outcome = NULL,
                                   subgroup_var,
                                   covariates = NULL,
                                   model_type = c("cox", "logistic", "linear",
                                                  "glm", "negbin"),
                                   family = "poisson",
                                   endpoint = NULL,
                                   ref_level = NULL) {

  model_type <- match.arg(model_type)

  # Parse GLM family (only used when model_type == "glm")
  family_obj <- if (model_type == "glm") .parse_glm_family(family) else NULL
  is_quasi   <- !is.null(family_obj) && grepl("^quasi", family_obj$family)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")

}

  if (!is.character(exposure) || length(exposure) != 1) {
    stop("'exposure' must be a single character string.")
  }

  if (!is.character(subgroup_var) || length(subgroup_var) != 1) {
    stop("'subgroup_var' must be a single character string.")
  }

  # Validate model-specific requirements
  if (model_type == "cox") {
    if (is.null(endpoint) || length(endpoint) != 2) {
      stop("For Cox models, 'endpoint' must be a character vector of length 2: c('time', 'status').")
    }
    required_cols <- c(exposure, subgroup_var, covariates, endpoint)
  } else {
    if (is.null(outcome)) {
      stop(sprintf(
        "'outcome' is required for %s models.", model_type
      ))
    }
    required_cols <- c(exposure, outcome, subgroup_var, covariates)
  }

  # Check all variables exist
  missing_vars <- setdiff(required_cols, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  # Convert to data.table for efficiency
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  # Ensure subgroup_var is a factor
  if (!is.factor(data[[subgroup_var]])) {
    data[[subgroup_var]] <- as.factor(data[[subgroup_var]])
  }

  # Set reference level if specified
  if (!is.null(ref_level)) {
    if (!ref_level %in% levels(data[[subgroup_var]])) {
      stop(sprintf("ref_level '%s' not found in subgroup_var levels.", ref_level))
    }
    data[[subgroup_var]] <- relevel(data[[subgroup_var]], ref = ref_level)
  }

  # Calculate interaction p-value
  p_interaction <- .calculate_interaction_pvalue(
    data = data,
    exposure = exposure,
    outcome = outcome,
    subgroup_var = subgroup_var,
    covariates = covariates,
    model_type = model_type,
    endpoint = endpoint,
    family_obj = family_obj,
    is_quasi = is_quasi
  )

  # Get subgroup levels
  subgroup_levels <- levels(data[[subgroup_var]])

  # Run analysis for each subgroup
  results_list <- lapply(subgroup_levels, function(level) {
    # Subset data
    subset_data <- data[data[[subgroup_var]] == level, ]

    # Get sample size
    n <- nrow(subset_data)

    # Get number of events
    n_event <- NA
    if (model_type == "cox") {
      n_event <- sum(subset_data[[endpoint[2]]] == 1, na.rm = TRUE)
    } else if (model_type %in% c("logistic", "glm", "negbin")) {
      n_event <- sum(subset_data[[outcome]] == 1, na.rm = TRUE)
    }

    # Exposure must vary within subgroup; otherwise effect is not estimable.
    exposure_values <- subset_data[[exposure]]
    n_unique_exposure <- length(unique(exposure_values[!is.na(exposure_values)]))
    if (n_unique_exposure < 2) {
      warning(sprintf(
        "Exposure '%s' has no variation in subgroup '%s'; returning NA for effect estimate.",
        exposure, level
      ))
      return(.na_subgroup_row(subgroup_var, level, n, n_event, p_interaction))
    }

    # Build formula
    if (model_type == "cox") {
      rhs <- exposure
      if (!is.null(covariates)) {
        rhs <- paste(c(exposure, covariates), collapse = " + ")
      }
      formula_str <- paste0("Surv(", endpoint[1], ", ", endpoint[2], ") ~ ", rhs)
    } else {
      rhs <- exposure
      if (!is.null(covariates)) {
        rhs <- paste(c(exposure, covariates), collapse = " + ")
      }
      formula_str <- paste(outcome, "~", rhs)
    }
    formula_obj <- stats::as.formula(formula_str)

    # Fit model
    tryCatch({
      if (model_type == "cox") {
        model <- coxph(formula_obj, data = subset_data)
        sum_model <- summary(model)
        coefs <- sum_model$coefficients
        conf <- sum_model$conf.int

        # Find the row for exposure using model-term mapping to avoid prefix mis-match.
        exp_row <- .find_exposure_row(model, coefs, exposure)
        if (is.na(exp_row)) {
          return(.na_subgroup_row(subgroup_var, level, n, n_event,
                                  p_interaction))
        }

        estimate <- conf[exp_row, "exp(coef)"]
        lower95 <- conf[exp_row, "lower .95"]
        upper95 <- conf[exp_row, "upper .95"]
        pvalue <- coefs[exp_row, "Pr(>|z|)"]

      } else if (model_type == "logistic") {
        model <- stats::glm(formula_obj, data = subset_data, family = stats::binomial())
        sum_model <- summary(model)
        coefs <- sum_model$coefficients
        ci <- suppressWarnings(stats::confint(model))

        exp_row <- .find_exposure_row(model, coefs, exposure)
        if (is.na(exp_row)) {
          return(.na_subgroup_row(subgroup_var, level, n, n_event,
                                  p_interaction))
        }

        estimate <- exp(coefs[exp_row, "Estimate"])
        lower95 <- exp(ci[exp_row, 1])
        upper95 <- exp(ci[exp_row, 2])
        pvalue <- coefs[exp_row, "Pr(>|z|)"]

      } else if (model_type == "glm") {
        model     <- stats::glm(formula_obj, data = subset_data,
                                family = family_obj)
        sum_model <- summary(model)
        coefs     <- sum_model$coefficients
        exp_row   <- .find_exposure_row(model, coefs, exposure)

        if (is.na(exp_row)) {
          return(.na_subgroup_row(subgroup_var, level, n, n_event,
                                  p_interaction))
        }

        b  <- coefs[exp_row, "Estimate"]
        se <- coefs[exp_row, "Std. Error"]

        if (is_quasi) {
          z95 <- stats::qnorm(0.975)
          ci_lo <- b - z95 * se
          ci_hi <- b + z95 * se
        } else {
          ci_mat <- suppressWarnings(stats::confint(model))
          ci_lo  <- ci_mat[exp_row, 1]
          ci_hi  <- ci_mat[exp_row, 2]
        }

        exp_scale <- family_obj$link %in% c("log", "logit", "cloglog")
        estimate  <- if (exp_scale) exp(b)    else b
        lower95   <- if (exp_scale) exp(ci_lo) else ci_lo
        upper95   <- if (exp_scale) exp(ci_hi) else ci_hi
        pval_col  <- if ("Pr(>|t|)" %in% colnames(coefs)) "Pr(>|t|)" else "Pr(>|z|)"
        pvalue    <- coefs[exp_row, pval_col]

      } else if (model_type == "negbin") {
        model     <- MASS::glm.nb(formula_obj, data = subset_data)
        sum_model <- summary(model)
        coefs     <- sum_model$coefficients
        exp_row   <- .find_exposure_row(model, coefs, exposure)

        if (is.na(exp_row)) {
          return(.na_subgroup_row(subgroup_var, level, n, n_event,
                                  p_interaction))
        }

        ci_mat   <- suppressWarnings(stats::confint(model))
        estimate <- exp(coefs[exp_row, "Estimate"])
        lower95  <- exp(ci_mat[exp_row, 1])
        upper95  <- exp(ci_mat[exp_row, 2])
        pvalue   <- coefs[exp_row, "Pr(>|z|)"]

      } else {
        # Linear model
        model <- stats::lm(formula_obj, data = subset_data)
        sum_model <- summary(model)
        coefs <- sum_model$coefficients
        ci <- stats::confint(model)

        exp_row <- .find_exposure_row(model, coefs, exposure)
        if (is.na(exp_row)) {
          return(.na_subgroup_row(subgroup_var, level, n, NA, p_interaction))
        }

        estimate <- coefs[exp_row, "Estimate"]
        lower95 <- ci[exp_row, 1]
        upper95 <- ci[exp_row, 2]
        pvalue <- coefs[exp_row, "Pr(>|t|)"]
      }

      data.frame(
        subgroup_var = subgroup_var,
        subgroup = level,
        n = n,
        n_event = n_event,
        estimate = round(estimate, 3),
        lower95 = round(lower95, 3),
        upper95 = round(upper95, 3),
        pvalue = signif(pvalue, 3),
        p_interaction = signif(p_interaction, 3),
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      warning(sprintf("Model fitting failed for subgroup '%s': %s",
                      level, e$message))
      .na_subgroup_row(subgroup_var, level, n, n_event, p_interaction)
    })
  })

  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL
  return(result_df)
}

#' @keywords internal
#' @noRd
.find_exposure_row <- function(model, coefs, exposure) {
  coef_names <- rownames(coefs)
  if (is.null(coef_names) || length(coef_names) == 0) {
    return(NA_integer_)
  }

  # 1) Prefer exact coefficient name match (continuous/numeric exposure)
  exact <- which(coef_names == exposure)
  if (length(exact) > 0) {
    return(exact[1])
  }

  # 2) Use model-matrix assign mapping (safe for factor exposure)
  mm <- tryCatch(stats::model.matrix(model), error = function(e) NULL)
  if (!is.null(mm)) {
    assign <- attr(mm, "assign")
    term_labels <- attr(stats::terms(model), "term.labels")
    exposure_terms <- c(exposure, paste0("`", exposure, "`"))
    exp_term_idx <- which(term_labels %in% exposure_terms)

    if (length(exp_term_idx) >= 1) {
      # If repeated terms appear, use the first exact term occurrence.
      exp_cols <- colnames(mm)[assign == exp_term_idx[1]]
      exp_rows <- match(exp_cols, coef_names)
      exp_rows <- exp_rows[!is.na(exp_rows)]
      if (length(exp_rows) > 0) {
        return(exp_rows[1])
      }
    }
  }

  # Do not use prefix regex fallback here to avoid x -> x2 mis-match.
  NA_integer_
}


#' Run Multiple Subgroup Analyses
#'
#' @description
#' Perform subgroup analyses across multiple subgroup variables.
#'
#' @inheritParams run_subgroup_analysis
#' @param subgroup_vars Character vector of subgroup variable names.
#'
#' @return A data.frame with results from all subgroup analyses combined.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' result <- run_multi_subgroup(
#'   data = lung,
#'   exposure = "ph.ecog",
#'   subgroup_vars = c("sex"),
#'   model_type = "cox",
#'   endpoint = c("time", "status")
#' )
#' }
#'
#' @export
run_multi_subgroup <- function(data,
                                exposure,
                                outcome = NULL,
                                subgroup_vars,
                                covariates = NULL,
                                model_type = c("cox", "logistic", "linear",
                                               "glm", "negbin"),
                                family = "poisson",
                                endpoint = NULL) {

  model_type <- match.arg(model_type)

  if (!is.character(subgroup_vars) || length(subgroup_vars) == 0) {
    stop("'subgroup_vars' must be a non-empty character vector.")
  }

  # Run subgroup analysis for each variable
  results_list <- lapply(subgroup_vars, function(sv) {
    message(sprintf("[run_multi_subgroup] Analyzing subgroup: %s", sv))
    run_subgroup_analysis(
      data = data,
      exposure = exposure,
      outcome = outcome,
      subgroup_var = sv,
      covariates = covariates,
      model_type = model_type,
      family = family,
      endpoint = endpoint
    )
  })

  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Calculate Interaction P-value
#'
#' @description
#' Internal function to calculate the p-value for the interaction between
#' exposure and subgroup variable.
#'
#' @inheritParams run_subgroup_analysis
#'
#' @return Numeric p-value for the interaction term.
#'
#' @keywords internal
#' @noRd
.calculate_interaction_pvalue <- function(data,
                                           exposure,
                                           outcome,
                                           subgroup_var,
                                           covariates,
                                           model_type,
                                           endpoint,
                                           family_obj = NULL,
                                           is_quasi   = FALSE) {

  .build_int_rhs <- function(int_term) {
    if (!is.null(covariates))
      paste(c(int_term, covariates), collapse = " + ")
    else
      int_term
  }
  .build_no_int_rhs <- function() {
    paste(c(exposure, subgroup_var, covariates), collapse = " + ")
  }
  int_term <- paste0(exposure, " * ", subgroup_var)

  tryCatch({
    # ── fit full (interaction) model ──────────────────────────────────────────
    if (model_type == "cox") {
      formula_obj <- stats::as.formula(paste0(
        "Surv(", endpoint[1], ", ", endpoint[2], ") ~ ",
        .build_int_rhs(int_term)
      ))
      model <- coxph(formula_obj, data = data)

    } else if (model_type == "logistic") {
      formula_obj <- stats::as.formula(
        paste(outcome, "~", .build_int_rhs(int_term))
      )
      model <- stats::glm(formula_obj, data = data, family = stats::binomial())

    } else if (model_type == "glm") {
      formula_obj <- stats::as.formula(
        paste(outcome, "~", .build_int_rhs(int_term))
      )
      model <- stats::glm(formula_obj, data = data, family = family_obj)

    } else if (model_type == "negbin") {
      formula_obj <- stats::as.formula(
        paste(outcome, "~", .build_int_rhs(int_term))
      )
      model <- MASS::glm.nb(formula_obj, data = data)

    } else {
      formula_obj <- stats::as.formula(
        paste(outcome, "~", .build_int_rhs(int_term))
      )
      model <- stats::lm(formula_obj, data = data)
    }

    # ── locate interaction coefficient rows ───────────────────────────────────
    coefs <- summary(model)$coefficients
    pat   <- paste0(exposure, ":", subgroup_var, "|",
                    subgroup_var, ":", exposure)
    interaction_rows <- grep(pat, rownames(coefs))

    if (length(interaction_rows) == 0) return(NA)

    # ── single interaction term → Wald p-value ────────────────────────────────
    if (length(interaction_rows) == 1) {
      pval_col <- if ("Pr(>|t|)" %in% colnames(coefs)) "Pr(>|t|)" else "Pr(>|z|)"
      return(coefs[interaction_rows[1], pval_col])
    }

    # ── multiple terms → LRT via anova ───────────────────────────────────────
    rhs_no_int <- .build_no_int_rhs()

    if (model_type == "cox") {
      formula_no_int <- stats::as.formula(paste0(
        "Surv(", endpoint[1], ", ", endpoint[2], ") ~ ", rhs_no_int
      ))
      model_no_int <- coxph(formula_no_int, data = data)
      av <- stats::anova(model_no_int, model)

    } else if (model_type == "logistic") {
      formula_no_int <- stats::as.formula(paste(outcome, "~", rhs_no_int))
      model_no_int   <- stats::glm(formula_no_int, data = data,
                                   family = stats::binomial())
      av <- stats::anova(model_no_int, model, test = "Chisq")

    } else if (model_type == "glm") {
      formula_no_int <- stats::as.formula(paste(outcome, "~", rhs_no_int))
      model_no_int   <- stats::glm(formula_no_int, data = data,
                                   family = family_obj)
      test_type <- if (is_quasi) "F" else "Chisq"
      av <- stats::anova(model_no_int, model, test = test_type)

    } else if (model_type == "negbin") {
      formula_no_int <- stats::as.formula(paste(outcome, "~", rhs_no_int))
      model_no_int   <- MASS::glm.nb(formula_no_int, data = data)
      av <- stats::anova(model_no_int, model, test = "Chisq")

    } else {
      formula_no_int <- stats::as.formula(paste(outcome, "~", rhs_no_int))
      model_no_int   <- stats::lm(formula_no_int, data = data)
      av <- stats::anova(model_no_int, model)
    }

    # Robustly extract p-value regardless of column name differences
    pval_candidates <- c("Pr(>|Chi|)", "Pr(>Chi)", "Pr(>F)", "Pr(Chi)")
    av_df <- as.data.frame(av)
    for (col in pval_candidates) {
      if (col %in% names(av_df)) return(av_df[[col]][2])
    }
    NA

  }, error = function(e) {
    warning(sprintf("Failed to calculate interaction p-value: %s", e$message))
    NA
  })
}


#' Build a one-row NA result for a subgroup
#' @keywords internal
#' @noRd
.na_subgroup_row <- function(subgroup_var, level, n, n_event, p_interaction) {
  data.frame(
    subgroup_var  = subgroup_var,
    subgroup      = level,
    n             = n,
    n_event       = n_event,
    estimate      = NA_real_,
    lower95       = NA_real_,
    upper95       = NA_real_,
    pvalue        = NA_real_,
    p_interaction = p_interaction,
    stringsAsFactors = FALSE
  )
}
