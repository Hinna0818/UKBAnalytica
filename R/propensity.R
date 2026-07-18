#' Estimate Propensity Score
#'
#' @description
#' Calculate propensity scores using logistic regression or gradient boosting.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param treatment Character string specifying the treatment variable name (binary 0/1).
#' @param covariates Character vector of covariate names used to estimate propensity scores.
#' @param method Character string specifying the estimation method: "logistic" (default) or "gbm".
#' @param formula Optional custom formula. If NULL, formula is built from treatment and covariates.
#'
#' @return A data.table with the original data plus:
#'   \describe{
#'     \item{ps}{Propensity score (probability of treatment)}
#'   }
#'
#' @import data.table
#' @importFrom stats glm binomial predict as.formula
#' @importFrom survival Surv coxph
#' @importFrom sandwich vcovHC
#' @importFrom lmtest coeftest
#' @export
estimate_propensity_score <- function(data,
                                       treatment,
                                       covariates,
                                       method = c("logistic", "gbm"),
                                       formula = NULL) {

  method <- match.arg(method)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!is.character(treatment) || length(treatment) != 1) {
    stop("'treatment' must be a single character string.")
  }

  if (!treatment %in% names(data)) {
    stop(sprintf("Treatment variable '%s' not found in data.", treatment))
  }

  # Check treatment is binary
  treat_vals <- data[[treatment]]
  treat_vals <- treat_vals[!is.na(treat_vals)]
  unique_vals <- sort(unique(treat_vals))
  if (!all(unique_vals %in% c(0, 1))) {
    stop(sprintf("Treatment variable '%s' must be binary (0/1). Found values: %s",
                 treatment, paste(unique_vals, collapse = ", ")))
  }

  if (!is.character(covariates) || length(covariates) == 0) {
    stop("'covariates' must be a non-empty character vector.")
  }

  missing_vars <- setdiff(covariates, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Covariates not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  # Convert to data.table
  if (!data.table::is.data.table(data)) {
    dt <- data.table::as.data.table(data)
  } else {
    dt <- data.table::copy(data)
  }

  message("[estimate_propensity_score] Estimating propensity scores...")
  message(sprintf("  Method: %s", method))
  message(sprintf("  Treatment: %s", treatment))
  message(sprintf("  Covariates: %s", paste(covariates, collapse = ", ")))

  # Build formula if not provided
  if (is.null(formula)) {
    formula <- stats::as.formula(paste(treatment, "~", paste(covariates, collapse = " + ")))
  }
  model_vars <- unique(all.vars(formula))
  missing_formula_vars <- setdiff(model_vars, names(dt))
  if (length(missing_formula_vars) > 0L) {
    stop(
      "Variables in `formula` were not found in data: ",
      paste(missing_formula_vars, collapse = ", "),
      call. = FALSE
    )
  }
  complete_idx <- stats::complete.cases(dt[, model_vars, with = FALSE])
  fit_rows <- which(complete_idx)
  if (length(fit_rows) == 0L) {
    stop("No complete rows are available for propensity-score estimation.", call. = FALSE)
  }
  fit_data <- as.data.frame(dt[fit_rows])
  if (length(unique(fit_data[[treatment]])) != 2L) {
    stop(
      "Both treatment groups must be present after complete-case filtering.",
      call. = FALSE
    )
  }

  if (method == "logistic") {
    # Logistic regression
    model <- stats::glm(formula, data = fit_data, family = stats::binomial())
    fitted_ps <- stats::predict(model, newdata = fit_data, type = "response")

  } else if (method == "gbm") {
    # Gradient boosting
    if (!requireNamespace("gbm", quietly = TRUE)) {
      stop("Package 'gbm' is required for method = 'gbm'. Please install it.")
    }

    # GBM requires numeric response
    model <- gbm::gbm(
      formula,
      data = fit_data,
      distribution = "bernoulli",
      n.trees = 1000,
      interaction.depth = 3,
      shrinkage = 0.01,
      bag.fraction = 0.5,
      cv.folds = 5,
      verbose = FALSE
    )

    best_iter <- gbm::gbm.perf(model, method = "cv", plot.it = FALSE)
    fitted_ps <- gbm::predict.gbm(
      model,
      newdata = fit_data,
      n.trees = best_iter,
      type = "response"
    )
  }

  # Ensure PS is within (0, 1)
  ps <- rep(NA_real_, nrow(dt))
  ps[fit_rows] <- pmax(0.001, pmin(0.999, as.numeric(fitted_ps)))
  dt[, "ps" := ps]
  attr(dt, "propensity_model") <- model
  attr(dt, "ps_complete_rows") <- fit_rows

  message(sprintf(
    "[estimate_propensity_score] Complete. Evaluated %d/%d rows; PS range: [%.3f, %.3f]",
    length(fit_rows), nrow(dt), min(fitted_ps), max(fitted_ps)
  ))

  return(dt)
}


#' Propensity Score Matching
#'
#' @description
#' Perform propensity score matching using nearest neighbor or optimal matching.
#'
#' @param data A data.table containing propensity scores.
#' @param ps_col Character string specifying the propensity score column name. Default "ps".
#' @param treatment Character string specifying the treatment variable name.
#' @param ratio Numeric matching ratio (1:ratio). Default 1 for 1:1 matching.
#' @param caliper Numeric caliper width in standard deviations of PS. Default 0.2.
#' @param method Character string specifying matching method: "nearest" or "optimal".
#' @param replace Logical; whether to match with replacement. Default FALSE.
#' @param exact_match Character vector of variable names for exact matching. Default NULL.
#'
#' @return A data.table with matched data, including:
#'   \describe{
#'     \item{match_id}{Matched pair identifier}
#'     \item{match_distance}{Distance between matched pairs}
#'   }
#'
#' @import data.table
#' @export
match_propensity <- function(data,
                              ps_col = "ps",
                              treatment,
                              ratio = 1,
                              caliper = 0.2,
                              method = c("nearest", "optimal"),
                              replace = FALSE,
                              exact_match = NULL) {

  method <- match.arg(method)

  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' is required for propensity score matching. Please install it.")
  }

  # Validate inputs
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  if (!ps_col %in% names(data)) {
    stop(sprintf("Propensity score column '%s' not found in data.", ps_col))
  }

  if (!treatment %in% names(data)) {
    stop(sprintf("Treatment variable '%s' not found in data.", treatment))
  }

  message("[match_propensity] Performing propensity score matching...")
  message(sprintf("  Method: %s", method))
  message(sprintf("  Ratio: 1:%d", ratio))
  message(sprintf("  Caliper: %.2f SD", caliper))

  n_treated <- sum(data[[treatment]] == 1, na.rm = TRUE)
  n_control <- sum(data[[treatment]] == 0, na.rm = TRUE)
  message(sprintf("  Before matching: %d treated, %d control", n_treated, n_control))

  # Prepare MatchIt call
  # Build formula using PS
  formula_str <- paste(treatment, "~", ps_col)
  if (!is.null(exact_match)) {
    if (!all(exact_match %in% names(data))) {
      stop("Some exact_match variables not found in data.")
    }
  }

  # Convert method name for MatchIt
  matchit_method <- ifelse(method == "optimal", "optimal", "nearest")

  # Run matching
  match_obj <- MatchIt::matchit(
    formula = stats::as.formula(formula_str),
    data = as.data.frame(data),
    method = matchit_method,
    distance = data[[ps_col]],
    ratio = ratio,
    caliper = caliper,
    std.caliper = TRUE,
    replace = replace,
    exact = exact_match
  )

  # Extract matched data
  matched_df <- MatchIt::match.data(match_obj)
  matched_dt <- data.table::as.data.table(matched_df)

  # Rename MatchIt columns
  if ("subclass" %in% names(matched_dt)) {
    data.table::setnames(matched_dt, "subclass", "match_id")
  }
  if ("distance" %in% names(matched_dt)) {
    data.table::setnames(matched_dt, "distance", "match_distance")
  }
  attr(matched_dt, "matching_caliper") <- list(
    value = caliper,
    scale = "standard_deviations_of_distance"
  )

  n_matched_treated <- sum(matched_dt[[treatment]] == 1, na.rm = TRUE)
  n_matched_control <- sum(matched_dt[[treatment]] == 0, na.rm = TRUE)

  message(sprintf("[match_propensity] Complete. After matching: %d treated, %d control",
                  n_matched_treated, n_matched_control))

  return(matched_dt)
}


#' Calculate IPTW Weights
#'
#' @description
#' Calculate inverse probability of treatment weights (IPTW) for causal inference.
#'
#' @param data A data.table containing propensity scores.
#' @param ps_col Character string specifying the propensity score column name. Default "ps".
#' @param treatment Character string specifying the treatment variable name.
#' @param weight_type Character string specifying weight type: "ATE", "ATT", or "ATC".
#' @param stabilized Logical; whether to use stabilized weights. Default TRUE.
#' @param truncate Numeric vector of length 2 specifying quantiles for weight truncation.
#'   Default c(0.01, 0.99).
#'
#' @return A data.table with the original data plus:
#'   \describe{
#'     \item{weight}{IPTW weight}
#'   }
#'
#' @details
#' Weight formulas:
#' \itemize{
#'   \item ATE: \code{T/PS + (1-T)/(1-PS)}
#'   \item ATT: \code{T + (1-T) * PS/(1-PS)}
#'   \item ATC: \code{T * (1-PS)/PS + (1-T)}
#' }
#'
#' Stabilized weights multiply by the marginal probability of treatment.
#'
#' @import data.table
#' @export
calculate_weights <- function(data,
                               ps_col = "ps",
                               treatment,
                               weight_type = c("ATE", "ATT", "ATC"),
                               stabilized = TRUE,
                               truncate = c(0.01, 0.99)) {

  weight_type <- match.arg(weight_type)

  # Validate inputs
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  } else {
    data <- data.table::copy(data)
  }

  if (!ps_col %in% names(data)) {
    stop(sprintf("Propensity score column '%s' not found in data.", ps_col))
  }

  if (!treatment %in% names(data)) {
    stop(sprintf("Treatment variable '%s' not found in data.", treatment))
  }

  message(sprintf("[calculate_weights] Calculating %s weights...", weight_type))

  ps <- data[[ps_col]]
  treat <- data[[treatment]]

  # Calculate weights based on type
  if (weight_type == "ATE") {
    # Average Treatment Effect
    w <- treat / ps + (1 - treat) / (1 - ps)
  } else if (weight_type == "ATT") {
    # Average Treatment Effect on Treated
    w <- treat + (1 - treat) * ps / (1 - ps)
  } else {
    # Average Treatment Effect on Controls
    w <- treat * (1 - ps) / ps + (1 - treat)
  }

  # Stabilize weights
  if (stabilized) {
    p_treat <- mean(treat, na.rm = TRUE)
    if (weight_type == "ATE") {
      w <- treat * p_treat / ps + (1 - treat) * (1 - p_treat) / (1 - ps)
    } else if (weight_type == "ATT") {
      w <- treat + (1 - treat) * ps / (1 - ps)  # ATT stabilized = same
    } else {
      w <- treat * (1 - ps) / ps + (1 - treat)  # ATC stabilized = same
    }
  }

  # Truncate extreme weights
  if (!is.null(truncate) && length(truncate) == 2) {
    lower <- stats::quantile(w, truncate[1], na.rm = TRUE)
    upper <- stats::quantile(w, truncate[2], na.rm = TRUE)
    w[w < lower] <- lower
    w[w > upper] <- upper
    message(sprintf("  Weights truncated to [%.2f, %.2f]", lower, upper))
  }

  data[, weight := w]

  message(sprintf("[calculate_weights] Complete. Weight range: [%.2f, %.2f], mean: %.2f",
                  min(data$weight, na.rm = TRUE),
                  max(data$weight, na.rm = TRUE),
                  mean(data$weight, na.rm = TRUE)))

  return(data)
}


#' Assess Covariate Balance
#'
#' @description
#' Assess balance of covariates between treatment groups before and after
#' matching or weighting.
#'
#' @param data A data.frame or data.table.
#' @param treatment Character string specifying the treatment variable name.
#' @param covariates Character vector of covariate names to assess.
#' @param method Character string specifying the data type: "unmatched", "matched", or "weighted".
#' @param weight_col Character string specifying the weight column name (for weighted method).
#' @param threshold Numeric threshold for SMD to determine balance. Default 0.1.
#'
#' @return A data.frame with balance statistics:
#'   \describe{
#'     \item{variable}{Variable name}
#'     \item{mean_treated}{Mean in treatment group}
#'     \item{mean_control}{Mean in control group}
#'     \item{smd}{Standardized mean difference}
#'     \item{variance_ratio}{Variance ratio (treated/control)}
#'     \item{balanced}{Whether SMD < threshold}
#'   }
#'
#' @importFrom stats weighted.mean var
#' @export
assess_balance <- function(data,
                            treatment,
                            covariates,
                            method = c("unmatched", "matched", "weighted"),
                            weight_col = NULL,
                            threshold = 0.1) {

  method <- match.arg(method)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!treatment %in% names(data)) {
    stop(sprintf("Treatment variable '%s' not found in data.", treatment))
  }

  missing_vars <- setdiff(covariates, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Covariates not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  if (method == "weighted" && is.null(weight_col)) {
    stop("'weight_col' is required when method = 'weighted'.")
  }

  if (method == "weighted" && !weight_col %in% names(data)) {
    stop(sprintf("Weight column '%s' not found in data.", weight_col))
  }

  analysis_data <- as.data.frame(data)
  treat_values <- analysis_data[[treatment]]
  if (!all(stats::na.omit(unique(treat_values)) %in% c(0, 1))) {
    stop("`treatment` must be coded as 0/1.", call. = FALSE)
  }

  # Split data by treatment
  treated <- analysis_data[!is.na(treat_values) & treat_values == 1, , drop = FALSE]
  control <- analysis_data[!is.na(treat_values) & treat_values == 0, , drop = FALSE]

  calculate_one <- function(x_t, x_c, variable, level = NA_character_) {
    if (method == "weighted") {
      w_t <- treated[[weight_col]]
      w_c <- control[[weight_col]]
      mean_t <- stats::weighted.mean(x_t, w_t, na.rm = TRUE)
      mean_c <- stats::weighted.mean(x_c, w_c, na.rm = TRUE)
      # Weighted variance
      var_t <- sum(w_t * (x_t - mean_t)^2, na.rm = TRUE) / sum(w_t, na.rm = TRUE)
      var_c <- sum(w_c * (x_c - mean_c)^2, na.rm = TRUE) / sum(w_c, na.rm = TRUE)
    } else {
      mean_t <- mean(x_t, na.rm = TRUE)
      mean_c <- mean(x_c, na.rm = TRUE)
      var_t <- stats::var(x_t, na.rm = TRUE)
      var_c <- stats::var(x_c, na.rm = TRUE)
    }

    pooled_sd <- sqrt((var_t + var_c) / 2)
    smd <- if (is.finite(pooled_sd) && pooled_sd > 0) {
      (mean_t - mean_c) / pooled_sd
    } else if (isTRUE(all.equal(mean_t, mean_c))) {
      0
    } else {
      sign(mean_t - mean_c) * Inf
    }
    var_ratio <- if (is.finite(var_c) && var_c > 0) var_t / var_c else NA_real_

    data.frame(
      variable = variable,
      level = level,
      mean_treated = round(mean_t, 3),
      mean_control = round(mean_c, 3),
      smd = round(smd, 3),
      variance_ratio = round(var_ratio, 3),
      balanced = is.finite(smd) && abs(smd) < threshold,
      stringsAsFactors = FALSE
    )
  }

  # Nominal variables are represented by one binary indicator per level.
  results <- lapply(covariates, function(var) {
    x_t <- treated[[var]]
    x_c <- control[[var]]
    if (is.factor(analysis_data[[var]]) || is.character(analysis_data[[var]])) {
      levels_all <- sort(unique(c(as.character(x_t), as.character(x_c))))
      levels_all <- levels_all[!is.na(levels_all)]
      rows <- lapply(levels_all, function(level) {
        calculate_one(
          as.numeric(as.character(x_t) == level),
          as.numeric(as.character(x_c) == level),
          variable = var,
          level = level
        )
      })
      return(do.call(rbind, rows))
    }
    calculate_one(as.numeric(x_t), as.numeric(x_c), variable = var)
  })

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL

  n_balanced <- sum(result_df$balanced, na.rm = TRUE)
  message(sprintf("[assess_balance] %d/%d covariates balanced (|SMD| < %.2f)",
                  n_balanced, nrow(result_df), threshold))

  return(result_df)
}


#' Run Weighted Analysis
#'
#' @description
#' Fit regression models using IPTW weights with robust standard errors.
#'
#' @param data A data.frame or data.table containing all variables and weights.
#' @param exposure Character string specifying the exposure variable name.
#' @param outcome Character string specifying the outcome variable name (for logistic/linear).
#' @param covariates Character vector of covariate names. Default NULL.
#' @param weight_col Character string specifying the weight column name. Default "weight".
#' @param model_type Character string specifying model type: "cox", "logistic", or "linear".
#' @param endpoint Character vector of length 2 for Cox models: c("time", "status").
#' @param robust_se Logical; whether to use robust standard errors. Default TRUE.
#'
#' @return A data.frame with effect estimates and confidence intervals.
#'
#' @importFrom stats as.formula glm binomial lm coef confint vcov qnorm
#' @export
run_weighted_analysis <- function(data,
                                   exposure,
                                   outcome = NULL,
                                   covariates = NULL,
                                   weight_col = "weight",
                                   model_type = c("cox", "logistic", "linear"),
                                   endpoint = NULL,
                                   robust_se = TRUE) {

  model_type <- match.arg(model_type)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!weight_col %in% names(data)) {
    stop(sprintf("Weight column '%s' not found in data.", weight_col))
  }

  # Validate model-specific requirements
  if (model_type == "cox") {
    if (is.null(endpoint) || length(endpoint) != 2) {
      stop("For Cox models, 'endpoint' must be c('time', 'status').")
    }
    required_cols <- c(exposure, covariates, endpoint, weight_col)
  } else {
    if (is.null(outcome)) {
      stop("'outcome' is required for logistic and linear models.")
    }
    required_cols <- c(exposure, outcome, covariates, weight_col)
  }

  missing_vars <- setdiff(required_cols, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variables not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  message(sprintf("[run_weighted_analysis] Fitting weighted %s model...", model_type))

  # Build formula
  rhs <- exposure
  if (!is.null(covariates)) {
    rhs <- paste(c(exposure, covariates), collapse = " + ")
  }

  weights <- data[[weight_col]]

  if (model_type == "cox") {
    formula_str <- paste0("Surv(", endpoint[1], ", ", endpoint[2], ") ~ ", rhs)
    formula_obj <- stats::as.formula(formula_str)
    model <- coxph(formula_obj, data = data, weights = weights)

    # Extract results
    sum_model <- summary(model)
    coefs <- sum_model$coefficients
    conf <- sum_model$conf.int

    exp_row <- grep(paste0("^", exposure), rownames(coefs))[1]
    if (is.na(exp_row)) {
      stop("Exposure variable not found in model results.")
    }

    if (robust_se) {
      # Robust SE for Cox
      robust_vcov <- vcovHC(model, type = "HC0")
      robust_se_val <- sqrt(diag(robust_vcov))[exp_row]
      coef_val <- coefs[exp_row, "coef"]
      z_val <- coef_val / robust_se_val
      p_val <- 2 * stats::pnorm(-abs(z_val))
      hr <- exp(coef_val)
      lower95 <- exp(coef_val - 1.96 * robust_se_val)
      upper95 <- exp(coef_val + 1.96 * robust_se_val)
    } else {
      hr <- conf[exp_row, "exp(coef)"]
      lower95 <- conf[exp_row, "lower .95"]
      upper95 <- conf[exp_row, "upper .95"]
      p_val <- coefs[exp_row, "Pr(>|z|)"]
    }

    result <- data.frame(
      variable = exposure,
      HR = round(hr, 3),
      lower95 = round(lower95, 3),
      upper95 = round(upper95, 3),
      pvalue = signif(p_val, 3),
      stringsAsFactors = FALSE
    )

  } else if (model_type == "logistic") {
    formula_obj <- stats::as.formula(paste(outcome, "~", rhs))
    model <- stats::glm(formula_obj, data = data, family = stats::binomial(), weights = weights)

    exp_row <- grep(paste0("^", exposure), names(stats::coef(model)))[1]
    coef_val <- stats::coef(model)[exp_row]

    if (robust_se) {
      robust_vcov <- vcovHC(model, type = "HC0")
      coef_test <- coeftest(model, vcov. = robust_vcov)
      se_val <- coef_test[exp_row, "Std. Error"]
      p_val <- coef_test[exp_row, "Pr(>|z|)"]
    } else {
      se_val <- summary(model)$coefficients[exp_row, "Std. Error"]
      p_val <- summary(model)$coefficients[exp_row, "Pr(>|z|)"]
    }

    or <- exp(coef_val)
    lower95 <- exp(coef_val - 1.96 * se_val)
    upper95 <- exp(coef_val + 1.96 * se_val)

    result <- data.frame(
      variable = exposure,
      OR = round(or, 3),
      lower95 = round(lower95, 3),
      upper95 = round(upper95, 3),
      pvalue = signif(p_val, 3),
      stringsAsFactors = FALSE
    )

  } else {
    # Linear model
    formula_obj <- stats::as.formula(paste(outcome, "~", rhs))
    model <- stats::lm(formula_obj, data = data, weights = weights)

    exp_row <- grep(paste0("^", exposure), names(stats::coef(model)))[1]
    coef_val <- stats::coef(model)[exp_row]

    if (robust_se) {
      robust_vcov <- vcovHC(model, type = "HC0")
      coef_test <- coeftest(model, vcov. = robust_vcov)
      se_val <- coef_test[exp_row, "Std. Error"]
      p_val <- coef_test[exp_row, "Pr(>|t|)"]
    } else {
      se_val <- summary(model)$coefficients[exp_row, "Std. Error"]
      p_val <- summary(model)$coefficients[exp_row, "Pr(>|t|)"]
    }

    lower95 <- coef_val - 1.96 * se_val
    upper95 <- coef_val + 1.96 * se_val

    result <- data.frame(
      variable = exposure,
      beta = round(coef_val, 3),
      lower95 = round(lower95, 3),
      upper95 = round(upper95, 3),
      pvalue = signif(p_val, 3),
      stringsAsFactors = FALSE
    )
  }

  message("[run_weighted_analysis] Complete.")
  return(result)
}
