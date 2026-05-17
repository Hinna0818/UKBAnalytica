#' Regression Extension Functions for UKBAnalytica
#'
#' @description
#' Additional regression helpers for common UK Biobank workflows, including
#' Cox proportional hazards diagnostics, grouped-exposure trend tests,
#' Fine-Gray competing-risk models, and lagged Cox sensitivity analyses.
#'
#' @name regression_extensions
#' @keywords internal
#' @importFrom survival Surv coxph finegray cox.zph
NULL

.coxext_require_survival <- function(caller) {
  invisible(TRUE)
}

.coxext_validate_logical_flag <- function(x, name) {
  if (!isTRUE(x) && !identical(x, FALSE)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
}

.coxext_validate_endpoint <- function(endpoint) {
  if (!is.character(endpoint) || length(endpoint) != 2 || anyNA(endpoint)) {
    stop("'endpoint' must be a character vector of length 2.", call. = FALSE)
  }
}

.coxext_make_rhs <- function(exposure, covariates = NULL) {
  rhs_terms <- c(exposure, covariates)
  paste(rhs_terms, collapse = " + ")
}

.coxext_fit_cox_model <- function(data,
                                  exposure,
                                  covariates = NULL,
                                  endpoint = c("time", "status"),
                                  ...) {
  .coxext_require_survival("runmulti_cox_zph")
  .coxext_validate_endpoint(endpoint)

  time_col <- endpoint[[1]]
  status_col <- endpoint[[2]]

  if (!is.numeric(data[[time_col]])) {
    stop(sprintf("Time column '%s' must be numeric.", time_col), call. = FALSE)
  }

  observed_status <- unique(data[[status_col]][!is.na(data[[status_col]])])
  if (!all(observed_status %in% c(0, 1))) {
    stop(sprintf(
      "Status column '%s' must use 0/1 coding for censored/event.",
      status_col
    ), call. = FALSE)
  }

  formula_obj <- stats::as.formula(
    paste0(
      "Surv(",
      time_col,
      ", ",
      status_col,
      ") ~ ",
      .coxext_make_rhs(exposure, covariates)
    )
  )

  needed_cols <- unique(c(time_col, status_col, exposure, covariates))
  model_data <- data[stats::complete.cases(data[, needed_cols, drop = FALSE]), needed_cols, drop = FALSE]

  if (nrow(model_data) == 0) {
    stop(sprintf("No complete-case rows available for exposure '%s'.", exposure), call. = FALSE)
  }

  fit_args <- c(
    list(formula = formula_obj, data = model_data),
    list(...)
  )
  if (is.null(fit_args$x)) {
    fit_args$x <- TRUE
  }
  if (is.null(fit_args$model)) {
    fit_args$model <- TRUE
  }
  fit <- do.call(coxph, fit_args)

  list(
    model = fit,
    formula = formula_obj,
    data = model_data,
    n = nrow(model_data),
    n_event = sum(model_data[[status_col]] == 1, na.rm = TRUE)
  )
}

.coxext_fit_standard_model <- function(model_type,
                                       data,
                                       exposure,
                                       outcome = NULL,
                                       covariates = NULL,
                                       endpoint = NULL,
                                       ...) {
  model_type <- match.arg(model_type, c("cox", "logistic", "linear"))

  if (model_type == "cox") {
    return(.coxext_fit_cox_model(
      data = data,
      exposure = exposure,
      covariates = covariates,
      endpoint = endpoint,
      ...
    ))
  }

  if (!is.character(outcome) || length(outcome) != 1 || is.na(outcome)) {
    stop("'outcome' must be a single character string.", call. = FALSE)
  }

  if (!outcome %in% names(data)) {
    stop(sprintf("Outcome variable '%s' not found in data.", outcome), call. = FALSE)
  }

  if (model_type == "logistic") {
    outcome_vals <- data[[outcome]]
    outcome_vals <- outcome_vals[!is.na(outcome_vals)]
    if (!all(sort(unique(outcome_vals)) %in% c(0, 1))) {
      stop(sprintf(
        "Outcome variable '%s' must be binary (0/1) for logistic regression.",
        outcome
      ), call. = FALSE)
    }
  }

  formula_obj <- stats::as.formula(
    paste(outcome, "~", .coxext_make_rhs(exposure, covariates))
  )

  needed_cols <- unique(c(outcome, exposure, covariates))
  model_data <- data[stats::complete.cases(data[, needed_cols, drop = FALSE]), needed_cols, drop = FALSE]

  if (nrow(model_data) == 0) {
    stop(sprintf("No complete-case rows available for exposure '%s'.", exposure), call. = FALSE)
  }

  fit <- switch(
    model_type,
    logistic = stats::glm(
      formula_obj,
      data = model_data,
      family = stats::binomial(),
      ...
    ),
    linear = stats::lm(
      formula_obj,
      data = model_data,
      ...
    )
  )

  list(
    model = fit,
    formula = formula_obj,
    data = model_data,
    n = nrow(model_data),
    n_event = if (model_type == "logistic") {
      sum(model_data[[outcome]] == 1, na.rm = TRUE)
    } else {
      NA_integer_
    }
  )
}

.coxext_get_term_coefficient_names <- function(model, term_name) {
  terms_obj <- stats::terms(model)
  mm <- if (!is.null(model$x)) {
    model$x
  } else {
    stats::model.matrix(model)
  }
  assign_vec <- attr(mm, "assign")
  term_labels <- attr(terms_obj, "term.labels")
  term_idx <- match(term_name, term_labels)

  if (is.na(term_idx)) {
    return(character(0))
  }

  colnames(mm)[assign_vec == term_idx]
}

.coxext_extract_effect_rows <- function(fit_info,
                                        exposure,
                                        model_type,
                                        null_value = NULL) {
  model_type <- match.arg(model_type, c("cox", "logistic", "linear", "finegray"))

  term_names <- .coxext_get_term_coefficient_names(fit_info$model, exposure)

  if (length(term_names) == 0) {
    stop(sprintf("No model term columns found for exposure '%s'.", exposure), call. = FALSE)
  }

  if (model_type %in% c("cox", "finegray")) {
    sum_model <- summary(fit_info$model)
    coef_mat <- sum_model$coefficients
    conf_mat <- sum_model$conf.int
    effect_col <- if (model_type == "cox") "HR" else "SHR"
    effect_values <- conf_mat[term_names, "exp(coef)"]
    lower_values <- conf_mat[term_names, "lower .95"]
    upper_values <- conf_mat[term_names, "upper .95"]
    p_values <- coef_mat[term_names, "Pr(>|z|)"]
  } else if (model_type == "logistic") {
    sum_model <- summary(fit_info$model)
    coef_mat <- sum_model$coefficients
    ci <- suppressWarnings(stats::confint(fit_info$model, parm = term_names, level = 0.95))
    if (!is.matrix(ci)) ci <- t(as.matrix(ci))
    effect_col <- "OR"
    effect_values <- exp(coef_mat[term_names, "Estimate"])
    lower_values <- exp(ci[, 1])
    upper_values <- exp(ci[, 2])
    p_values <- coef_mat[term_names, "Pr(>|z|)"]
  } else {
    sum_model <- summary(fit_info$model)
    coef_mat <- sum_model$coefficients
    ci <- stats::confint(fit_info$model, parm = term_names, level = 0.95)
    if (!is.matrix(ci)) ci <- t(as.matrix(ci))
    effect_col <- "beta"
    effect_values <- coef_mat[term_names, "Estimate"]
    lower_values <- ci[, 1]
    upper_values <- ci[, 2]
    p_values <- coef_mat[term_names, "Pr(>|t|)"]
  }

  out <- data.frame(
    variable = exposure,
    level = ifelse(term_names == exposure, NA_character_, sub(paste0("^", exposure), "", term_names)),
    n = fit_info$n,
    n_event = fit_info$n_event,
    lower95 = round(lower_values, 3),
    upper95 = round(upper_values, 3),
    pvalue = signif(p_values, 3),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out[[effect_col]] <- round(effect_values, 3)
  out <- out[, c("variable", "level", "n", "n_event", effect_col, "lower95", "upper95", "pvalue")]

  if (!is.null(null_value)) {
    ref_level <- levels(fit_info$data[[exposure]])[[1]]
    ref_row <- data.frame(
      variable = exposure,
      level = ref_level,
      n = fit_info$n,
      n_event = fit_info$n_event,
      lower95 = NA_real_,
      upper95 = NA_real_,
      pvalue = NA_real_,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    ref_row[[effect_col]] <- null_value
    ref_row <- ref_row[, c("variable", "level", "n", "n_event", effect_col, "lower95", "upper95", "pvalue")]
    out <- rbind(ref_row, out)
    rownames(out) <- NULL
  }

  out
}

.coxext_make_trend_score <- function(x,
                                     score_method = c("integer", "median", "custom"),
                                     custom_scores = NULL,
                                     variable = NULL) {
  score_method <- match.arg(score_method)

  if (isTRUE(all(is.na(x)))) {
    return(rep(NA_real_, length(x)))
  }

  x_factor <- if (is.factor(x)) {
    x
  } else {
    factor(x, ordered = is.ordered(x))
  }

  lvls <- levels(x_factor)
  if (length(lvls) < 2) {
    stop(sprintf("Exposure '%s' must contain at least 2 levels for trend analysis.", variable), call. = FALSE)
  }

  if (score_method == "median") {
    stop(
      "'score_method = \"median\"' is reserved for a later version. Please use 'integer' or 'custom'.",
      call. = FALSE
    )
  }

  if (score_method == "integer") {
    score_map <- stats::setNames(seq_along(lvls), lvls)
  } else {
    if (is.null(custom_scores) || is.null(variable) || !variable %in% names(custom_scores)) {
      stop(sprintf(
        "Custom scores for exposure '%s' were not provided.",
        variable
      ), call. = FALSE)
    }
    score_map <- custom_scores[[variable]]
    if (is.null(names(score_map))) {
      stop(sprintf(
        "Custom scores for exposure '%s' must be a named vector keyed by exposure levels.",
        variable
      ), call. = FALSE)
    }
    if (!all(lvls %in% names(score_map))) {
      missing_levels <- setdiff(lvls, names(score_map))
      stop(sprintf(
        "Custom scores for exposure '%s' are missing levels: %s",
        variable,
        paste(missing_levels, collapse = ", ")
      ), call. = FALSE)
    }
    score_map <- score_map[lvls]
  }

  as.numeric(score_map[as.character(x_factor)])
}

.coxext_fit_finegray_model <- function(data,
                                       exposure,
                                       covariates = NULL,
                                       time_col,
                                       status_col,
                                       ...) {
  .coxext_require_survival("runmulti_competing")

  formula_fg <- stats::as.formula(
    paste0(
      "Surv(",
      time_col,
      ", ",
      status_col,
      ") ~ ",
      .coxext_make_rhs(exposure, covariates)
    )
  )

  fg_data <- finegray(
    formula_fg,
    data = data,
    etype = "event"
  )

  fit_formula <- stats::as.formula(
    paste0(
      "Surv(fgstart, fgstop, fgstatus) ~ ",
      .coxext_make_rhs(exposure, covariates)
    )
  )

  fit_args <- c(
    list(
      formula = fit_formula,
      data = fg_data,
      weights = fg_data$fgwt
    ),
    list(...)
  )
  if (is.null(fit_args$x)) {
    fit_args$x <- TRUE
  }
  if (is.null(fit_args$model)) {
    fit_args$model <- TRUE
  }
  fit <- do.call(coxph, fit_args)

  list(
    model = fit,
    formula = fit_formula,
    data = fg_data,
    n = attr(data, "coxext_n_input"),
    n_event = attr(data, "coxext_n_event"),
    n_compete = attr(data, "coxext_n_compete")
  )
}

.coxext_normalize_competing_status <- function(data,
                                               event_col,
                                               compete_col = NULL,
                                               event_value = 1,
                                               compete_value = 2) {
  if (!event_col %in% names(data)) {
    stop(sprintf("Event column '%s' not found in data.", event_col), call. = FALSE)
  }

  if (is.null(compete_col)) {
    raw <- data[[event_col]]
    status <- ifelse(
      is.na(raw),
      NA_character_,
      ifelse(raw == event_value, "event",
        ifelse(raw == compete_value, "compete", "censor")
      )
    )
  } else {
    if (!compete_col %in% names(data)) {
      stop(sprintf("Competing-event column '%s' not found in data.", compete_col), call. = FALSE)
    }
    event_raw <- data[[event_col]]
    compete_raw <- data[[compete_col]]

    event_obs <- unique(event_raw[!is.na(event_raw)])
    compete_obs <- unique(compete_raw[!is.na(compete_raw)])
    if (!all(event_obs %in% c(0, 1)) || !all(compete_obs %in% c(0, 1))) {
      stop("Binary event and competing-event columns must use 0/1 coding.", call. = FALSE)
    }

    status <- ifelse(
      is.na(event_raw) | is.na(compete_raw),
      NA_character_,
      ifelse(event_raw == 1, "event",
        ifelse(compete_raw == 1, "compete", "censor")
      )
    )
  }

  factor(status, levels = c("censor", "event", "compete"))
}

.coxext_bind_results <- function(results, caller) {
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) {
    stop(sprintf("No valid results were produced by %s().", caller), call. = FALSE)
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

.coxext_resolve_sensitivity_info <- function(info,
                                             n_input_fallback,
                                             n_output_fallback) {
  if (is.list(info) &&
      !is.null(names(info)) &&
      all(c("n_input", "n_removed") %in% names(info))) {
    resolved <- info
  } else if (is.list(info) && length(info) > 0 && is.list(info[[length(info)]])) {
    resolved <- info[[length(info)]]
  } else {
    resolved <- list()
  }

  n_input <- resolved[["n_input"]]
  if (is.null(n_input) || length(n_input) != 1 || is.na(n_input)) {
    n_input <- n_input_fallback
  }

  n_removed <- resolved[["n_removed"]]
  if (is.null(n_removed) || length(n_removed) != 1 || is.na(n_removed)) {
    n_removed <- n_input - n_output_fallback
  }

  n_output <- resolved[["n_output"]]
  if (is.null(n_output) || length(n_output) != 1 || is.na(n_output)) {
    n_output <- n_output_fallback
  }

  list(
    n_input = as.integer(n_input),
    n_removed = as.integer(n_removed),
    n_output = as.integer(n_output)
  )
}

#' Diagnose Proportional Hazards Assumptions for a Cox Model
#'
#' @param model A fitted \code{coxph()} model.
#' @param transform Character scalar passed to \code{cox.zph()}.
#' @param terms Logical; keep term-level rows.
#' @param global Logical; keep the GLOBAL row.
#' @param alpha Numeric threshold for flagging PH violations.
#' @param return_object Logical; if TRUE, include the raw \code{cox.zph} object.
#'
#' @return A list containing a tidy diagnostics table, the global p-value, and
#'   optionally the raw \code{cox.zph} object.
#' @export
ukb_cox_diagnostics <- function(model,
                                transform = c("km", "rank", "identity"),
                                terms = TRUE,
                                global = TRUE,
                                alpha = 0.05,
                                return_object = TRUE) {
  .coxext_require_survival("ukb_cox_diagnostics")
  .coxext_validate_logical_flag(terms, "terms")
  .coxext_validate_logical_flag(global, "global")
  .coxext_validate_logical_flag(return_object, "return_object")

  transform <- match.arg(transform)

  if (missing(model) || is.null(model) || !"coxph" %in% class(model)) {
    stop("'model' must be a fitted 'coxph' object.", call. = FALSE)
  }

  zph <- cox.zph(model, transform = transform, terms = isTRUE(terms), global = isTRUE(global))
  zph_table <- as.data.frame(zph$table, stringsAsFactors = FALSE)
  zph_table$term <- rownames(zph$table)
  rownames(zph_table) <- NULL

  if ("p" %in% names(zph_table)) {
    names(zph_table)[names(zph_table) == "p"] <- "pvalue"
  }
  if ("chisq" %in% names(zph_table)) {
    zph_table$chisq <- round(zph_table$chisq, 3)
  }
  if ("pvalue" %in% names(zph_table)) {
    zph_table$pvalue <- signif(zph_table$pvalue, 3)
    zph_table$ph_violation <- zph_table$pvalue < alpha
  } else {
    zph_table$ph_violation <- NA
  }

  global_pvalue <- NA_real_
  if ("GLOBAL" %in% zph_table$term && "pvalue" %in% names(zph_table)) {
    global_pvalue <- zph_table$pvalue[zph_table$term == "GLOBAL"][[1]]
  }

  out <- list(
    table = zph_table,
    global_pvalue = global_pvalue
  )
  if (isTRUE(return_object)) {
    out$cox_zph <- zph
  }

  out
}

#' Run Multiple Cox Models with PH Diagnostics
#'
#' @param data A data.frame or data.table.
#' @param main_var Character vector of exposure variable names.
#' @param covariates Optional character vector of covariate names.
#' @param endpoint Character vector of length 2: \code{c(time, status)}.
#' @param transform Character scalar passed to \code{cox.zph()}.
#' @param alpha Numeric threshold for flagging PH violations.
#' @param keep_models Logical; if TRUE, attach fitted models as an attribute.
#' @param ... Additional arguments passed to \code{coxph()}.
#'
#' @return A data.frame with effect estimates and PH-diagnostic columns.
#' @export
runmulti_cox_zph <- function(data,
                             main_var,
                             covariates = NULL,
                             endpoint = c("time", "status"),
                             transform = c("km", "rank", "identity"),
                             alpha = 0.05,
                             keep_models = FALSE,
                             ...) {
  .coxext_require_survival("runmulti_cox_zph")
  .coxext_validate_endpoint(endpoint)
  .coxext_validate_logical_flag(keep_models, "keep_models")
  .validate_regression_inputs(data, main_var, covariates, endpoint)

  transform <- match.arg(transform)
  results <- vector("list", length(main_var))
  model_store <- list()

  for (i in seq_along(main_var)) {
    var <- main_var[[i]]

    fit_info <- tryCatch(
      .coxext_fit_cox_model(
        data = data,
        exposure = var,
        covariates = covariates,
        endpoint = endpoint,
        ...
      ),
      error = function(e) {
        warning(sprintf("[runmulti_cox_zph] Skipping '%s': %s", var, e$message), call. = FALSE)
        NULL
      }
    )

    if (is.null(fit_info)) {
      next
    }

    effect_df <- tryCatch(
      .coxext_extract_effect_rows(
        fit_info = fit_info,
        exposure = var,
        model_type = "cox"
      ),
      error = function(e) {
        warning(sprintf("[runmulti_cox_zph] Skipping '%s': %s", var, e$message), call. = FALSE)
        NULL
      }
    )
    if (is.null(effect_df)) {
      next
    }

    diag_info <- ukb_cox_diagnostics(
      model = fit_info$model,
      transform = transform,
      terms = TRUE,
      global = TRUE,
      alpha = alpha,
      return_object = FALSE
    )

    diag_table <- diag_info$table
    exposure_pvalue <- if (var %in% diag_table$term) {
      diag_table$pvalue[diag_table$term == var][[1]]
    } else {
      NA_real_
    }

    effect_df$zph_pvalue <- exposure_pvalue
    effect_df$zph_global_pvalue <- diag_info$global_pvalue
    effect_df$ph_violation <- !is.na(effect_df$zph_pvalue) & effect_df$zph_pvalue < alpha
    effect_df$ph_global_violation <- !is.na(effect_df$zph_global_pvalue) & effect_df$zph_global_pvalue < alpha

    results[[i]] <- effect_df
    model_store[[var]] <- fit_info$model
  }

  out <- .coxext_bind_results(results, "runmulti_cox_zph")
  if (isTRUE(keep_models)) {
    attr(out, "cox_models") <- model_store
  }
  out
}

#' Run Grouped-Exposure Trend Tests
#'
#' @param data A data.frame or data.table.
#' @param main_var Character vector of grouped exposure variable names.
#' @param outcome Outcome column for logistic or linear models.
#' @param covariates Optional character vector of covariates.
#' @param model_type One of \code{"cox"}, \code{"logistic"}, or \code{"linear"}.
#' @param endpoint Character vector of length 2 for Cox models.
#' @param ref_level Optional reference level applied to every grouped exposure.
#' @param score_method One of \code{"integer"}, \code{"median"}, or \code{"custom"}.
#' @param custom_scores Optional named list of custom score mappings.
#' @param include_level_estimates Logical; if TRUE, include category-specific estimates.
#' @param ... Additional arguments passed to the fitted model.
#'
#' @return A data.frame containing grouped-effect estimates and a repeated
#'   \code{p_trend} column for each exposure.
#' @export
runmulti_trend <- function(data,
                           main_var,
                           outcome = NULL,
                           covariates = NULL,
                           model_type = c("cox", "logistic", "linear"),
                           endpoint = NULL,
                           ref_level = NULL,
                           score_method = c("integer", "median", "custom"),
                           custom_scores = NULL,
                           include_level_estimates = TRUE,
                           ...) {
  model_type <- match.arg(model_type)
  score_method <- match.arg(score_method)
  .coxext_validate_logical_flag(include_level_estimates, "include_level_estimates")

  required_cols <- if (model_type == "cox") {
    .coxext_require_survival("runmulti_trend")
    .coxext_validate_endpoint(endpoint)
    endpoint
  } else {
    outcome
  }
  .validate_regression_inputs(data, main_var, covariates, required_cols)

  results <- vector("list", length(main_var))

  for (i in seq_along(main_var)) {
    var <- main_var[[i]]
    work_dt <- if (data.table::is.data.table(data)) {
      data.table::copy(data)
    } else {
      data
    }

    exposure_vec <- work_dt[[var]]
    if (is.null(exposure_vec)) {
      warning(sprintf("[runmulti_trend] Skipping '%s': variable not found.", var), call. = FALSE)
      next
    }

    exposure_factor <- if (is.factor(exposure_vec)) {
      exposure_vec
    } else {
      factor(exposure_vec, ordered = is.ordered(exposure_vec))
    }

    if (length(levels(exposure_factor)) < 2) {
      warning(sprintf("[runmulti_trend] Skipping '%s': fewer than 2 non-missing levels.", var), call. = FALSE)
      next
    }

    if (!is.null(ref_level)) {
      if (!ref_level %in% levels(exposure_factor)) {
        warning(sprintf(
          "[runmulti_trend] Skipping '%s': ref_level '%s' not found in exposure levels.",
          var, ref_level
        ), call. = FALSE)
        next
      }
      exposure_factor <- stats::relevel(exposure_factor, ref = ref_level)
    }

    work_dt$.coxext_group_factor <- exposure_factor
    work_dt$.coxext_group_score <- tryCatch(
      .coxext_make_trend_score(
        x = exposure_factor,
        score_method = score_method,
        custom_scores = custom_scores,
        variable = var
      ),
      error = function(e) {
        warning(sprintf("[runmulti_trend] Skipping '%s': %s", var, e$message), call. = FALSE)
        NULL
      }
    )

    if (is.null(work_dt$.coxext_group_score)) {
      next
    }

    factor_fit <- tryCatch(
      .coxext_fit_standard_model(
        model_type = model_type,
        data = work_dt,
        exposure = ".coxext_group_factor",
        outcome = outcome,
        covariates = covariates,
        endpoint = endpoint,
        ...
      ),
      error = function(e) {
        warning(sprintf("[runmulti_trend] Skipping '%s': %s", var, e$message), call. = FALSE)
        NULL
      }
    )
    if (is.null(factor_fit)) {
      next
    }

    trend_fit <- tryCatch(
      .coxext_fit_standard_model(
        model_type = model_type,
        data = work_dt,
        exposure = ".coxext_group_score",
        outcome = outcome,
        covariates = covariates,
        endpoint = endpoint,
        ...
      ),
      error = function(e) {
        warning(sprintf("[runmulti_trend] Skipping '%s': %s", var, e$message), call. = FALSE)
        NULL
      }
    )
    if (is.null(trend_fit)) {
      next
    }

    level_df <- if (isTRUE(include_level_estimates)) {
      .coxext_extract_effect_rows(
        fit_info = factor_fit,
        exposure = ".coxext_group_factor",
        model_type = model_type,
        null_value = if (model_type == "linear") 0 else 1
      )
    } else {
      data.frame(
        variable = var,
        level = NA_character_,
        n = factor_fit$n,
        n_event = factor_fit$n_event,
        estimate = NA_real_,
        lower95 = NA_real_,
        upper95 = NA_real_,
        pvalue = NA_real_,
        stringsAsFactors = FALSE
      )
    }

    if (include_level_estimates) {
      level_df$variable <- var
    }

    trend_df <- .coxext_extract_effect_rows(
      fit_info = trend_fit,
      exposure = ".coxext_group_score",
      model_type = model_type
    )

    if (nrow(trend_df) != 1) {
      warning(sprintf(
        "[runmulti_trend] Skipping '%s': trend model returned an unexpected number of rows.",
        var
      ), call. = FALSE)
      next
    }

    trend_effect_col <- intersect(c("HR", "OR", "beta"), names(trend_df))
    if (length(trend_effect_col) != 1) {
      warning(sprintf(
        "[runmulti_trend] Skipping '%s': could not determine trend effect column.",
        var
      ), call. = FALSE)
      next
    }

    trend_estimate <- trend_df[[trend_effect_col]]
    level_df$p_trend <- trend_df$pvalue[[1]]
    level_df$trend_estimate <- trend_estimate[[1]]
    level_df$trend_lower95 <- trend_df$lower95[[1]]
    level_df$trend_upper95 <- trend_df$upper95[[1]]
    level_df$score_method <- score_method

    if ("HR" %in% names(level_df)) {
      level_df$estimate <- level_df$HR
      level_df$HR <- NULL
    } else if ("OR" %in% names(level_df)) {
      level_df$estimate <- level_df$OR
      level_df$OR <- NULL
    } else if ("beta" %in% names(level_df)) {
      level_df$estimate <- level_df$beta
      level_df$beta <- NULL
    }

    level_df <- level_df[, c(
      "variable", "level", "n", "n_event", "estimate",
      "lower95", "upper95", "pvalue", "p_trend",
      "trend_estimate", "trend_lower95", "trend_upper95", "score_method"
    )]

    results[[i]] <- level_df
  }

  .coxext_bind_results(results, "runmulti_trend")
}

#' Run Multiple Fine-Gray Competing-Risk Models
#'
#' @param data A data.frame or data.table.
#' @param main_var Character vector of exposure variable names.
#' @param covariates Optional character vector of covariates.
#' @param time_col Follow-up time column.
#' @param event_col Event-status column, or the primary-event column in dual-column mode.
#' @param compete_col Optional competing-event column in dual-column mode.
#' @param event_value Event code used in single-column mode.
#' @param compete_value Competing-event code used in single-column mode.
#' @param conf_level Confidence level, reserved for future use.
#' @param ... Additional arguments passed to the weighted Cox fit.
#'
#' @return A data.frame with subdistribution hazard ratios.
#' @export
runmulti_competing <- function(data,
                               main_var,
                               covariates = NULL,
                               time_col,
                               event_col,
                               compete_col = NULL,
                               event_value = 1,
                               compete_value = 2,
                               conf_level = 0.95,
                               ...) {
  .coxext_require_survival("runmulti_competing")

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.character(main_var) || length(main_var) == 0) {
    stop("'main_var' must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(time_col) || length(time_col) != 1 || is.na(time_col)) {
    stop("'time_col' must be a single character string.", call. = FALSE)
  }
  if (!is.character(event_col) || length(event_col) != 1 || is.na(event_col)) {
    stop("'event_col' must be a single character string.", call. = FALSE)
  }
  if (!is.null(compete_col) && (!is.character(compete_col) || length(compete_col) != 1 || is.na(compete_col))) {
    stop("'compete_col' must be NULL or a single character string.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number between 0 and 1.", call. = FALSE)
  }

  needed_cols <- unique(c(main_var, covariates, time_col, event_col, compete_col))
  missing_cols <- setdiff(needed_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "The following variables are not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.numeric(data[[time_col]])) {
    stop(sprintf("Time column '%s' must be numeric.", time_col), call. = FALSE)
  }

  status_factor <- .coxext_normalize_competing_status(
    data = data,
    event_col = event_col,
    compete_col = compete_col,
    event_value = event_value,
    compete_value = compete_value
  )

  results <- vector("list", length(main_var))

  for (i in seq_along(main_var)) {
    var <- main_var[[i]]
    subset_cols <- unique(c(var, covariates, time_col))
    work_dt <- if (data.table::is.data.table(data)) {
      data.table::copy(data)
    } else {
      data
    }
    work_dt$.coxext_fg_status <- status_factor
    model_dt <- work_dt[stats::complete.cases(work_dt[, c(subset_cols, ".coxext_fg_status"), drop = FALSE]), , drop = FALSE]

    if (nrow(model_dt) == 0) {
      warning(sprintf("[runmulti_competing] Skipping '%s': no complete-case rows available.", var), call. = FALSE)
      next
    }

    model_dt <- model_dt[!is.na(model_dt[[time_col]]) & model_dt[[time_col]] >= 0, , drop = FALSE]
    if (nrow(model_dt) == 0) {
      warning(sprintf("[runmulti_competing] Skipping '%s': no rows with valid non-negative follow-up time.", var), call. = FALSE)
      next
    }

    n_event <- sum(model_dt$.coxext_fg_status == "event", na.rm = TRUE)
    n_compete <- sum(model_dt$.coxext_fg_status == "compete", na.rm = TRUE)

    if (n_event == 0) {
      warning(sprintf("[runmulti_competing] Skipping '%s': no primary events in analysis set.", var), call. = FALSE)
      next
    }

    attr(model_dt, "coxext_n_input") <- nrow(model_dt)
    attr(model_dt, "coxext_n_event") <- n_event
    attr(model_dt, "coxext_n_compete") <- n_compete

    fit_info <- tryCatch(
      .coxext_fit_finegray_model(
        data = model_dt,
        exposure = var,
        covariates = covariates,
        time_col = time_col,
        status_col = ".coxext_fg_status",
        ...
      ),
      error = function(e) {
        warning(sprintf("[runmulti_competing] Skipping '%s': %s", var, e$message), call. = FALSE)
        NULL
      }
    )
    if (is.null(fit_info)) {
      next
    }

    effect_df <- tryCatch(
      .coxext_extract_effect_rows(
        fit_info = fit_info,
        exposure = var,
        model_type = "finegray"
      ),
      error = function(e) {
        warning(sprintf("[runmulti_competing] Skipping '%s': %s", var, e$message), call. = FALSE)
        NULL
      }
    )
    if (is.null(effect_df)) {
      next
    }

    effect_df$n_compete <- fit_info$n_compete
    effect_df <- effect_df[, c("variable", "level", "n", "n_event", "n_compete", "SHR", "lower95", "upper95", "pvalue")]
    results[[i]] <- effect_df
  }

  .coxext_bind_results(results, "runmulti_competing")
}

#' Run Lagged Cox Sensitivity Analyses
#'
#' @param data A data.frame or data.table.
#' @param main_var Character vector of exposure variable names.
#' @param covariates Optional character vector of covariates.
#' @param endpoint Character vector of length 2: \code{c(time, status)}.
#' @param lag_years Numeric vector of lag windows in years. \code{0} means no filtering.
#' @param verbose Logical; print progress messages.
#' @param ... Additional arguments passed to \code{coxph()}.
#'
#' @return A data.frame containing lag-specific hazard-ratio estimates.
#' @export
runmulti_cox_lag <- function(data,
                             main_var,
                             covariates = NULL,
                             endpoint = c("time", "status"),
                             lag_years = c(0, 1, 2, 5),
                             verbose = TRUE,
                             ...) {
  .coxext_require_survival("runmulti_cox_lag")
  .coxext_validate_endpoint(endpoint)
  .coxext_validate_logical_flag(verbose, "verbose")
  .validate_regression_inputs(data, main_var, covariates, endpoint)

  if (!is.numeric(lag_years) || length(lag_years) == 0 || anyNA(lag_years) ||
      any(!is.finite(lag_years)) || any(lag_years < 0)) {
    stop("'lag_years' must be a non-empty numeric vector of non-negative values.", call. = FALSE)
  }

  lag_years <- sort(unique(lag_years))
  results <- vector("list", length(lag_years))

  for (i in seq_along(lag_years)) {
    lag_n <- lag_years[[i]]

    if (isTRUE(verbose)) {
      message(sprintf("[runmulti_cox_lag] Running lag = %s years", lag_n))
    }

    if (lag_n == 0) {
      lag_dt <- data
      info <- list(
        n_input = nrow(data),
        n_removed = 0L,
        n_output = nrow(data)
      )
    } else {
      lag_dt <- sensitivity_exclude_early_events(
        data = data,
        endpoint = endpoint,
        n_years = lag_n,
        copy = TRUE,
        verbose = FALSE
      )
      info <- .coxext_resolve_sensitivity_info(
        info = attr(lag_dt, "sensitivity_info", exact = TRUE),
        n_input_fallback = nrow(data),
        n_output_fallback = nrow(lag_dt)
      )
    }

    lag_results <- vector("list", length(main_var))
    for (j in seq_along(main_var)) {
      var <- main_var[[j]]
      fit_info <- tryCatch(
        .coxext_fit_cox_model(
          data = lag_dt,
          exposure = var,
          covariates = covariates,
          endpoint = endpoint,
          ...
        ),
        error = function(e) {
          warning(sprintf(
            "[runmulti_cox_lag] Skipping '%s' at lag %s: %s",
            var, lag_n, e$message
          ), call. = FALSE)
          NULL
        }
      )
      if (is.null(fit_info)) {
        next
      }

      effect_df <- tryCatch(
        .coxext_extract_effect_rows(
          fit_info = fit_info,
          exposure = var,
          model_type = "cox"
        ),
        error = function(e) {
          warning(sprintf(
            "[runmulti_cox_lag] Skipping '%s' at lag %s: %s",
            var, lag_n, e$message
          ), call. = FALSE)
          NULL
        }
      )
      if (is.null(effect_df)) {
        next
      }

      effect_df$lag_years <- lag_n
      effect_df$n_input <- info$n_input
      effect_df$n_removed <- info$n_removed
      effect_df <- effect_df[, c("lag_years", "variable", "level", "n_input", "n_removed", "n", "n_event", "HR", "lower95", "upper95", "pvalue")]
      lag_results[[j]] <- effect_df
    }

    lag_results <- Filter(Negate(is.null), lag_results)
    if (length(lag_results) == 0) {
      warning(sprintf(
        "[runmulti_cox_lag] No valid results were produced for lag %s years.",
        lag_n
      ), call. = FALSE)
      results[[i]] <- NULL
    } else {
      results[[i]] <- do.call(rbind, lag_results)
      rownames(results[[i]]) <- NULL
    }
  }

  .coxext_bind_results(results, "runmulti_cox_lag")
}
