#' Run a Cox sensitivity-analysis suite
#'
#' @description
#' Fit a primary Cox model and common sensitivity models using the same
#' endpoint, exposure, and covariate structure. The suite currently supports
#' complete-case filtering, exclusion of early events, and additional covariate
#' adjustment sets.
#'
#' @param data A data.frame or data.table.
#' @param exposure Character vector of exposure variables.
#' @param covariates Optional character vector of primary adjustment covariates.
#' @param endpoint Character vector of length 2 giving survival time and status.
#' @param early_event_years Optional numeric vector of lag periods used to
#'   exclude events occurring at or before each cut point.
#' @param complete_case_covariates Optional covariates for a complete-case
#'   sensitivity dataset.
#' @param additional_covariate_sets Optional named list of extra covariate
#'   vectors. Each set is added to the primary covariates and refitted.
#' @param conf_level Confidence level for hazard-ratio intervals.
#' @param verbose Logical. If `TRUE`, print a compact summary.
#'
#' @return A list with class `ukb_sensitivity_suite` containing model objects,
#'   flow metadata, and a tidy summary table.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   time = rexp(100, 0.1),
#'   status = rbinom(100, 1, 0.3),
#'   exposure = rnorm(100),
#'   age = rnorm(100, 60, 5),
#'   sex = rbinom(100, 1, 0.5)
#' )
#' res <- ukb_sensitivity_suite(
#'   dat,
#'   exposure = "exposure",
#'   covariates = c("age", "sex"),
#'   endpoint = c("time", "status"),
#'   early_event_years = 1,
#'   verbose = FALSE
#' )
ukb_sensitivity_suite <- function(data,
                                  exposure,
                                  covariates = NULL,
                                  endpoint = c("outcome_surv_time", "outcome_status"),
                                  early_event_years = c(2, 4, 6),
                                  complete_case_covariates = NULL,
                                  additional_covariate_sets = NULL,
                                  conf_level = 0.95,
                                  verbose = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.character(exposure) || length(exposure) == 0 || anyNA(exposure)) {
    stop("`exposure` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.null(covariates) && (!is.character(covariates) || anyNA(covariates))) {
    stop("`covariates` must be NULL or a character vector.", call. = FALSE)
  }
  if (!is.character(endpoint) || length(endpoint) != 2 || anyNA(endpoint)) {
    stop("`endpoint` must be a character vector of length 2.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a number between 0 and 1.", call. = FALSE)
  }

  required <- unique(c(endpoint, exposure, covariates, complete_case_covariates, unlist(additional_covariate_sets)))
  missing_vars <- setdiff(required, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in `data`: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  status_values <- unique(data[[endpoint[[2]]]][!is.na(data[[endpoint[[2]]]])])
  if (!all(status_values %in% c(0, 1))) {
    stop("The status column must use 0/1 coding.", call. = FALSE)
  }

  models <- list()
  flows <- list()
  summaries <- list()

  fit_one <- function(label, fit_data, fit_covariates, sensitivity_type, detail = NA_character_) {
    fit <- .ukb_fit_sensitivity_cox(
      data = fit_data,
      exposure = exposure,
      covariates = fit_covariates,
      endpoint = endpoint,
      conf_level = conf_level
    )
    list(
      model = fit$model,
      summary = cbind(
        data.frame(
          analysis = label,
          sensitivity_type = sensitivity_type,
          detail = detail,
          stringsAsFactors = FALSE
        ),
        fit$summary
      ),
      label = label
    )
  }

  primary_fit <- fit_one("primary", data, covariates, "primary")
  models[[primary_fit$label]] <- primary_fit$model
  summaries[[primary_fit$label]] <- primary_fit$summary

  if (!is.null(complete_case_covariates)) {
    cc_data <- sensitivity_exclude_missing_covariates(
      data = data,
      covariates = complete_case_covariates,
      stepwise = TRUE,
      verbose = FALSE
    )
    flows[["complete_case"]] <- attr(cc_data, "complete_case_flow")
    cc_fit <- fit_one(
      "complete_case",
      cc_data,
      covariates,
      "complete_case",
      paste(complete_case_covariates, collapse = "+")
    )
    models[[cc_fit$label]] <- cc_fit$model
    summaries[[cc_fit$label]] <- cc_fit$summary
  }

  if (!is.null(early_event_years) && length(early_event_years) > 0) {
    for (yr in early_event_years) {
      lag_data <- sensitivity_exclude_early_events(
        data = data,
        endpoint = endpoint,
        n_years = yr,
        verbose = FALSE
      )
      label <- paste0("exclude_events_le_", yr, "y")
      flows[[label]] <- attr(lag_data, "sensitivity_info")
      lag_fit <- fit_one(label, lag_data, covariates, "early_event_exclusion", paste0(yr, " years"))
      models[[lag_fit$label]] <- lag_fit$model
      summaries[[lag_fit$label]] <- lag_fit$summary
    }
  }

  if (!is.null(additional_covariate_sets)) {
    if (!is.list(additional_covariate_sets)) {
      stop("`additional_covariate_sets` must be a named list.", call. = FALSE)
    }
    if (is.null(names(additional_covariate_sets)) || any(!nzchar(names(additional_covariate_sets)))) {
      names(additional_covariate_sets) <- paste0("additional_set_", seq_along(additional_covariate_sets))
    }
    for (nm in names(additional_covariate_sets)) {
      extra <- additional_covariate_sets[[nm]]
      if (!is.character(extra) || anyNA(extra)) {
        stop("Each additional covariate set must be a character vector.", call. = FALSE)
      }
      adj_fit <- fit_one(
        paste0("adjust_", nm),
        data,
        unique(c(covariates, extra)),
        "additional_adjustment",
        paste(extra, collapse = "+")
      )
      models[[adj_fit$label]] <- adj_fit$model
      summaries[[adj_fit$label]] <- adj_fit$summary
    }
  }

  summary_df <- do.call(rbind, summaries)
  rownames(summary_df) <- NULL
  out <- list(
    summary = summary_df,
    models = models,
    flows = flows,
    settings = list(
      exposure = exposure,
      covariates = covariates,
      endpoint = endpoint,
      early_event_years = early_event_years,
      complete_case_covariates = complete_case_covariates,
      additional_covariate_sets = additional_covariate_sets,
      conf_level = conf_level
    )
  )
  class(out) <- c("ukb_sensitivity_suite", class(out))

  if (isTRUE(verbose)) {
    print(out)
  }
  out
}

#' @export
print.ukb_sensitivity_suite <- function(x, ...) {
  cat("UKB sensitivity suite\n")
  cat("  analyses: ", length(unique(x$summary$analysis)), "\n", sep = "")
  cat("  exposure: ", paste(x$settings$exposure, collapse = ", "), "\n", sep = "")
  cat("  endpoint: ", paste(x$settings$endpoint, collapse = ", "), "\n", sep = "")
  print(x$summary[, c("analysis", "term", "HR", "lower95", "upper95", "pvalue", "n", "n_event")], row.names = FALSE)
  invisible(x)
}

.ukb_fit_sensitivity_cox <- function(data,
                                     exposure,
                                     covariates,
                                     endpoint,
                                     conf_level) {
  rhs <- paste(unique(c(exposure, covariates)), collapse = " + ")
  formula_obj <- as.formula(paste0("Surv(", endpoint[[1]], ", ", endpoint[[2]], ") ~ ", rhs))
  model <- coxph(formula_obj, data = data)
  sm <- summary(model, conf.int = conf_level)
  coefs <- sm$coefficients
  conf <- sm$conf.int
  lower_col <- grep("^lower", colnames(conf), value = TRUE)[1]
  upper_col <- grep("^upper", colnames(conf), value = TRUE)[1]

  data_terms <- data.frame(
    term = rownames(coefs),
    coef = unname(coefs[, "coef"]),
    se = unname(coefs[, "se(coef)"]),
    z = unname(coefs[, "z"]),
    HR = unname(conf[, "exp(coef)"]),
    lower95 = unname(conf[, lower_col]),
    upper95 = unname(conf[, upper_col]),
    pvalue = unname(coefs[, "Pr(>|z|)"]),
    n = unname(sm$n),
    n_event = unname(sm$nevent),
    stringsAsFactors = FALSE
  )
  list(model = model, summary = data_terms)
}
