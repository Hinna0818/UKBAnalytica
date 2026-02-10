#' Run multiple Cox proportional hazards models
#'
#' @description
#' Fit Cox proportional hazards models for each main variable separately.
#' When \code{covariates} is \code{NULL}, univariate models are fitted.
#' Otherwise, multivariate models adjusting for the specified covariates are fitted.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param main_var A character vector of main variable names to test.
#' @param covariates A character vector of covariate names to adjust for. Default \code{NULL} (univariate).
#' @param endpoint A character vector of length 2: \code{c("time", "status")}, indicating survival time and event columns.
#' @param ... Additional arguments passed to \code{survival::coxph()}.
#'
#' @importFrom stats as.formula confint
#' @return A data.frame with columns: \code{variable}, \code{HR}, \code{lower95}, \code{upper95}, \code{pvalue}.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' # Univariate Cox
#' res <- runmulti_cox(lung, main_var = c("age", "sex"), endpoint = c("time", "status"))
#'
#' # Multivariate Cox
#' res <- runmulti_cox(lung, main_var = c("age", "sex"),
#'                     covariates = c("ph.ecog"), endpoint = c("time", "status"))
#' }
#'
#' @export
runmulti_cox <- function(data,
                         main_var,
                         covariates = NULL,
                         endpoint = c("time", "status"),
                         ...) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for runmulti_cox(). Please install it with: install.packages('survival')")
  }

  stopifnot(length(endpoint) == 2)
  .validate_regression_inputs(data, main_var, covariates, c(endpoint))

  results <- list()

  for (var in main_var) {
    # construct formula
    rhs <- var
    if (!is.null(covariates)) {
      rhs <- paste(c(var, covariates), collapse = " + ")
    }
    formula_str <- paste0("survival::Surv(", endpoint[1], ", ", endpoint[2], ") ~ ", rhs)
    formula_obj <- stats::as.formula(formula_str)

    # fit model
    model <- survival::coxph(formula_obj, data = data, ...)

    # extract results
    sum_model <- summary(model)
    coefs <- sum_model$coefficients
    conf <- sum_model$conf.int
    main_row <- coefs[rownames(coefs) == var, , drop = FALSE]
    conf_row <- conf[rownames(conf) == var, , drop = FALSE]

    results[[var]] <- data.frame(
      variable = var,
      HR = round(conf_row[, "exp(coef)"], 3),
      lower95 = round(conf_row[, "lower .95"], 3),
      upper95 = round(conf_row[, "upper .95"], 3),
      pvalue = signif(main_row[, "Pr(>|z|)"], 3),
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Run multiple linear regression models
#'
#' @description
#' Fit linear regression models (\code{lm}) for each main variable separately.
#' When \code{covariates} is \code{NULL}, univariate models are fitted.
#' Otherwise, multivariate models adjusting for the specified covariates are fitted.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param main_var A character vector of main variable names to test.
#' @param covariates A character vector of covariate names to adjust for. Default \code{NULL} (univariate).
#' @param outcome A character string specifying the outcome (dependent) variable name.
#' @param ... Additional arguments passed to \code{stats::lm()}.
#'
#' @importFrom stats as.formula lm confint coef
#' @return A data.frame with columns: \code{variable}, \code{beta}, \code{lower95}, \code{upper95}, \code{pvalue}.
#'
#' @examples
#' \dontrun{
#' # Univariate linear regression
#' res <- runmulti_lm(mtcars, main_var = c("hp", "wt"), outcome = "mpg")
#'
#' # Multivariate linear regression
#' res <- runmulti_lm(mtcars, main_var = c("hp", "wt"),
#'                    covariates = c("cyl"), outcome = "mpg")
#' }
#'
#' @export
runmulti_lm <- function(data,
                        main_var,
                        covariates = NULL,
                        outcome,
                        ...) {

  .validate_regression_inputs(data, main_var, covariates, outcome)

  results <- list()

  for (var in main_var) {
    # construct formula
    rhs <- var
    if (!is.null(covariates)) {
      rhs <- paste(c(var, covariates), collapse = " + ")
    }
    formula_obj <- stats::as.formula(paste(outcome, "~", rhs))

    # fit model
    model <- stats::lm(formula_obj, data = data, ...)

    # extract results
    sum_model <- summary(model)
    coefs <- sum_model$coefficients
    ci <- stats::confint(model, parm = var, level = 0.95)
    if (!is.matrix(ci)) ci <- t(as.matrix(ci))
    main_row <- coefs[rownames(coefs) == var, , drop = FALSE]

    results[[var]] <- data.frame(
      variable = var,
      beta = round(main_row[, "Estimate"], 3),
      lower95 = round(ci[, 1], 3),
      upper95 = round(ci[, 2], 3),
      pvalue = signif(main_row[, "Pr(>|t|)"], 3),
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Run multiple logistic regression models
#'
#' @description
#' Fit logistic regression models (\code{glm} with \code{family = binomial}) for each main variable separately.
#' When \code{covariates} is \code{NULL}, univariate models are fitted.
#' Otherwise, multivariate models adjusting for the specified covariates are fitted.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param main_var A character vector of main variable names to test.
#' @param covariates A character vector of covariate names to adjust for. Default \code{NULL} (univariate).
#' @param outcome A character string specifying the binary outcome (dependent) variable name (0/1).
#' @param ... Additional arguments passed to \code{stats::glm()}.
#'
#' @importFrom stats as.formula glm binomial confint coef
#' @return A data.frame with columns: \code{variable}, \code{OR}, \code{lower95}, \code{upper95}, \code{pvalue}.
#'
#' @examples
#' \dontrun{
#' # Create binary outcome
#' mtcars$am_bin <- ifelse(mtcars$am == 1, 1, 0)
#'
#' # Univariate logistic regression
#' res <- runmulti_logit(mtcars, main_var = c("hp", "wt"), outcome = "am_bin")
#'
#' # Multivariate logistic regression
#' res <- runmulti_logit(mtcars, main_var = c("hp", "wt"),
#'                       covariates = c("cyl"), outcome = "am_bin")
#' }
#'
#' @export
runmulti_logit <- function(data,
                           main_var,
                           covariates = NULL,
                           outcome,
                           ...) {

  .validate_regression_inputs(data, main_var, covariates, outcome)

  # validate binary outcome
  outcome_vals <- data[[outcome]]
  outcome_vals <- outcome_vals[!is.na(outcome_vals)]
  unique_vals <- sort(unique(outcome_vals))
  if (!all(unique_vals %in% c(0, 1))) {
    stop(sprintf("Outcome variable '%s' must be binary (0/1). Found values: %s",
                 outcome, paste(unique_vals, collapse = ", ")))
  }

  results <- list()

  for (var in main_var) {
    # construct formula
    rhs <- var
    if (!is.null(covariates)) {
      rhs <- paste(c(var, covariates), collapse = " + ")
    }
    formula_obj <- stats::as.formula(paste(outcome, "~", rhs))

    # fit model
    model <- stats::glm(formula_obj, data = data, family = stats::binomial(), ...)

    # extract results
    sum_model <- summary(model)
    coefs <- sum_model$coefficients
    ci <- suppressWarnings(stats::confint(model, parm = var, level = 0.95))
    if (!is.matrix(ci)) ci <- t(as.matrix(ci))
    main_row <- coefs[rownames(coefs) == var, , drop = FALSE]

    # OR = exp(beta), CI = exp(CI of beta)
    results[[var]] <- data.frame(
      variable = var,
      OR = round(exp(main_row[, "Estimate"]), 3),
      lower95 = round(exp(ci[, 1]), 3),
      upper95 = round(exp(ci[, 2]), 3),
      pvalue = signif(main_row[, "Pr(>|z|)"], 3),
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Validate regression function inputs
#' @param data The input data.
#' @param main_var The main variables.
#' @param covariates The covariates (can be NULL).
#' @param required_cols Column names that must exist in data.
#' @return NULL (invisibly). Stops with an error if validation fails.
#' @keywords internal
#' @noRd
.validate_regression_inputs <- function(data, main_var, covariates, required_cols) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or data.table.")
  }

  if (!is.character(main_var) || length(main_var) == 0) {
    stop("'main_var' must be a non-empty character vector.")
  }

  all_vars <- c(main_var, covariates, required_cols)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("The following variables are not found in data: %s",
                 paste(missing_vars, collapse = ", ")))
  }

  invisible(NULL)
}