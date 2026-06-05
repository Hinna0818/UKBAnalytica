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
#' @param ... Additional arguments passed to \code{coxph()}.
#'
#' @importFrom stats as.formula confint
#' @importFrom survival Surv coxph
#' @return A data.frame with columns: \code{variable}, \code{coef},
#'   \code{se}, \code{z}, \code{HR}, \code{lower95}, \code{upper95},
#'   \code{pvalue}, \code{n}, and \code{n_event}.
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
  stopifnot(length(endpoint) == 2)
  .validate_regression_inputs(data, main_var, covariates, c(endpoint))

  results <- list()

  for (var in main_var) {
    # construct formula
    rhs <- var
    if (!is.null(covariates)) {
      rhs <- paste(c(var, covariates), collapse = " + ")
    }
    formula_str <- paste0("Surv(", endpoint[1], ", ", endpoint[2], ") ~ ", rhs)
    formula_obj <- stats::as.formula(formula_str)

    # fit model
    model <- coxph(formula_obj, data = data, ...)

    # extract results
    sum_model <- summary(model)
    coefs <- sum_model$coefficients
    conf <- sum_model$conf.int
    main_row <- coefs[rownames(coefs) == var, , drop = FALSE]
    conf_row <- conf[rownames(conf) == var, , drop = FALSE]

    results[[var]] <- data.frame(
      variable = var,
      coef = unname(main_row[, "coef"]),
      se = unname(main_row[, "se(coef)"]),
      z = unname(main_row[, "z"]),
      HR = unname(conf_row[, "exp(coef)"]),
      lower95 = unname(conf_row[, "lower .95"]),
      upper95 = unname(conf_row[, "upper .95"]),
      pvalue = unname(main_row[, "Pr(>|z|)"]),
      n = unname(sum_model$n),
      n_event = unname(sum_model$nevent),
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


#' Run multiple generalised linear models
#'
#' @description
#' Fit GLMs for each main variable separately using any \code{stats} family.
#' When \code{covariates} is \code{NULL}, univariate models are fitted.
#' Otherwise, multivariate models are fitted.
#'
#' Quasi-families (\code{quasipoisson}, \code{quasibinomial}) use Wald
#' confidence intervals because profile-likelihood CIs are not available for
#' quasi-likelihood models.  All other families use profile-likelihood CIs via
#' \code{stats::confint}.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param main_var A character vector of main variable names to test.
#' @param family A GLM family.  Accepted forms:
#'   \itemize{
#'     \item A character string naming a \code{stats} family function, e.g.
#'           \code{"poisson"}, \code{"Gamma"}, \code{"gaussian"},
#'           \code{"quasipoisson"}, \code{"quasibinomial"},
#'           \code{"inverse.gaussian"}.
#'     \item A family function, e.g. \code{stats::poisson}.
#'     \item A family object, e.g. \code{stats::poisson(link = "sqrt")}.
#'   }
#' @param outcome A character string specifying the outcome column.
#' @param covariates A character vector of covariate names. Default \code{NULL}.
#' @param ... Additional arguments passed to \code{stats::glm()}.
#'
#' @importFrom stats as.formula glm confint qnorm model.frame
#' @return A data.frame with columns: \code{variable}, \code{family},
#'   \code{link}, \code{beta}, \code{lower95}, \code{upper95}, \code{pvalue},
#'   \code{n}.  For log- or logit-link families \code{exp(beta)} gives the
#'   ratio-scale effect (IRR, rate ratio, etc.).
#'
#' @examples
#' \dontrun{
#' # Poisson regression (IRR = exp(beta))
#' mtcars$count <- round(mtcars$mpg)
#' res <- runmulti_glm(mtcars, main_var = c("hp", "wt"),
#'                     family = "poisson", outcome = "count")
#'
#' # Quasi-Poisson (overdispersed counts, Wald CI)
#' res <- runmulti_glm(mtcars, main_var = c("hp", "wt"),
#'                     family = "quasipoisson", outcome = "count")
#'
#' # Gamma with log link
#' res <- runmulti_glm(mtcars, main_var = c("hp", "wt"),
#'                     family = stats::Gamma(link = "log"), outcome = "mpg")
#' }
#'
#' @export
runmulti_glm <- function(data,
                         main_var,
                         family = "poisson",
                         outcome,
                         covariates = NULL,
                         ...) {

  family_obj <- .parse_glm_family(family)
  .validate_regression_inputs(data, main_var, covariates, outcome)

  is_quasi <- grepl("^quasi", family_obj$family)

  results <- list()

  for (var in main_var) {
    rhs <- var
    if (!is.null(covariates)) {
      rhs <- paste(c(var, covariates), collapse = " + ")
    }
    formula_obj <- stats::as.formula(paste(outcome, "~", rhs))

    model <- stats::glm(formula_obj, data = data, family = family_obj, ...)

    sum_model <- summary(model)
    coefs     <- sum_model$coefficients
    main_row  <- coefs[rownames(coefs) == var, , drop = FALSE]

    if (is_quasi) {
      se_val <- main_row[, "Std. Error"]
      b_val  <- main_row[, "Estimate"]
      z95    <- stats::qnorm(0.975)
      ci     <- matrix(c(b_val - z95 * se_val, b_val + z95 * se_val),
                       nrow = 1,
                       dimnames = list(var, c("2.5 %", "97.5 %")))
    } else {
      ci <- suppressWarnings(
        stats::confint(model, parm = var, level = 0.95)
      )
      if (!is.matrix(ci)) ci <- t(as.matrix(ci))
    }

    pval_col <- if ("Pr(>|t|)" %in% colnames(coefs)) "Pr(>|t|)" else "Pr(>|z|)"

    results[[var]] <- data.frame(
      variable = var,
      family   = family_obj$family,
      link     = family_obj$link,
      beta     = round(main_row[, "Estimate"], 4),
      lower95  = round(ci[, 1], 4),
      upper95  = round(ci[, 2], 4),
      pvalue   = signif(main_row[, pval_col], 3),
      n        = nrow(stats::model.frame(model)),
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Run multiple negative-binomial regression models
#'
#' @description
#' Fit negative-binomial GLMs (\code{MASS::glm.nb}) for each main variable
#' separately.  This is the standard approach for overdispersed count outcomes
#' where the Poisson variance assumption is violated.
#'
#' The overdispersion parameter \code{theta} is estimated per model and
#' reported alongside the effect estimate.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param main_var A character vector of main variable names to test.
#' @param outcome A character string specifying the count outcome column.
#' @param covariates A character vector of covariate names. Default \code{NULL}.
#' @param ... Additional arguments passed to \code{MASS::glm.nb()}.
#'
#' @importFrom stats as.formula confint model.frame
#' @return A data.frame with columns: \code{variable}, \code{IRR},
#'   \code{lower95}, \code{upper95}, \code{pvalue}, \code{theta}, \code{n}.
#'   \code{IRR} is the incidence rate ratio (\code{exp(beta)}).
#'   \code{theta} is the estimated negative-binomial dispersion parameter
#'   (larger values indicate less overdispersion).
#'
#' @examples
#' \dontrun{
#' mtcars$count <- round(mtcars$mpg)
#'
#' # Univariate
#' res <- runmulti_negbin(mtcars, main_var = c("hp", "wt"), outcome = "count")
#'
#' # Adjusted
#' res <- runmulti_negbin(mtcars, main_var = c("hp", "wt"),
#'                        covariates = "cyl", outcome = "count")
#' }
#'
#' @export
runmulti_negbin <- function(data,
                            main_var,
                            outcome,
                            covariates = NULL,
                            ...) {

  .validate_regression_inputs(data, main_var, covariates, outcome)

  results <- list()

  for (var in main_var) {
    rhs <- var
    if (!is.null(covariates)) {
      rhs <- paste(c(var, covariates), collapse = " + ")
    }
    formula_obj <- stats::as.formula(paste(outcome, "~", rhs))

    model <- MASS::glm.nb(formula_obj, data = data, ...)

    sum_model <- summary(model)
    coefs    <- sum_model$coefficients
    main_row <- coefs[rownames(coefs) == var, , drop = FALSE]

    ci <- suppressWarnings(stats::confint(model, parm = var, level = 0.95))
    if (!is.matrix(ci)) ci <- t(as.matrix(ci))

    results[[var]] <- data.frame(
      variable = var,
      IRR      = round(exp(main_row[, "Estimate"]), 4),
      lower95  = round(exp(ci[, 1]), 4),
      upper95  = round(exp(ci[, 2]), 4),
      pvalue   = signif(main_row[, "Pr(>|z|)"], 3),
      theta    = round(model$theta, 4),
      n        = nrow(stats::model.frame(model)),
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Run multiple generalised additive models
#'
#' @description
#' Fit GAMs (\code{mgcv::gam}) for each main variable separately.  By default
#' each main variable enters the model as a penalised thin-plate regression
#' spline \code{s(var)}, allowing non-linear dose-response relationships to be
#' detected.
#'
#' When \code{smooth = TRUE} (default) the returned table reports the smooth
#' term's estimated degrees of freedom (\code{edf}), F-statistic, and p-value —
#' useful for screening whether an association exists and whether it is
#' non-linear (\code{edf > 1}).  When \code{smooth = FALSE} the main variable
#' enters as a parametric linear term and the output mirrors \code{runmulti_glm}
#' (\code{beta}, Wald CI, p-value).
#'
#' @param data A data.frame or data.table containing all variables.
#' @param main_var A character vector of main variable names to test.
#' @param outcome A character string specifying the outcome column.
#' @param covariates A character vector of covariate names added as parametric
#'   linear terms.  Default \code{NULL}.
#' @param smooth Logical.  If \code{TRUE} (default) the main variable is
#'   modelled as \code{s(var)}.  If \code{FALSE} it is treated as a linear
#'   term.
#' @param family A GLM family controlling the response distribution.  Accepts
#'   the same forms as \code{\link{runmulti_glm}}: character string, function,
#'   or family object.  Default \code{"gaussian"}.
#' @param k Integer.  Basis dimension for each smooth term.  \code{-1} (default)
#'   lets \code{mgcv} choose automatically.
#' @param ... Additional arguments passed to \code{mgcv::gam()}.
#'
#' @importFrom stats as.formula model.frame qnorm
#' @return When \code{smooth = TRUE}: a data.frame with columns
#'   \code{variable}, \code{edf}, \code{ref_df}, \code{F}, \code{pvalue},
#'   \code{family}, \code{link}, \code{n}.
#'   When \code{smooth = FALSE}: \code{variable}, \code{beta}, \code{lower95},
#'   \code{upper95}, \code{pvalue}, \code{family}, \code{link}, \code{n}.
#'
#' @examples
#' \dontrun{
#' # Smooth GAM (test non-linearity)
#' res <- runmulti_gam(mtcars, main_var = c("hp", "wt"), outcome = "mpg")
#'
#' # Linear GAM (same as runmulti_glm with family = "gaussian")
#' res <- runmulti_gam(mtcars, main_var = c("hp", "wt"), outcome = "mpg",
#'                     smooth = FALSE)
#'
#' # Poisson GAM for count outcomes
#' mtcars$count <- round(mtcars$mpg)
#' res <- runmulti_gam(mtcars, main_var = "hp", outcome = "count",
#'                     family = "poisson", covariates = "cyl")
#' }
#'
#' @export
runmulti_gam <- function(data,
                         main_var,
                         outcome,
                         covariates = NULL,
                         smooth = TRUE,
                         family = "gaussian",
                         k = -1,
                         ...) {

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required. Install with: install.packages('mgcv')")
  }

  family_obj <- .parse_glm_family(family)
  .validate_regression_inputs(data, main_var, covariates, outcome)

  results <- list()

  for (var in main_var) {
    main_term <- if (smooth) {
      if (k == -1) sprintf("s(%s)", var) else sprintf("s(%s, k = %d)", var, k)
    } else {
      var
    }

    rhs <- main_term
    if (!is.null(covariates)) {
      rhs <- paste(c(main_term, covariates), collapse = " + ")
    }
    formula_obj <- stats::as.formula(paste(outcome, "~", rhs))

    model     <- mgcv::gam(formula_obj, data = data, family = family_obj, ...)
    sum_model <- summary(model)
    n         <- nrow(stats::model.frame(model))

    if (smooth) {
      smt     <- sum_model$s.table
      row_idx <- grep(var, rownames(smt), fixed = TRUE)
      if (length(row_idx) == 0) {
        warning(sprintf("Smooth term for '%s' not found in GAM summary.", var))
        next
      }
      smt_row   <- smt[row_idx[1], , drop = FALSE]
      stat_col  <- if ("F" %in% colnames(smt)) "F" else "Chi.sq"
      stat_val  <- round(smt_row[, stat_col], 3)

      results[[var]] <- data.frame(
        variable  = var,
        edf       = round(smt_row[, "edf"],    3),
        ref_df    = round(smt_row[, "Ref.df"], 3),
        stat      = stat_val,
        stat_type = stat_col,
        pvalue    = signif(smt_row[, "p-value"], 3),
        family    = family_obj$family,
        link      = family_obj$link,
        n         = n,
        stringsAsFactors = FALSE
      )
    } else {
      pt       <- sum_model$p.table
      main_row <- pt[rownames(pt) == var, , drop = FALSE]
      b        <- main_row[, "Estimate"]
      se       <- main_row[, "Std. Error"]
      z95      <- stats::qnorm(0.975)
      pval_col <- if ("Pr(>|t|)" %in% colnames(pt)) "Pr(>|t|)" else "Pr(>|z|)"

      results[[var]] <- data.frame(
        variable = var,
        beta     = round(b, 4),
        lower95  = round(b - z95 * se, 4),
        upper95  = round(b + z95 * se, 4),
        pvalue   = signif(main_row[, pval_col], 3),
        family   = family_obj$family,
        link     = family_obj$link,
        n        = n,
        stringsAsFactors = FALSE
      )
    }
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


#' Run a regression model (unified interface)
#'
#' @description
#' A unified wrapper around \code{runmulti_cox}, \code{runmulti_lm},
#' \code{runmulti_logit}, \code{runmulti_glm}, \code{runmulti_negbin}, and
#' \code{runmulti_gam}.  Select the model family with \code{type}.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param main_var A character vector of main variable names to test.
#' @param type One of \code{"cox"}, \code{"lm"}, \code{"logit"}, \code{"glm"},
#'   \code{"negbin"}, or \code{"gam"}.
#' @param outcome For all types except \code{"cox"}: a character string naming
#'   the outcome column.
#' @param endpoint For \code{"cox"}: a character vector of length 2
#'   \code{c("time", "status")}.  Ignored for other types.
#' @param covariates A character vector of covariate names. Default \code{NULL}.
#' @param covariate_sets Optional named list of covariate sets for nested
#'   epidemiological models. Each element must be \code{NULL} or a character
#'   vector of covariate names. When supplied, \code{run_regression()} runs the
#'   same exposure-outcome model once per covariate set and returns one stacked
#'   table with a \code{model} column.
#' @param family For \code{"glm"} and \code{"gam"}: the model family.  Accepts
#'   a character string, function, or family object.  Default \code{"poisson"}
#'   for \code{"glm"} and \code{"gaussian"} for \code{"gam"}.  See
#'   \code{\link{runmulti_glm}} for details.
#' @param smooth For \code{"gam"} only: logical.  Use a smooth spline term
#'   (\code{TRUE}, default) or a linear term (\code{FALSE}).
#' @param ... Additional arguments forwarded to the underlying fitting function.
#'
#' @return A data.frame whose columns depend on \code{type}:
#' \describe{
#'   \item{cox}{variable, coef, se, z, HR, lower95, upper95, pvalue, n, n_event}
#'   \item{lm}{variable, beta, lower95, upper95, pvalue}
#'   \item{logit}{variable, OR, lower95, upper95, pvalue}
#'   \item{glm}{variable, family, link, beta, lower95, upper95, pvalue, n}
#'   \item{negbin}{variable, IRR, lower95, upper95, pvalue, theta, n}
#'   \item{gam (smooth)}{variable, edf, ref_df, F, pvalue, family, link, n}
#'   \item{gam (linear)}{variable, beta, lower95, upper95, pvalue, family, link, n}
#' }
#'
#' @examples
#' \dontrun{
#' library(survival)
#' run_regression(lung, main_var = c("age", "sex"),
#'                type = "cox", endpoint = c("time", "status"))
#'
#' run_regression(mtcars, main_var = c("hp", "wt"),
#'                type = "lm", outcome = "mpg")
#'
#' mtcars$am_bin <- as.integer(mtcars$am == 1)
#' run_regression(mtcars, main_var = c("hp", "wt"),
#'                type = "logit", outcome = "am_bin")
#'
#' run_regression(mtcars, main_var = "hp", type = "lm", outcome = "mpg",
#'                covariate_sets = list(
#'                  crude = NULL,
#'                  model1 = c("cyl"),
#'                  model2 = c("cyl", "wt")
#'                ))
#'
#' mtcars$count <- round(mtcars$mpg)
#' run_regression(mtcars, main_var = c("hp", "wt"),
#'                type = "glm", family = "poisson", outcome = "count")
#'
#' run_regression(mtcars, main_var = c("hp", "wt"),
#'                type = "negbin", outcome = "count")
#'
#' run_regression(mtcars, main_var = c("hp", "wt"),
#'                type = "gam", outcome = "mpg")
#' }
#'
#' @export
run_regression <- function(data,
                           main_var,
                           type = c("cox", "lm", "logit", "glm",
                                    "negbin", "gam"),
                           outcome  = NULL,
                           endpoint = c("time", "status"),
                           covariates = NULL,
                           covariate_sets = NULL,
                           family = NULL,
                           smooth = TRUE,
                           ...) {
  type <- match.arg(type)

  if (!is.null(covariate_sets)) {
    if (!is.null(covariates)) {
      stop("Use either 'covariates' or 'covariate_sets', not both.", call. = FALSE)
    }
    if (!is.list(covariate_sets)) {
      stop("'covariate_sets' must be a named list.", call. = FALSE)
    }
    if (length(covariate_sets) == 0) {
      stop("'covariate_sets' must contain at least one covariate set.", call. = FALSE)
    }
    if (is.null(names(covariate_sets)) || any(!nzchar(names(covariate_sets)))) {
      names(covariate_sets) <- paste0("model", seq_along(covariate_sets))
    }

    nested_results <- vector("list", length(covariate_sets))
    for (i in seq_along(covariate_sets)) {
      set_name <- names(covariate_sets)[[i]]
      set_covariates <- covariate_sets[[i]]
      if (!is.null(set_covariates) && !is.character(set_covariates)) {
        stop(sprintf(
          "Covariate set '%s' must be NULL or a character vector.",
          set_name
        ), call. = FALSE)
      }

      res <- run_regression(
        data = data,
        main_var = main_var,
        type = type,
        outcome = outcome,
        endpoint = endpoint,
        covariates = set_covariates,
        family = family,
        smooth = smooth,
        ...
      )
      res$model <- set_name
      res$covariates <- if (length(set_covariates) == 0) {
        NA_character_
      } else {
        paste(set_covariates, collapse = ", ")
      }
      res$n_covariates <- length(set_covariates)
      first_cols <- c("model", "covariates", "n_covariates")
      nested_results[[i]] <- res[, c(first_cols, setdiff(names(res), first_cols)), drop = FALSE]
    }

    out <- do.call(rbind, nested_results)
    rownames(out) <- NULL
    class(out) <- c("ukb_nested_regression", class(out))
    return(out)
  }

  .require_outcome <- function() {
    if (is.null(outcome))
      stop(sprintf("'outcome' must be specified when type = '%s'.", type))
  }

  switch(type,
    cox = {
      runmulti_cox(data = data, main_var = main_var,
                   covariates = covariates, endpoint = endpoint, ...)
    },
    lm = {
      .require_outcome()
      runmulti_lm(data = data, main_var = main_var,
                  covariates = covariates, outcome = outcome, ...)
    },
    logit = {
      .require_outcome()
      runmulti_logit(data = data, main_var = main_var,
                     covariates = covariates, outcome = outcome, ...)
    },
    glm = {
      .require_outcome()
      fam <- if (is.null(family)) "poisson" else family
      runmulti_glm(data = data, main_var = main_var,
                   family = fam, covariates = covariates,
                   outcome = outcome, ...)
    },
    negbin = {
      .require_outcome()
      runmulti_negbin(data = data, main_var = main_var,
                      covariates = covariates, outcome = outcome, ...)
    },
    gam = {
      .require_outcome()
      fam <- if (is.null(family)) "gaussian" else family
      runmulti_gam(data = data, main_var = main_var,
                   outcome = outcome, covariates = covariates,
                   smooth = smooth, family = fam, ...)
    }
  )
}


#' Parse GLM family argument
#' @param family Character string, function, or family object.
#' @return A family object.
#' @keywords internal
#' @noRd
.parse_glm_family <- function(family) {
  if (inherits(family, "family")) return(family)
  if (is.function(family))       return(family())
  if (is.character(family)) {
    fn <- tryCatch(
      get(family, envir = asNamespace("stats"), inherits = FALSE),
      error = function(e) {
        stop(sprintf("Unknown GLM family: '%s'. Must be a name in stats::",
                     family))
      }
    )
    if (!is.function(fn)) {
      stop(sprintf("'%s' is not a function in the stats package.", family))
    }
    return(fn())
  }
  stop("'family' must be a character string, a function, or a family object.")
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
