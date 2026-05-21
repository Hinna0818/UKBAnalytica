#' @importFrom stats AIC as.formula coef vcov quantile model.matrix predict
#'   delete.response terms anova qnorm model.frame
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_hline geom_vline
#'   geom_point geom_rect geom_rug annotate labs theme_classic theme element_text
#'   element_line scale_x_continuous scale_y_continuous sec_axis expansion
#'   element_blank
#' @importFrom grid unit
#' @importFrom rlang .data
#' @importFrom survival Surv coxph
#' @importFrom splines ns
NULL

# Internal helpers

.ukb_rcs_label_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  if (p < 0.01)  return(sprintf("= %.3f", p))
  sprintf("= %.3f", p)
}

.ukb_rcs_validate <- function(data, exposure, covariates, model_type, endpoint, outcome) {
  if (!is.data.frame(data) && !inherits(data, "data.table")) {
    stop("'data' must be a data.frame.")
  }
  if (length(exposure) != 1 || !is.character(exposure)) {
    stop("'exposure' must be a single character string.")
  }
  if (!exposure %in% names(data)) {
    stop(sprintf("Exposure column '%s' not found in data.", exposure))
  }
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, names(data))
    if (length(missing_cov) > 0) {
      stop(sprintf("Covariate columns not found: %s", paste(missing_cov, collapse = ", ")))
    }
  }
  model_type <- match.arg(model_type, c("cox", "logistic", "linear"))
  if (model_type == "cox") {
    if (is.null(endpoint) || length(endpoint) != 2) {
      stop("'endpoint' must be a character vector of length 2 for model_type = 'cox'.")
    }
    missing_ep <- setdiff(endpoint, names(data))
    if (length(missing_ep) > 0) {
      stop(sprintf("Endpoint columns not found: %s", paste(missing_ep, collapse = ", ")))
    }
  } else {
    if (is.null(outcome)) {
      stop("'outcome' must be specified for model_type = 'logistic' or 'linear'.")
    }
    if (!outcome %in% names(data)) {
      stop(sprintf("Outcome column '%s' not found in data.", outcome))
    }
  }
  invisible(TRUE)
}

.ukb_rcs_prepare_data <- function(data, exposure, covariates, endpoint, outcome,
                                   trim_quantiles) {
  keep_cols <- exposure
  if (!is.null(covariates)) keep_cols <- c(keep_cols, covariates)
  if (!is.null(endpoint))   keep_cols <- c(keep_cols, endpoint)
  if (!is.null(outcome))    keep_cols <- c(keep_cols, outcome)
  keep_cols <- unique(keep_cols)

  df <- as.data.frame(data)[, keep_cols, drop = FALSE]
  df <- df[stats::complete.cases(df), ]

  if (!is.null(trim_quantiles) && length(trim_quantiles) == 2) {
    qs <- stats::quantile(df[[exposure]], probs = trim_quantiles, na.rm = TRUE)
    df <- df[df[[exposure]] >= qs[1] & df[[exposure]] <= qs[2], ]
  }
  df
}

.ukb_rcs_build_formula <- function(exposure, covariates, model_type, endpoint,
                                    outcome, spline_term) {
  rhs <- spline_term
  if (!is.null(covariates) && length(covariates) > 0) {
    rhs <- paste(c(spline_term, covariates), collapse = " + ")
  }
  lhs <- if (model_type == "cox") {
    sprintf("survival::Surv(%s, %s)", endpoint[1], endpoint[2])
  } else {
    outcome
  }
  stats::as.formula(paste(lhs, "~", rhs))
}

# rms backend 

.rms_datadist_setup <- function(model_data, exposure, ref_value) {
  dd <- rms::datadist(model_data)
  dd[["limits"]]["Adjust to", exposure] <- ref_value
  old_opt <- getOption("datadist")
  list(datadist = dd, old_opt = old_opt)
}

.rms_datadist_teardown <- function(setup) {
  options(datadist = setup$old_opt)
}

.ukb_rcs_aic_rms <- function(model_data, exposure, covariates, model_type,
                               endpoint, outcome, knot_range) {
  vapply(knot_range, function(k) {
    spline_term <- sprintf("rms::rcs(%s, %d)", exposure, k)
    f <- .ukb_rcs_build_formula(exposure, covariates, model_type, endpoint,
                                  outcome, spline_term)
    m <- tryCatch({
      if (model_type == "cox") {
        rms::cph(f, data = model_data, x = TRUE, y = TRUE)
      } else if (model_type == "logistic") {
        rms::lrm(f, data = model_data, x = TRUE, y = TRUE)
      } else {
        rms::ols(f, data = model_data, x = TRUE, y = TRUE)
      }
    }, error = function(e) NULL)
    if (is.null(m)) return(NA_real_)
    tryCatch(AIC(m), error = function(e) NA_real_)
  }, numeric(1))
}

.ukb_rcs_fit_rms <- function(model_data, exposure, covariates, model_type,
                               endpoint, outcome, knots, ref_value, conf_level,
                               grid_size, knot_range) {
  if (!requireNamespace("rms", quietly = TRUE)) {
    stop("Package 'rms' is required for backend = 'rms'. Install with install.packages('rms').")
  }

  # AIC-based knot selection if needed
  if (is.null(knots)) {
    aic_vals <- .ukb_rcs_aic_rms(model_data, exposure, covariates, model_type,
                                   endpoint, outcome, knot_range)
    best_idx <- which.min(aic_vals)
    knots <- knot_range[best_idx]
    aic_table <- data.frame(knots = knot_range, AIC = aic_vals,
                             stringsAsFactors = FALSE)
  } else {
    aic_table <- data.frame(knots = knots, AIC = NA_real_,
                             stringsAsFactors = FALSE)
  }

  # Set up datadist
  dd_setup <- .rms_datadist_setup(model_data, exposure, ref_value)
  options(datadist = dd_setup$datadist)
  on.exit(.rms_datadist_teardown(dd_setup), add = TRUE)

  # Fit final model
  spline_term <- sprintf("rms::rcs(%s, %d)", exposure, knots)
  f <- .ukb_rcs_build_formula(exposure, covariates, model_type, endpoint,
                                outcome, spline_term)
  model <- if (model_type == "cox") {
    rms::cph(f, data = model_data, x = TRUE, y = TRUE)
  } else if (model_type == "logistic") {
    rms::lrm(f, data = model_data, x = TRUE, y = TRUE)
  } else {
    rms::ols(f, data = model_data, x = TRUE, y = TRUE)
  }

  # Extract P values from anova
  aov_mat <- tryCatch({
    aov <- stats::anova(model)
    as.matrix(aov)
  }, error = function(e) NULL)

  p_overall   <- NA_real_
  p_nonlinear <- NA_real_
  if (!is.null(aov_mat)) {
    p_col <- intersect(c("P", "Pr(Chi)"), colnames(aov_mat))
    if (length(p_col) > 0) {
      rn <- rownames(aov_mat)
      exp_row <- grep(exposure, rn, fixed = TRUE)
      nl_row  <- grep("nonlinear", rn, ignore.case = TRUE)
      if (length(exp_row) > 0) p_overall   <- aov_mat[exp_row[1], p_col[1]]
      if (length(nl_row)  > 0) p_nonlinear <- aov_mat[nl_row[1],  p_col[1]]
    }
  }

  # Predict over exposure range
  use_exp <- model_type != "linear"
  pred_args <- c(
    list(model,
         fun     = if (use_exp) exp else NULL,
         ref.zero = TRUE,
         np      = grid_size),
    stats::setNames(list(NA), exposure)
  )
  pred_obj <- do.call(rms::Predict, pred_args)
  pred_df  <- as.data.frame(pred_obj)

  # Standardise column names (rms uses yhat, lower, upper)
  col_map <- c(yhat = "estimate", lower = "lower95", upper = "upper95")
  for (old in names(col_map)) {
    if (old %in% names(pred_df)) names(pred_df)[names(pred_df) == old] <- col_map[old]
  }

  prediction <- data.frame(
    x        = pred_df[[exposure]],
    estimate = pred_df$estimate,
    lower95  = pred_df$lower95,
    upper95  = pred_df$upper95,
    stringsAsFactors = FALSE
  )

  # n and n_event
  n       <- nrow(model_data)
  n_event <- if (model_type == "cox") sum(model_data[[endpoint[2]]], na.rm = TRUE) else NA_integer_

  list(
    model        = model,
    knots        = knots,
    ref          = ref_value,
    n            = n,
    n_event      = n_event,
    p_overall    = p_overall,
    p_nonlinear  = p_nonlinear,
    prediction   = prediction,
    aic_table    = aic_table
  )
}

# ns backend 

.ukb_rcs_lrt_p <- function(null_m, alt_m, model_type) {
  tryCatch({
    if (model_type == "cox") {
      tbl <- anova(null_m, alt_m)
    } else if (model_type == "logistic") {
      tbl <- anova(null_m, alt_m, test = "LRT")
    } else {
      tbl <- anova(null_m, alt_m)
    }
    p_col <- grep("Pr\\(|P\\(", colnames(tbl), value = TRUE)
    if (length(p_col) == 0 || nrow(tbl) < 2) {
      return(NA_real_)
    }
    p_val <- tbl[nrow(tbl), p_col[[1]]]
    if (length(p_val) == 0 || is.null(p_val)) NA_real_ else as.numeric(p_val)
  }, error = function(e) NA_real_)
}

.ukb_rcs_predict_ns <- function(model, model_data, exposure, covariates,
                                  model_type, ref_idx, conf_level) {
  alpha <- 1 - conf_level
  z_crit <- stats::qnorm(1 - alpha / 2)

  # Build model matrix for prediction data
  trms <- stats::delete.response(stats::terms(model))
  mm <- tryCatch(
    stats::model.matrix(trms, data = model_data),
    error = function(e) NULL
  )
  if (is.null(mm)) return(NULL)

  # Cox: remove intercept column
  if (model_type == "cox") {
    int_col <- which(colnames(mm) == "(Intercept)")
    if (length(int_col) > 0) mm <- mm[, -int_col, drop = FALSE]
  }

  beta <- stats::coef(model)
  V    <- stats::vcov(model)

  # Keep only matched columns (in case of NA coefs from rank deficiency)
  valid <- intersect(names(beta), colnames(mm))
  mm   <- mm[, valid, drop = FALSE]
  beta <- beta[valid]
  V    <- V[valid, valid, drop = FALSE]

  lp      <- as.numeric(mm %*% beta)
  lp_ref  <- lp[ref_idx]
  rel_lp  <- lp - lp_ref

  # Variance via delta method on contrasts
  C      <- mm - mm[rep(ref_idx, nrow(mm)), , drop = FALSE]
  se_vec <- sqrt(pmax(0, apply(C, 1, function(d) as.numeric(t(d) %*% V %*% d))))

  if (model_type == "linear") {
    list(
      estimate = rel_lp,
      lower95  = rel_lp - z_crit * se_vec,
      upper95  = rel_lp + z_crit * se_vec
    )
  } else {
    list(
      estimate = exp(rel_lp),
      lower95  = exp(rel_lp - z_crit * se_vec),
      upper95  = exp(rel_lp + z_crit * se_vec)
    )
  }
}

.ukb_rcs_fit_ns <- function(model_data, exposure, covariates, model_type,
                              endpoint, outcome, knots, ref_value, conf_level,
                              grid_size, knot_range) {
  # AIC-based knot selection via df (knots = df + 1 for ns)
  fit_ns_model <- function(k, pred_data = model_data) {
    df_ns <- k - 1L
    spline_term <- sprintf("splines::ns(%s, df = %d)", exposure, df_ns)
    f <- .ukb_rcs_build_formula(exposure, covariates, model_type, endpoint,
                                  outcome, spline_term)
    if (model_type == "cox") {
      survival::coxph(f, data = pred_data)
    } else if (model_type == "logistic") {
      stats::glm(f, data = pred_data, family = stats::binomial())
    } else {
      stats::lm(f, data = pred_data)
    }
  }

  if (is.null(knots)) {
    aic_vals <- vapply(knot_range, function(k) {
      m <- tryCatch(fit_ns_model(k), error = function(e) NULL)
      if (is.null(m)) NA_real_ else tryCatch(AIC(m), error = function(e) NA_real_)
    }, numeric(1))
    knots     <- knot_range[which.min(aic_vals)]
    aic_table <- data.frame(knots = knot_range, AIC = aic_vals, stringsAsFactors = FALSE)
  } else {
    aic_table <- data.frame(knots = knots, AIC = NA_real_, stringsAsFactors = FALSE)
  }

  # Fit final spline model
  spline_model <- fit_ns_model(knots)

  # Linear model (for nonlinear P)
  lin_term <- exposure
  f_lin    <- .ukb_rcs_build_formula(exposure, covariates, model_type, endpoint,
                                      outcome, lin_term)
  linear_model <- tryCatch(
    if (model_type == "cox") survival::coxph(f_lin, data = model_data)
    else if (model_type == "logistic") stats::glm(f_lin, data = model_data, family = stats::binomial())
    else stats::lm(f_lin, data = model_data),
    error = function(e) NULL
  )

  # Null model (no exposure) for p_overall
  cov_rhs <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    "1"
  }
  lhs <- if (model_type == "cox") {
    sprintf("survival::Surv(%s, %s)", endpoint[1], endpoint[2])
  } else {
    outcome
  }
  f_null <- stats::as.formula(paste(lhs, "~", cov_rhs))
  null_model <- tryCatch(
    if (model_type == "cox") survival::coxph(f_null, data = model_data)
    else if (model_type == "logistic") stats::glm(f_null, data = model_data, family = stats::binomial())
    else stats::lm(f_null, data = model_data),
    error = function(e) NULL
  )

  p_overall   <- if (!is.null(null_model)) .ukb_rcs_lrt_p(null_model, spline_model, model_type) else NA_real_
  p_nonlinear <- if (!is.null(linear_model)) .ukb_rcs_lrt_p(linear_model, spline_model, model_type) else NA_real_

  # Build prediction grid over exposure range
  x_seq <- seq(
    min(model_data[[exposure]], na.rm = TRUE),
    max(model_data[[exposure]], na.rm = TRUE),
    length.out = grid_size
  )
  pred_df <- data.frame(x_seq)
  names(pred_df)[1] <- exposure

  # Fill covariates at their medians / modes
  if (!is.null(covariates)) {
    for (cov in covariates) {
      cv <- model_data[[cov]]
      pred_df[[cov]] <- if (is.numeric(cv)) {
        stats::median(cv, na.rm = TRUE)
      } else if (is.factor(cv)) {
        mode_value <- names(which.max(table(cv)))
        factor(mode_value, levels = levels(cv), ordered = is.ordered(cv))
      } else {
        names(which.max(table(cv)))
      }
    }
  }
  if (model_type == "cox") {
    pred_df[[endpoint[1]]] <- stats::median(model_data[[endpoint[1]]], na.rm = TRUE)
    pred_df[[endpoint[2]]] <- 0L
  }

  # Reference index in the grid
  ref_idx <- which.min(abs(x_seq - ref_value))

  ci_list <- .ukb_rcs_predict_ns(
    spline_model, pred_df, exposure, covariates,
    model_type, ref_idx, conf_level
  )
  if (is.null(ci_list)) {
    stop("Failed to compute predictions for ns backend.")
  }

  prediction <- data.frame(
    x        = x_seq,
    estimate = ci_list$estimate,
    lower95  = ci_list$lower95,
    upper95  = ci_list$upper95,
    stringsAsFactors = FALSE
  )

  n       <- nrow(model_data)
  n_event <- if (model_type == "cox") sum(model_data[[endpoint[2]]], na.rm = TRUE) else NA_integer_

  list(
    model        = spline_model,
    knots        = knots,
    ref          = ref_value,
    n            = n,
    n_event      = n_event,
    p_overall    = p_overall,
    p_nonlinear  = p_nonlinear,
    prediction   = prediction,
    aic_table    = aic_table
  )
}


#' Fit a restricted cubic spline exposure-response model
#'
#' @description
#' Fits a restricted cubic spline (RCS) model to characterise nonlinear
#' exposure-response relationships. Supports Cox, logistic, and linear
#' regression. Returns prediction curves, confidence intervals, overall and
#' nonlinear P values, and the AIC-selected knot count. The returned object
#' is passed directly to \code{plot_rcs()} for publication-ready figures.
#'
#' @param data A data.frame containing all required columns.
#' @param exposure Character. Name of the continuous exposure variable.
#' @param covariates Character vector of covariate names, or \code{NULL}.
#' @param model_type One of \code{"cox"}, \code{"logistic"}, \code{"linear"}.
#' @param endpoint Character vector of length 2 giving \code{c(time, status)}
#'   column names. Required when \code{model_type = "cox"}.
#' @param outcome Character. Outcome column name. Required for
#'   \code{"logistic"} and \code{"linear"}.
#' @param knots Integer. Number of knots (3-7). If \code{NULL}, the knot count
#'   with the lowest AIC within \code{knot_range} is chosen automatically.
#' @param knot_range Integer vector of candidate knot counts for AIC selection.
#'   Default \code{3:7}.
#' @param ref Numeric. Reference value for the exposure. If \code{NULL},
#'   \code{ref_quantile} is used.
#' @param ref_quantile Numeric (0-1). Quantile of the exposure used as the
#'   reference when \code{ref} is \code{NULL}. Default \code{0.5} (median).
#' @param conf_level Numeric. Confidence level for intervals. Default \code{0.95}.
#' @param trim_quantiles Numeric vector of length 2. Exposure values outside
#'   these quantiles are excluded before fitting. Default \code{c(0.01, 0.99)}.
#' @param grid_size Integer. Number of points in the prediction grid. Default \code{200}.
#' @param backend One of \code{"rms"} (default, requires the \pkg{rms} package)
#'   or \code{"ns"} (base-R natural cubic splines, no additional dependencies).
#'
#' @return An object of class \code{c("ukb_rcs", "list")} with elements:
#' \describe{
#'   \item{model}{The fitted model object.}
#'   \item{model_type}{Character. One of \code{"cox"}, \code{"logistic"}, \code{"linear"}.}
#'   \item{backend}{Character. \code{"rms"} or \code{"ns"}.}
#'   \item{exposure}{Character. Name of the exposure variable.}
#'   \item{covariates}{Character vector of covariate names.}
#'   \item{endpoint}{Character vector. Cox endpoint columns.}
#'   \item{outcome}{Character. Outcome column name.}
#'   \item{knots}{Integer. Number of knots used.}
#'   \item{ref}{Numeric. Reference exposure value.}
#'   \item{n}{Integer. Number of observations in the fitted model.}
#'   \item{n_event}{Integer. Number of events (Cox only, else \code{NA}).}
#'   \item{p_overall}{Numeric. Overall P value for the exposure term.}
#'   \item{p_nonlinear}{Numeric. P value for the nonlinear component.}
#'   \item{prediction}{data.frame with columns \code{x}, \code{estimate},
#'     \code{lower95}, \code{upper95}.}
#'   \item{distribution}{data.frame with column \code{x} (untrimmed exposure values).}
#'   \item{aic_table}{data.frame with columns \code{knots} and \code{AIC}.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 500
#' df <- data.frame(
#'   time   = rexp(n, 0.05),
#'   status = rbinom(n, 1, 0.4),
#'   no2    = rnorm(n, 40, 10),
#'   age    = rnorm(n, 55, 10),
#'   sex    = rbinom(n, 1, 0.5)
#' )
#' fit <- run_rcs(df, exposure = "no2",
#'                covariates = c("age", "sex"),
#'                model_type = "cox",
#'                endpoint   = c("time", "status"))
#' fit$p_overall
#' fit$aic_table
#' plot_rcs(fit, xlab = "NO2 (ug/m3)", ylab = "HR (95% CI)")
#' }
#'
#' @importFrom stats quantile complete.cases
#' @export
run_rcs <- function(data,
                    exposure,
                    covariates       = NULL,
                    model_type       = c("cox", "logistic", "linear"),
                    endpoint         = NULL,
                    outcome          = NULL,
                    knots            = NULL,
                    knot_range       = 3:7,
                    ref              = NULL,
                    ref_quantile     = 0.5,
                    conf_level       = 0.95,
                    trim_quantiles   = c(0.01, 0.99),
                    grid_size        = 200L,
                    backend          = c("rms", "ns")) {

  model_type <- match.arg(model_type)
  backend    <- match.arg(backend)

  .ukb_rcs_validate(data, exposure, covariates, model_type, endpoint, outcome)

  model_data <- .ukb_rcs_prepare_data(data, exposure, covariates, endpoint,
                                        outcome, trim_quantiles)
  if (nrow(model_data) < 20) {
    stop("Fewer than 20 complete observations remain after trimming.")
  }

  # Reference value
  ref_value <- if (!is.null(ref)) {
    as.numeric(ref)
  } else {
    as.numeric(stats::quantile(model_data[[exposure]], probs = ref_quantile,
                                na.rm = TRUE))
  }

  # Store full (untrimmed but complete) distribution for the plot
  dist_data <- as.data.frame(data)[[exposure]]
  dist_data <- dist_data[!is.na(dist_data)]

  fit_result <- if (backend == "rms") {
    .ukb_rcs_fit_rms(model_data, exposure, covariates, model_type, endpoint,
                      outcome, knots, ref_value, conf_level, grid_size, knot_range)
  } else {
    .ukb_rcs_fit_ns(model_data, exposure, covariates, model_type, endpoint,
                     outcome, knots, ref_value, conf_level, grid_size, knot_range)
  }

  out <- c(
    list(
      model_type   = model_type,
      backend      = backend,
      exposure     = exposure,
      covariates   = covariates,
      endpoint     = endpoint,
      outcome      = outcome,
      distribution = data.frame(x = dist_data, stringsAsFactors = FALSE)
    ),
    fit_result
  )

  class(out) <- c("ukb_rcs", "list")
  out
}

# plot_rcs

#' Plot a restricted cubic spline exposure-response curve
#'
#' @description
#' Produces a publication-ready ggplot2 figure from a \code{ukb_rcs} object
#' returned by \code{\link{run_rcs}}. The main panel shows the estimated
#' effect curve with a 95\% confidence ribbon. An optional distribution layer
#' (histogram, density, or rug) is drawn behind the curve to show exposure
#' density. P values and the knot count are annotated by default.
#'
#' @param x A \code{ukb_rcs} object from \code{\link{run_rcs}}.
#' @param show_distribution Logical. Whether to overlay an exposure distribution
#'   layer. Default \code{TRUE}.
#' @param distribution One of \code{"histogram"} (default), \code{"density"},
#'   or \code{"rug"}.
#' @param show_ref Logical. Whether to mark the reference value with a point.
#'   Default \code{TRUE}.
#' @param show_p Logical. Whether to annotate P-overall and P-nonlinear.
#'   Default \code{TRUE}.
#' @param show_knots Logical. Whether to annotate the knot count. Default \code{TRUE}.
#' @param curve_color Character. Hex color for the main curve and ribbon.
#'   Default \code{"#2166AC"} (deep blue).
#' @param dist_color Character. Fill color for the distribution layer.
#'   Default \code{"#AECDE8"}.
#' @param title Character. Plot title. Default \code{NULL} (no title).
#' @param xlab Character. x-axis label. Default: the exposure variable name.
#' @param ylab Character. y-axis label. Default is chosen from model type.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' fit <- run_rcs(df, exposure = "no2", model_type = "cox",
#'                endpoint = c("time", "status"))
#' plot_rcs(fit, xlab = "NO2 (ug/m3)")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_hline geom_point
#'   geom_rect geom_rug annotate labs theme_classic theme element_text
#'   element_line element_blank scale_y_continuous sec_axis expansion
#' @importFrom grid unit
#' @importFrom rlang .data
#' @export
plot_rcs <- function(x, ...) {
  UseMethod("plot_rcs")
}

#' @rdname plot_rcs
#' @export
plot.ukb_rcs <- function(x, ...) {
  plot_rcs(x, ...)
}

#' @rdname plot_rcs
#' @export
plot_rcs.ukb_rcs <- function(x,
                               show_distribution = TRUE,
                               distribution      = c("histogram", "density", "rug"),
                               show_ref          = TRUE,
                               show_p            = TRUE,
                               show_knots        = TRUE,
                               curve_color       = "#2166AC",
                               dist_color        = "#AECDE8",
                               title             = NULL,
                               xlab              = NULL,
                               ylab              = NULL,
                               ...) {

  distribution <- match.arg(distribution)

  pred <- x$prediction
  dist <- x$distribution$x

  # Axis labels
  if (is.null(xlab)) xlab <- x$exposure
  if (is.null(ylab)) {
    ylab <- switch(x$model_type,
                   cox      = "Hazard ratio (95% CI)",
                   logistic = "Odds ratio (95% CI)",
                   linear   = "Difference (95% CI)")
  }

  # Null reference level
  null_y <- if (x$model_type == "linear") 0 else 1

  # Base plot: CI ribbon + curve
  p <- ggplot(pred, aes(x = .data$x, y = .data$estimate)) +
    geom_ribbon(
      aes(ymin = .data$lower95, ymax = .data$upper95),
      fill  = curve_color,
      alpha = 0.15
    ) +
    geom_hline(
      yintercept = null_y,
      linetype   = "dashed",
      color      = "grey55",
      linewidth  = 0.5
    ) +
    geom_line(
      color     = curve_color,
      linewidth = 0.9
    )

  # Reference point
  if (show_ref) {
    ref_row <- pred[which.min(abs(pred$x - x$ref)), ]
    p <- p + geom_point(
      data  = ref_row,
      color = "#B2182B",
      size  = 2.5,
      shape = 19
    )
  }

  # Distribution layer
  if (show_distribution && length(dist) > 0) {
    p <- .rcs_add_distribution(p, dist, pred, distribution, dist_color)
  }

  # P-value annotation
  if (show_p) {
    p_lab <- paste0(
      "P-overall ", .ukb_rcs_label_p(x$p_overall), "\n",
      "P-nonlinear ", .ukb_rcs_label_p(x$p_nonlinear)
    )
    x_pos <- min(pred$x) + diff(range(pred$x)) * 0.60
    y_max <- max(pred$upper95, na.rm = TRUE)
    y_min <- min(pred$lower95, na.rm = TRUE)
    y_pos <- y_max - diff(c(y_min, y_max)) * 0.04

    p <- p + annotate(
      "text",
      x     = x_pos,
      y     = y_pos,
      label = p_lab,
      hjust = 0,
      vjust = 1,
      size  = 3.0,
      color = "black"
    )
  }

  # Knot annotation
  if (show_knots) {
    k_lab <- sprintf("Knots: %d", x$knots)
    x_kpos <- min(pred$x) + diff(range(pred$x)) * 0.02
    y_max  <- max(pred$upper95, na.rm = TRUE)
    y_min  <- min(pred$lower95, na.rm = TRUE)
    y_kpos <- y_max - diff(c(y_min, y_max)) * 0.04
    p <- p + annotate(
      "text",
      x     = x_kpos,
      y     = y_kpos,
      label = k_lab,
      hjust = 0,
      vjust = 1,
      size  = 3.0,
      color = "grey40"
    )
  }

  # Theme
  p <- p +
    labs(title = title, x = xlab, y = ylab) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      axis.line     = element_line(color = "black", linewidth = 0.4),
      axis.ticks    = element_line(color = "black", linewidth = 0.3),
      panel.grid    = element_blank(),
      axis.text     = element_text(color = "black"),
      axis.title    = element_text(color = "black")
    )

  p
}

# Helper: add distribution layer with secondary axis for histogram/density
.rcs_add_distribution <- function(p, dist_vals, pred, distribution, dist_color) {

  if (distribution == "rug") {
    rug_df <- data.frame(x = dist_vals)
    return(
      p + geom_rug(
        data         = rug_df,
        aes(x = .data$x),
        sides        = "b",
        alpha        = 0.25,
        length       = unit(0.025, "npc"),
        inherit.aes  = FALSE
      )
    )
  }

  # y-range from prediction bounds
  y_min  <- min(pred$lower95, na.rm = TRUE)
  y_max  <- max(pred$upper95, na.rm = TRUE)
  y_span <- y_max - y_min

  if (distribution == "histogram") {
    h <- graphics::hist(dist_vals, breaks = 30, plot = FALSE)
    bar_df <- data.frame(
      xmid  = h$mids,
      count = h$counts,
      w     = diff(h$breaks)[1],
      stringsAsFactors = FALSE
    )
    scale_f <- (y_span * 0.22) / max(bar_df$count)
    y_bot   <- y_min - y_span * 0.08

    bar_df$y_top <- bar_df$count * scale_f + y_bot
    bar_df$y_bot <- y_bot

    p <- p +
      geom_rect(
        data = bar_df,
        aes(
          xmin = .data$xmid - .data$w / 2,
          xmax = .data$xmid + .data$w / 2,
          ymin = .data$y_bot,
          ymax = .data$y_top
        ),
        fill        = dist_color,
        alpha       = 0.55,
        inherit.aes = FALSE
      ) +
      scale_y_continuous(
        expand   = expansion(mult = c(0.1, 0.05)),
        sec.axis = sec_axis(
          transform = ~ (. - y_bot) / scale_f,
          name      = "Count"
        )
      )

  } else if (distribution == "density") {
    dens   <- stats::density(dist_vals, n = 256)
    den_df <- data.frame(x = dens$x, dens = dens$y, stringsAsFactors = FALSE)
    # clip to the observed range
    den_df <- den_df[den_df$x >= min(pred$x) & den_df$x <= max(pred$x), ]

    scale_f <- (y_span * 0.22) / max(den_df$dens)
    y_bot   <- y_min - y_span * 0.08

    den_df$y_bot <- y_bot
    den_df$y_scaled <- den_df$dens * scale_f + y_bot

    p <- p +
      geom_ribbon(
        data = den_df,
        aes(x = .data$x, ymin = .data$y_bot, ymax = .data$y_scaled),
        fill        = dist_color,
        alpha       = 0.50,
        inherit.aes = FALSE
      ) +
      scale_y_continuous(
        expand   = expansion(mult = c(0.1, 0.05)),
        sec.axis = sec_axis(
          transform = ~ (. - y_bot) / scale_f,
          name      = "Density"
        )
      )
  }

  p
}
