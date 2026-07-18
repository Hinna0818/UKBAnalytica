#' Survival Machine Learning Module
#'
#' @description
#' Machine learning models for survival analysis including
#' random survival forests and gradient boosting survival models.
#'
#' @name ml_survival
#' @keywords internal
#' @importFrom survival Surv coxph basehaz concordance
NULL

# Main Survival ML Function

#' Train Survival Machine Learning Model
#'
#' @description
#' Deprecated legacy interface for training machine learning models for
#' survival analysis. New analyses should use
#' \code{\link{ukb_ml_survival_workflow}}, which freezes the test set before
#' feature selection, tuning, final refit, and evaluation.
#'
#' @param formula Survival formula (e.g., Surv(time, event) ~ x1 + x2)
#' @param data Data frame
#' @param model Model type: "rsf" (random survival forest), "gbm_surv" (gradient boosting),
#'   "coxnet" (regularized Cox)
#' @param split_ratio Train/test split ratio (default 0.8)
#' @param seed Random seed
#' @param params List of model-specific parameters
#' @param verbose Print progress
#' @param ... Additional arguments
#'
#' @return A ukb_ml_surv object containing:
#' \itemize{
#'   \item model: Fitted survival model
#'   \item c_index: Harrell's C-index on test data
#'   \item train_data, test_data: Split datasets
#' }
#'
#' @export
ukb_ml_survival <- function(formula,
                            data,
                            model = c("rsf", "gbm_surv", "coxnet"),
                            split_ratio = 0.8,
                            seed = NULL,
                            params = list(),
                            verbose = TRUE,
                            ...) {
  if (exists(".ukb_ml_deprecated", mode = "function")) {
    .ukb_ml_deprecated("ukb_ml_survival", "ukb_ml_survival_workflow()")
  }
  
  model <- match.arg(model)
  
  .check_ml_package("survival")
  
  if (!is.null(seed)) set.seed(seed)
  
  # Parse formula
  terms_obj <- terms(formula, data = data)
  response_vars <- all.vars(formula[[2]])  # Surv(time, event)
  predictor_vars <- attr(terms_obj, "term.labels")
  
  # Handle "." notation
  if (length(predictor_vars) == 1 && predictor_vars == ".") {
    predictor_vars <- setdiff(names(data), response_vars)
  }
  
  # Get time and event variables
  time_var <- response_vars[1]
  event_var <- response_vars[2]
  
  # Remove incomplete cases
  all_vars <- c(time_var, event_var, predictor_vars)
  complete_idx <- complete.cases(data[, all_vars, drop = FALSE])
  clean_data <- data[complete_idx, , drop = FALSE]
  
  if (nrow(clean_data) < nrow(data)) {
    message(sprintf("Removed %d rows with missing values", nrow(data) - nrow(clean_data)))
  }
  
  # Split data
  n <- nrow(clean_data)
  train_idx <- sample(n, round(n * split_ratio))
  test_idx <- setdiff(seq_len(n), train_idx)
  
  train_data <- clean_data[train_idx, ]
  test_data <- clean_data[test_idx, ]
  
  if (verbose) {
    message(sprintf("Training %s survival model", model))
    message(sprintf("Train: %d, Test: %d observations", nrow(train_data), nrow(test_data)))
  }
  
  # Fit model
  fitted <- switch(model,
    rsf = .fit_rsf(formula, train_data, params, verbose, ...),
    gbm_surv = .fit_gbm_surv(formula, train_data, time_var, event_var, predictor_vars, params, verbose, ...),
    coxnet = .fit_coxnet(formula, train_data, time_var, event_var, predictor_vars, params, verbose, ...)
  )
  
  # Create result object
  result <- list(
    model = fitted,
    model_type = model,
    formula = formula,
    time_var = time_var,
    event_var = event_var,
    predictors = predictor_vars,
    train_data = train_data,
    test_data = test_data,
    train_idx = train_idx,
    test_idx = test_idx,
    seed = seed,
    call = match.call()
  )
  
  class(result) <- "ukb_ml_surv"
  
  # Calculate C-index
  result$c_index <- .calculate_c_index(result, test_data)
  
  if (verbose) {
    message(sprintf("Test C-index: %.3f", result$c_index))
  }
  
  result
}

# Model Fitting Functions

#' Fit Random Survival Forest
#' @keywords internal
.fit_rsf <- function(formula, data, params, verbose, ...) {
  .check_ml_package("randomForestSRC")
  
  default_params <- list(
    ntree = 500,
    mtry = NULL,
    nodesize = 15,
    importance = TRUE
  )
  
  model_params <- modifyList(default_params, params)
  
  do.call(randomForestSRC::rfsrc, c(
    list(formula = formula, data = data),
    model_params
  ))
}

#' Fit GBM Survival
#' @keywords internal
.fit_gbm_surv <- function(formula, data, time_var, event_var, predictor_vars, params, verbose, ...) {
  .check_ml_package("gbm")
  
  default_params <- list(
    distribution = "coxph",
    n.trees = 500,
    interaction.depth = 3,
    shrinkage = 0.01,
    n.minobsinnode = 10,
    cv.folds = 0
  )
  
  model_params <- modifyList(default_params, params)
  
  # GBM requires Surv object as response
  gbm_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ ",
    paste(predictor_vars, collapse = " + ")
  ))
  
  do.call(gbm::gbm, c(
    list(formula = gbm_formula, data = data),
    model_params
  ))
}

#' Fit Cox with Elastic Net
#' @keywords internal
.fit_coxnet <- function(formula, data, time_var, event_var, predictor_vars, params, verbose, ...) {
  .check_ml_package("glmnet")
  
  default_params <- list(
    family = "cox",
    alpha = 1  # LASSO
  )
  
  model_params <- modifyList(default_params, params)
  
  # Prepare X matrix and Surv object
  X <- model.matrix(~ . - 1, data = data[, predictor_vars, drop = FALSE])
  y <- Surv(data[[time_var]], data[[event_var]])
  
  do.call(glmnet::cv.glmnet, c(
    list(x = X, y = y),
    model_params
  ))
}

# Prediction Function

#' Predict from Survival ML Model
#'
#' @description
#' Generate predictions from a survival ML model or survival ML workflow.
#'
#' @param object A \code{ukb_ml_survival_workflow},
#'   \code{ukb_ml_survival_final}, or legacy \code{ukb_ml_surv} object.
#' @param newdata Optional new data
#' @param times Time points for survival prediction
#' @param type Prediction type: "risk", "survival", "chf" (cumulative hazard)
#' @param ... Additional arguments
#'
#' @return Matrix of predictions (observations x time points)
#'
#' @export
ukb_ml_survival_predict <- function(object,
                                    newdata = NULL,
                                    times = c(1, 3, 5, 10),
                                    type = c("survival", "risk", "chf"),
                                    ...) {
  
  type <- match.arg(type)

  if (inherits(object, "ukb_ml_survival_workflow")) {
    if (is.null(object$final_model)) {
      stop("The workflow does not contain a final model.", call. = FALSE)
    }
    if (is.null(newdata)) newdata <- object$split$test
    return(.sml_predict_core(object$final_model, newdata = newdata, times = times, type = type))
  }

  if (inherits(object, "ukb_ml_survival_final")) {
    if (is.null(newdata)) newdata <- object$train_data
    return(.sml_predict_core(object, newdata = newdata, times = times, type = type))
  }
  
  if (is.null(newdata)) {
    newdata <- object$test_data
  }
  
  model <- object$model
  model_type <- object$model_type
  
  pred <- switch(model_type,
    rsf = .predict_rsf(model, newdata, times, type),
    gbm_surv = .predict_gbm_surv(model, newdata, times, type, object),
    coxnet = .predict_coxnet(model, newdata, times, type, object)
  )
  
  pred
}

#' @keywords internal
.predict_rsf <- function(model, newdata, times, type) {
  pred <- predict(model, newdata = newdata)
  
  # Get survival probabilities at specified times
  time_idx <- sapply(times, function(t) {
    which.min(abs(pred$time.interest - t))
  })
  
  if (type == "survival") {
    pred$survival[, time_idx, drop = FALSE]
  } else if (type == "chf") {
    pred$chf[, time_idx, drop = FALSE]
  } else {
    # Risk scores (higher = worse prognosis)
    pred$predicted
  }
}

#' @keywords internal
.predict_gbm_surv <- function(model, newdata, times, type, object) {
  .check_ml_package("gbm")
  
  # GBM returns relative risk
  risk <- predict(model, newdata = newdata, n.trees = model$n.trees, type = "link")
  
  if (type == "risk") {
    return(risk)
  }
  
  # For survival probabilities, need baseline hazard estimation
  # Using Cox model baseline
  .check_ml_package("survival")
  
  train_data <- object$train_data
  train_risk <- predict(model, newdata = train_data, n.trees = model$n.trees, type = "link")
  
  # Fit baseline hazard
  surv_obj <- Surv(train_data[[object$time_var]], train_data[[object$event_var]])
  basehaz_fit <- basehaz(coxph(surv_obj ~ offset(train_risk), data = train_data))
  
  # Interpolate baseline hazard at requested times
  H0 <- approx(basehaz_fit$time, basehaz_fit$hazard, times, rule = 2)$y
  
  # Calculate survival or CHF
  surv_mat <- matrix(NA, nrow = length(risk), ncol = length(times))
  
  for (i in seq_along(times)) {
    if (type == "survival") {
      surv_mat[, i] <- exp(-H0[i] * exp(risk))
    } else {
      surv_mat[, i] <- H0[i] * exp(risk)
    }
  }
  
  colnames(surv_mat) <- paste0("t", times)
  surv_mat
}

#' @keywords internal
.predict_coxnet <- function(model, newdata, times, type, object) {
  .check_ml_package("glmnet")
  .check_ml_package("survival")
  
  X <- model.matrix(~ . - 1, data = newdata[, object$predictors, drop = FALSE])
  risk <- as.numeric(predict(model, newx = X, s = "lambda.min", type = "link"))
  
  if (type == "risk") {
    return(risk)
  }
  
  # Estimate baseline hazard from training data
  train_X <- model.matrix(~ . - 1, data = object$train_data[, object$predictors, drop = FALSE])
  train_risk <- as.numeric(predict(model, newx = train_X, s = "lambda.min", type = "link"))
  
  surv_obj <- Surv(object$train_data[[object$time_var]], object$train_data[[object$event_var]])
  basehaz_fit <- basehaz(coxph(surv_obj ~ offset(train_risk), data = object$train_data))
  
  H0 <- approx(basehaz_fit$time, basehaz_fit$hazard, times, rule = 2)$y
  
  surv_mat <- matrix(NA, nrow = length(risk), ncol = length(times))
  
  for (i in seq_along(times)) {
    if (type == "survival") {
      surv_mat[, i] <- exp(-H0[i] * exp(risk))
    } else {
      surv_mat[, i] <- H0[i] * exp(risk)
    }
  }
  
  colnames(surv_mat) <- paste0("t", times)
  surv_mat
}

# C-index Calculation

#' @keywords internal
.calculate_c_index <- function(object, test_data) {
  .check_ml_package("survival")
  
  # Get risk scores
  risk <- switch(object$model_type,
    rsf = predict(object$model, newdata = test_data)$predicted,
    gbm_surv = predict(object$model, newdata = test_data, n.trees = object$model$n.trees, type = "link"),
    coxnet = {
      X <- model.matrix(~ . - 1, data = test_data[, object$predictors, drop = FALSE])
      as.numeric(predict(object$model, newx = X, s = "lambda.min", type = "link"))
    }
  )
  
  .sml_c_index(
    test_data[[object$time_var]],
    test_data[[object$event_var]],
    risk
  )
}

# Frozen-test Survival ML Workflow -----------------------------------------

#' @keywords internal
.sml_parse_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a survival formula.", call. = FALSE)
  }
  response_vars <- all.vars(formula[[2]])
  if (length(response_vars) < 2L) {
    stop("Survival formula must use Surv(time, event) on the left-hand side.", call. = FALSE)
  }
  time_var <- response_vars[1]
  event_var <- response_vars[2]
  missing_response <- setdiff(c(time_var, event_var), names(data))
  if (length(missing_response) > 0L) {
    stop("Survival response column(s) not found: ", paste(missing_response, collapse = ", "), call. = FALSE)
  }

  terms_obj <- stats::terms(formula, data = data)
  predictors <- attr(terms_obj, "term.labels")
  if (length(predictors) == 1L && predictors == ".") {
    predictors <- setdiff(names(data), c(time_var, event_var))
  }
  if (length(predictors) == 0L) {
    stop("Formula must include at least one predictor.", call. = FALSE)
  }

  missing_predictors <- setdiff(all.vars(stats::delete.response(terms_obj)), names(data))
  if (length(missing_predictors) > 0L) {
    stop("Predictor column(s) not found: ", paste(missing_predictors, collapse = ", "), call. = FALSE)
  }

  list(
    time_var = time_var,
    event_var = event_var,
    predictors = predictors,
    terms = terms_obj
  )
}

#' @keywords internal
.sml_formula_from_features <- function(time_var, event_var, features) {
  if (length(features) == 0L) {
    stop("No features available for survival model fitting.", call. = FALSE)
  }
  stats::as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ ",
    paste(features, collapse = " + ")
  ))
}

#' @keywords internal
.sml_model_frame <- function(formula, data) {
  data <- as.data.frame(data)
  parsed <- .sml_parse_formula(formula, data)
  vars <- unique(c(parsed$time_var, parsed$event_var, all.vars(stats::delete.response(parsed$terms))))
  complete_idx <- stats::complete.cases(data[, vars, drop = FALSE])
  mf <- data[complete_idx, vars, drop = FALSE]

  if (nrow(mf) == 0L) {
    stop("No complete-case rows available for survival ML fitting/evaluation.", call. = FALSE)
  }
  if (any(mf[[parsed$time_var]] <= 0, na.rm = TRUE)) {
    stop("Survival time must be positive.", call. = FALSE)
  }
  event_vals <- sort(unique(mf[[parsed$event_var]][!is.na(mf[[parsed$event_var]])]))
  if (!all(event_vals %in% c(0, 1, FALSE, TRUE))) {
    stop("Survival event column must be coded as 0/1 or FALSE/TRUE.", call. = FALSE)
  }
  mf[[parsed$time_var]] <- as.numeric(mf[[parsed$time_var]])
  mf[[parsed$event_var]] <- as.integer(mf[[parsed$event_var]] == 1L | mf[[parsed$event_var]] == TRUE)

  list(
    data = mf,
    time_var = parsed$time_var,
    event_var = parsed$event_var,
    predictors = parsed$predictors,
    complete_idx = which(complete_idx),
    terms = parsed$terms
  )
}

#' @keywords internal
.sml_model_matrix <- function(formula, data) {
  mf <- .sml_model_frame(formula, data)
  x_terms <- stats::delete.response(stats::terms(formula, data = mf$data))
  X <- stats::model.matrix(x_terms, data = mf$data)
  contrasts <- attr(X, "contrasts")
  assign <- attr(X, "assign")
  if ("(Intercept)" %in% colnames(X)) {
    keep <- colnames(X) != "(Intercept)"
    X <- X[, keep, drop = FALSE]
    assign <- assign[keep]
  }
  storage.mode(X) <- "numeric"
  term_labels <- attr(x_terms, "term.labels")
  feature <- colnames(X)
  positive_assign <- assign > 0
  feature[positive_assign] <- term_labels[assign[positive_assign]]
  feature_map <- data.frame(
    matrix_column = colnames(X),
    feature = feature,
    stringsAsFactors = FALSE
  )

  list(
    X = X,
    y = Surv(mf$data[[mf$time_var]], mf$data[[mf$event_var]]),
    model_frame = mf$data,
    time_var = mf$time_var,
    event_var = mf$event_var,
    predictors = mf$predictors,
    complete_idx = mf$complete_idx,
    x_terms = x_terms,
    contrasts = contrasts,
    xlevels = stats::.getXlevels(x_terms, mf$data),
    feature_map = feature_map
  )
}

#' @keywords internal
.sml_check_model <- function(model) {
  match.arg(model, c("cox", "rsf", "gbm_surv", "coxnet"))
}

#' @keywords internal
.sml_default_param_grid <- function(model) {
  switch(
    model,
    cox = list(.none = 1),
    rsf = list(ntree = c(300, 500), nodesize = c(10, 15)),
    gbm_surv = list(n.trees = c(200, 500), interaction.depth = c(1, 3), shrinkage = c(0.01, 0.05)),
    coxnet = list(alpha = c(0, 0.5, 1))
  )
}

#' @keywords internal
.sml_row_to_params <- function(row) {
  row <- as.list(row)
  row$.score <- NULL
  row$.none <- NULL
  row[!vapply(row, function(x) length(x) == 1L && is.na(x), logical(1))]
}

#' @keywords internal
.sml_param_grid <- function(model, param_grid = NULL, param_space = NULL, search = "grid", n_iter = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(param_grid) && is.null(param_space)) {
    param_grid <- .sml_default_param_grid(model)
  }

  if (search == "grid") {
    grid_source <- if (!is.null(param_grid)) param_grid else param_space
    grid <- if (is.data.frame(grid_source)) {
      grid_source
    } else {
      do.call(expand.grid, c(grid_source, stringsAsFactors = FALSE))
    }
  } else {
    space <- if (!is.null(param_space)) param_space else param_grid
    if (is.null(space)) stop("'param_space' or 'param_grid' is required for random search.", call. = FALSE)
    n_iter <- if (is.null(n_iter)) 10L else as.integer(n_iter)
    sampled <- lapply(space, function(vals) sample(vals, n_iter, replace = TRUE))
    grid <- unique(as.data.frame(sampled, stringsAsFactors = FALSE))
  }

  if (".none" %in% names(grid)) grid$.none <- NULL
  if (ncol(grid) == 0L) grid <- data.frame(.none = 1)
  grid
}

#' @keywords internal
.sml_fit_core <- function(formula,
                          data,
                          model,
                          params = list(),
                          seed = NULL,
                          verbose = TRUE) {
  .check_ml_package("survival")
  model <- .sml_check_model(model)
  if (!is.null(seed)) set.seed(seed)

  prep <- .sml_model_matrix(formula, data)
  mf <- prep$model_frame
  time_var <- prep$time_var
  event_var <- prep$event_var
  predictors <- prep$predictors

  fitted <- switch(
    model,
    cox = {
      default_params <- list(x = TRUE, model = TRUE)
      fit_params <- utils::modifyList(default_params, params)
      do.call(coxph, c(list(formula = formula, data = mf), fit_params))
    },
    rsf = .fit_rsf(formula, mf, params, verbose),
    gbm_surv = .fit_gbm_surv(formula, mf, time_var, event_var, predictors, params, verbose),
    coxnet = {
      .check_ml_package("glmnet")
      default_params <- list(family = "cox", alpha = 1)
      fit_params <- utils::modifyList(default_params, params)
      do.call(glmnet::cv.glmnet, c(list(x = prep$X, y = prep$y), fit_params))
    }
  )

  out <- list(
    fitted_model = fitted,
    model = model,
    formula = formula,
    time_var = time_var,
    event_var = event_var,
    predictors = predictors,
    params = params,
    train_data = mf,
    x_terms = prep$x_terms,
    contrasts = prep$contrasts,
    xlevels = prep$xlevels,
    feature_names = colnames(prep$X),
    feature_map = prep$feature_map,
    train_rows = nrow(mf),
    call = match.call()
  )
  class(out) <- c("ukb_ml_survival_final", "ukb_ml_surv_core")
  out
}

#' @keywords internal
.sml_predict_matrix <- function(object, newdata) {
  newdata <- as.data.frame(newdata)
  required <- unique(all.vars(object$x_terms))
  missing <- setdiff(required, names(newdata))
  if (length(missing) > 0L) {
    stop(
      "Prediction data are missing required variable(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  for (variable in names(object$xlevels)) {
    values <- as.character(newdata[[variable]])
    unknown <- setdiff(unique(values[!is.na(values)]), object$xlevels[[variable]])
    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "Prediction variable '%s' contains unseen level(s): %s.",
          variable,
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    newdata[[variable]] <- factor(values, levels = object$xlevels[[variable]])
  }
  if (any(!stats::complete.cases(newdata[, required, drop = FALSE]))) {
    stop(
      "Survival prediction data contain missing predictor values; provide complete predictors.",
      call. = FALSE
    )
  }

  X <- stats::model.matrix(
    object$x_terms,
    data = newdata,
    contrasts.arg = object$contrasts,
    xlev = object$xlevels
  )
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  }
  missing_columns <- setdiff(object$feature_names, colnames(X))
  if (length(missing_columns) > 0L) {
    zeros <- matrix(0, nrow = nrow(X), ncol = length(missing_columns))
    colnames(zeros) <- missing_columns
    X <- cbind(X, zeros)
  }
  X <- X[, object$feature_names, drop = FALSE]
  storage.mode(X) <- "numeric"
  X
}

#' @keywords internal
.sml_baseline_survival <- function(time, event, risk, times) {
  surv_obj <- Surv(time, event)
  fit <- coxph(surv_obj ~ offset(risk))
  basehaz_fit <- basehaz(fit, centered = FALSE)
  H0 <- stats::approx(basehaz_fit$time, basehaz_fit$hazard, times, rule = 2)$y
  H0
}

#' @keywords internal
.sml_predict_core <- function(object,
                              newdata,
                              times = c(1, 3, 5, 10),
                              type = c("survival", "risk", "chf")) {
  type <- match.arg(type)
  if (!inherits(object, "ukb_ml_survival_final")) {
    stop("'object' must be a ukb_ml_survival_final object.", call. = FALSE)
  }
  newdata <- as.data.frame(newdata)

  risk <- switch(
    object$model,
    cox = as.numeric(stats::predict(object$fitted_model, newdata = newdata, type = "lp", reference = "zero")),
    rsf = {
      pred <- predict(object$fitted_model, newdata = newdata)
      if (type == "risk") return(as.numeric(pred$predicted))
      time_idx <- vapply(times, function(t) which.min(abs(pred$time.interest - t)), integer(1))
      out <- if (type == "survival") pred$survival[, time_idx, drop = FALSE] else pred$chf[, time_idx, drop = FALSE]
      colnames(out) <- paste0("t", times)
      return(out)
    },
    gbm_surv = as.numeric(stats::predict(object$fitted_model, newdata = newdata, n.trees = object$fitted_model$n.trees, type = "link")),
    coxnet = {
      .check_ml_package("glmnet")
      X <- .sml_predict_matrix(object, newdata)
      as.numeric(stats::predict(object$fitted_model, newx = X, s = "lambda.min", type = "link"))
    }
  )

  if (type == "risk") return(risk)

  if (object$model == "cox") {
    basehaz_fit <- basehaz(object$fitted_model, centered = FALSE)
    H0 <- stats::approx(basehaz_fit$time, basehaz_fit$hazard, times, rule = 2)$y
    chf <- outer(exp(risk), H0, "*")
    colnames(chf) <- paste0("t", times)
    if (type == "chf") return(chf)
    return(exp(-chf))
  }

  train_risk <- switch(
    object$model,
    gbm_surv = as.numeric(stats::predict(object$fitted_model, newdata = object$train_data, n.trees = object$fitted_model$n.trees, type = "link")),
    coxnet = {
      .check_ml_package("glmnet")
      train_X <- .sml_predict_matrix(object, object$train_data)
      as.numeric(stats::predict(object$fitted_model, newx = train_X, s = "lambda.min", type = "link"))
    }
  )
  H0 <- .sml_baseline_survival(
    object$train_data[[object$time_var]],
    object$train_data[[object$event_var]],
    train_risk,
    times
  )
  chf <- outer(exp(risk), H0, "*")
  colnames(chf) <- paste0("t", times)
  if (type == "chf") return(chf)
  exp(-chf)
}

#' @keywords internal
.sml_c_index <- function(time, event, risk) {
  .check_ml_package("survival")
  keep <- stats::complete.cases(time, event, risk)
  if (sum(keep) < 2L || sum(event[keep] == 1L) == 0L) return(NA_real_)
  conc <- concordance(
    Surv(time[keep], event[keep]) ~ risk[keep],
    reverse = TRUE
  )
  as.numeric(conc$concordance)
}

#' @keywords internal
.sml_censoring_km <- function(time, event) {
  observed_times <- sort(unique(time))
  censoring_survival <- 1
  rows <- vector("list", length(observed_times))
  for (index in seq_along(observed_times)) {
    current_time <- observed_times[[index]]
    at_risk <- sum(time >= current_time)
    censored <- sum(time == current_time & event == 0L)
    before <- censoring_survival
    if (at_risk > 0L && censored > 0L) {
      censoring_survival <- censoring_survival * (1 - censored / at_risk)
    }
    rows[[index]] <- data.frame(
      time = current_time,
      before = before,
      after = censoring_survival
    )
  }
  do.call(rbind, rows)
}

#' @keywords internal
.sml_censoring_probability <- function(km, times, left_limit = FALSE) {
  vapply(times, function(current_time) {
    eligible <- if (isTRUE(left_limit)) {
      km$time < current_time
    } else {
      km$time <= current_time
    }
    if (!any(eligible)) return(1)
    km$after[max(which(eligible))]
  }, numeric(1))
}

#' @keywords internal
.sml_brier_at_times <- function(time, event, survival_prob, times) {
  time <- as.numeric(time)
  event <- as.integer(event)
  survival_prob <- as.matrix(survival_prob)
  if (nrow(survival_prob) != length(time) || ncol(survival_prob) != length(times)) {
    stop(
      "`survival_prob` must have one row per participant and one column per time.",
      call. = FALSE
    )
  }
  keep <- stats::complete.cases(time, event) & time > 0 & event %in% c(0L, 1L)
  time <- time[keep]
  event <- event[keep]
  survival_prob <- survival_prob[keep, , drop = FALSE]
  if (length(time) == 0L) {
    return(stats::setNames(rep(NA_real_, length(times)), paste0("brier_t", times)))
  }
  if (any(!is.na(survival_prob) &
    (!is.finite(survival_prob) | survival_prob < 0 | survival_prob > 1))) {
    stop("Survival probabilities must be between 0 and 1.", call. = FALSE)
  }

  censoring_km <- .sml_censoring_km(time, event)
  stats <- numeric(length(times))
  names(stats) <- paste0("brier_t", times)
  for (i in seq_along(times)) {
    t <- times[[i]]
    predicted_survival <- survival_prob[, i]
    valid_prediction <- !is.na(predicted_survival)
    event_before_t <- valid_prediction & time <= t & event == 1L
    event_free_at_t <- valid_prediction & time > t
    contribution <- rep(0, length(time))

    if (any(event_before_t)) {
      g_before_event <- .sml_censoring_probability(
        censoring_km,
        time[event_before_t],
        left_limit = TRUE
      )
      if (any(g_before_event <= 0)) {
        stats[[i]] <- NA_real_
        next
      }
      contribution[event_before_t] <-
        predicted_survival[event_before_t]^2 / g_before_event
    }

    if (any(event_free_at_t)) {
      g_at_t <- .sml_censoring_probability(censoring_km, t)
      if (!is.finite(g_at_t) || g_at_t <= 0) {
        stats[[i]] <- NA_real_
        next
      }
      contribution[event_free_at_t] <-
        (1 - predicted_survival[event_free_at_t])^2 / g_at_t
    }
    stats[[i]] <- sum(contribution) / sum(valid_prediction)
  }
  stats
}

#' Standardize Manual Survival ML Train/Test Splits
#'
#' @description
#' Register user-provided survival train/test datasets as a frozen split object.
#' This is the survival analogue of \code{\link{ukb_ml_as_split}}.
#'
#' @param train_data Training/development data.
#' @param test_data Frozen test data.
#' @param validation_data Optional validation data.
#' @param time Survival time column.
#' @param event Event indicator column coded 0/1.
#' @param id_col Optional participant ID column used to check overlap.
#' @param check_overlap Logical. Check duplicated and overlapping IDs.
#'
#' @return A \code{ukb_ml_survival_split} object.
#'
#' @export
ukb_ml_survival_as_split <- function(train_data,
                                     test_data,
                                     validation_data = NULL,
                                     time,
                                     event,
                                     id_col = NULL,
                                     check_overlap = TRUE) {
  split <- ukb_ml_as_split(
    train_data = train_data,
    test_data = test_data,
    validation_data = validation_data,
    id_col = id_col,
    check_overlap = check_overlap,
    outcome = event,
    outcome_type = "binary"
  )
  for (nm in c("train", "test", if (!is.null(split$validation)) "validation")) {
    obj <- split[[nm]]
    missing_cols <- setdiff(c(time, event), names(obj))
    if (length(missing_cols) > 0L) {
      stop(sprintf("%s data is missing survival column(s): %s", nm, paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
  }
  split$time_var <- time
  split$event_var <- event
  split$outcome <- event
  split$outcome_type <- "survival"
  class(split) <- c("ukb_ml_survival_split", "ukb_ml_split")
  split
}

#' Split Data into Frozen Survival ML Train/Test Sets
#'
#' @description
#' Creates a frozen train/test or train/validation/test split for time-to-event
#' machine learning. Event status is used for stratification by default.
#'
#' @param df A data.frame or data.table.
#' @param time Survival time column.
#' @param event Event indicator column coded 0/1.
#' @param split Either \code{"train_test"} or \code{"train_valid_test"}.
#' @param train_ratio Training proportion.
#' @param validation_ratio Validation proportion for train/validation/test.
#' @param test_ratio Test proportion.
#' @param stratify_by \code{"event"}, \code{"custom"}, \code{"none"}, or a
#'   column name.
#' @param stratify_col Column used when \code{stratify_by = "custom"}.
#' @param seed Optional random seed.
#' @param verbose Logical. Print split summary.
#'
#' @return A \code{ukb_ml_survival_split} object.
#'
#' @export
ukb_ml_survival_split_data <- function(df,
                                       time,
                                       event,
                                       split = c("train_test", "train_valid_test"),
                                       train_ratio = 0.7,
                                       validation_ratio = 0.1,
                                       test_ratio = 0.2,
                                       stratify_by = c("event", "custom", "none"),
                                       stratify_col = NULL,
                                       seed = NULL,
                                       verbose = TRUE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame or data.table.", call. = FALSE)
  }
  missing_cols <- setdiff(c(time, event), names(df))
  if (length(missing_cols) > 0L) {
    stop("Survival column(s) not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  stratify_by <- stratify_by[[1]]
  ml_stratify <- if (stratify_by == "event") {
    "outcome"
  } else if (stratify_by %in% c("custom", "none")) {
    stratify_by
  } else {
    stratify_col <- stratify_by
    "custom"
  }

  split_obj <- ukb_ml_split_data(
    df = df,
    outcome = event,
    outcome_type = "binary",
    split = match.arg(split),
    train_ratio = train_ratio,
    validation_ratio = validation_ratio,
    test_ratio = test_ratio,
    stratify_by = ml_stratify,
    stratify_col = stratify_col,
    seed = seed,
    verbose = FALSE
  )
  split_obj$time_var <- time
  split_obj$event_var <- event
  split_obj$outcome_type <- "survival"
  split_obj$split_info$stratify_by <- stratify_by
  class(split_obj) <- c("ukb_ml_survival_split", "ukb_ml_split")

  if (isTRUE(verbose)) {
    event_rate <- mean(split_obj$train[[event]] == 1, na.rm = TRUE)
    message(sprintf(
      "[ukb_ml_survival_split_data] train=%d, validation=%d, test=%d, train_event_rate=%.3f",
      nrow(split_obj$train),
      if (!is.null(split_obj$validation)) nrow(split_obj$validation) else 0L,
      nrow(split_obj$test),
      event_rate
    ))
  }
  split_obj
}

#' Select Features for Survival ML Workflows
#'
#' @description
#' Performs optional feature selection using only training data. The test set is
#' never used.
#'
#' @param split A \code{ukb_ml_survival_split} object.
#' @param formula Survival formula.
#' @param method \code{"none"}, \code{"filter"}, or \code{"glmnet"}.
#' @param max_features Optional maximum number of selected features.
#' @param seed Optional random seed.
#' @param verbose Logical.
#'
#' @return A \code{ukb_ml_survival_feature} object.
#'
#' @export
ukb_ml_survival_feature_select <- function(split,
                                           formula,
                                           method = c("none", "filter", "glmnet"),
                                           max_features = NULL,
                                           seed = NULL,
                                           verbose = TRUE) {
  if (!inherits(split, "ukb_ml_survival_split")) {
    stop("'split' must be a ukb_ml_survival_split object.", call. = FALSE)
  }
  method <- match.arg(method)
  parsed <- .sml_parse_formula(formula, split$train)
  if (!is.null(seed)) set.seed(seed)

  features <- parsed$predictors
  selected <- features
  status <- data.frame(feature = features, status = "Selected", stringsAsFactors = FALSE)
  info <- list()

  if (method == "filter") {
    scores <- vapply(features, function(v) {
      f <- .sml_formula_from_features(parsed$time_var, parsed$event_var, v)
      fit <- try(coxph(f, data = as.data.frame(split$train)), silent = TRUE)
      if (inherits(fit, "try-error")) return(0)
      z <- summary(fit)$wald["test"]
      ifelse(is.finite(z), as.numeric(z), 0)
    }, numeric(1))
    scores[!is.finite(scores)] <- 0
    status <- data.frame(feature = names(scores), score = as.numeric(scores), status = "Ranked", stringsAsFactors = FALSE)
    status <- status[order(status$score, decreasing = TRUE), , drop = FALSE]
    selected <- status$feature
  } else if (method == "glmnet") {
    .check_ml_package("glmnet")
    fit <- .sml_fit_core(formula, split$train, model = "coxnet", seed = seed, verbose = FALSE)
    coefs <- as.matrix(stats::coef(fit$fitted_model, s = "lambda.min"))
    nonzero <- rownames(coefs)[coefs[, 1] != 0]
    selected <- unique(fit$feature_map$feature[fit$feature_map$matrix_column %in% nonzero])
    selected <- intersect(features, selected)
    status <- data.frame(feature = features, status = ifelse(features %in% selected, "Selected", "Rejected"), stringsAsFactors = FALSE)
    info$glmnet <- fit$fitted_model
  }

  if (!is.null(max_features) && length(selected) > max_features) {
    selected <- selected[seq_len(max_features)]
  }
  if (length(selected) == 0L) {
    stop("Survival feature selection returned zero features.", call. = FALSE)
  }

  out <- list(
    method = method,
    selected_features = selected,
    selected_status = status,
    formula = .sml_formula_from_features(parsed$time_var, parsed$event_var, selected),
    time_var = parsed$time_var,
    event_var = parsed$event_var,
    feature_info = info
  )
  class(out) <- "ukb_ml_survival_feature"
  if (isTRUE(verbose)) {
    message(sprintf("[ukb_ml_survival_feature_select] method=%s, selected=%d feature(s)", method, length(selected)))
  }
  out
}

#' Tune Survival ML Hyperparameters Without Touching the Test Set
#'
#' @description
#' Tunes survival ML models using validation data or cross-validation inside the
#' training set. The frozen test set is never used.
#'
#' @param split A \code{ukb_ml_survival_split} object.
#' @param formula Survival formula.
#' @param model \code{"cox"}, \code{"rsf"}, \code{"gbm_surv"}, or
#'   \code{"coxnet"}.
#' @param search \code{"grid"} or \code{"random"}.
#' @param param_grid List or data.frame of candidate parameters.
#' @param param_space Parameter space for random search.
#' @param n_iter Number of random-search iterations.
#' @param resampling \code{"cv"} or \code{"validation"}.
#' @param folds Number of CV folds.
#' @param metric Currently \code{"c_index"}.
#' @param maximize Logical. Whether higher metric values are better.
#' @param seed Optional random seed.
#' @param verbose Logical.
#' @param ... Reserved for future extensions.
#'
#' @return A \code{ukb_ml_survival_tune} object.
#'
#' @export
ukb_ml_survival_tune <- function(split,
                                 formula,
                                 model,
                                 search = c("grid", "random"),
                                 param_grid = NULL,
                                 param_space = NULL,
                                 n_iter = NULL,
                                 resampling = c("cv", "validation"),
                                 folds = 5,
                                 metric = "c_index",
                                 maximize = TRUE,
                                 seed = NULL,
                                 verbose = TRUE,
                                 ...) {
  if (!inherits(split, "ukb_ml_survival_split")) {
    stop("'split' must be a ukb_ml_survival_split object.", call. = FALSE)
  }
  model <- .sml_check_model(model)
  search <- match.arg(search)
  resampling <- match.arg(resampling)
  if (resampling == "validation" && is.null(split$validation)) {
    if (isTRUE(verbose)) message("[ukb_ml_survival_tune] No validation set; using CV instead.")
    resampling <- "cv"
  }
  if (!identical(metric, "c_index")) {
    stop("Only metric = 'c_index' is currently supported for survival ML tuning.", call. = FALSE)
  }

  parsed <- .sml_parse_formula(formula, split$train)
  grid <- .sml_param_grid(model, param_grid, param_space, search, n_iter, seed)
  if (!is.null(seed)) set.seed(seed)

  results <- vector("list", nrow(grid))
  best_score_seen <- if (isTRUE(maximize)) -Inf else Inf
  prep_train <- NULL
  fold_idx <- NULL
  if (resampling == "cv") {
    prep_train <- .sml_model_frame(formula, split$train)
    fold_idx <- .mlw_create_folds(
      prep_train$data[[parsed$event_var]],
      folds,
      outcome_type = "binary"
    )
  }

  for (i in seq_len(nrow(grid))) {
    params <- .sml_row_to_params(grid[i, , drop = FALSE])
    scores <- numeric(0)

    if (resampling == "validation") {
      fit <- .sml_fit_core(formula, split$train, model, params = params, seed = seed, verbose = FALSE)
      risk <- .sml_predict_core(fit, split$validation, type = "risk")
      scores <- .sml_c_index(split$validation[[parsed$time_var]], split$validation[[parsed$event_var]], risk)
    } else {
      for (fold in seq_along(fold_idx)) {
        valid_idx <- fold_idx[[fold]]
        train_idx <- setdiff(seq_len(nrow(prep_train$data)), valid_idx)
        fit <- .sml_fit_core(formula, prep_train$data[train_idx, , drop = FALSE], model, params = params, seed = seed, verbose = FALSE)
        risk <- .sml_predict_core(fit, prep_train$data[valid_idx, , drop = FALSE], type = "risk")
        scores <- c(scores, .sml_c_index(
          prep_train$data[[parsed$time_var]][valid_idx],
          prep_train$data[[parsed$event_var]][valid_idx],
          risk
        ))
      }
    }

    score <- mean(scores, na.rm = TRUE)
    results[[i]] <- cbind(grid[i, , drop = FALSE], data.frame(.score = score, stringsAsFactors = FALSE))
    if (isTRUE(verbose)) {
      message(sprintf("[ukb_ml_survival_tune] %s %d/%d: c_index = %.4f", search, i, nrow(grid), score))
    }
    is_better <- if (isTRUE(maximize)) score > best_score_seen else score < best_score_seen
    if (isTRUE(is_better)) best_score_seen <- score
  }

  results_df <- do.call(rbind, results)
  best_idx <- if (isTRUE(maximize)) which.max(results_df$.score) else which.min(results_df$.score)
  out <- list(
    model = model,
    search = search,
    metric = metric,
    maximize = maximize,
    results = results_df,
    best_params = .sml_row_to_params(results_df[best_idx, , drop = FALSE]),
    best_score = results_df$.score[[best_idx]],
    split_info = split$split_info,
    tuning_info = list(
      resampling = resampling,
      folds = folds,
      fold_indices = fold_idx,
      complete_rows = if (!is.null(prep_train)) prep_train$complete_idx else NULL,
      n_candidates = nrow(grid)
    )
  )
  class(out) <- "ukb_ml_survival_tune"
  out
}

#' Refit Final Survival ML Model
#'
#' @description
#' Refits a survival ML model on training plus validation data when available,
#' leaving the frozen test set untouched.
#'
#' @param split A \code{ukb_ml_survival_split} object.
#' @param formula Survival formula.
#' @param model Survival model type.
#' @param best_params Model parameters.
#' @param feature_spec Optional feature-selection result.
#' @param seed Optional random seed.
#' @param verbose Logical.
#' @param ... Additional arguments passed to the fitter.
#'
#' @return A \code{ukb_ml_survival_final} object.
#'
#' @export
ukb_ml_survival_fit_final <- function(split,
                                      formula,
                                      model,
                                      best_params = list(),
                                      feature_spec = NULL,
                                      seed = NULL,
                                      verbose = TRUE,
                                      ...) {
  if (!inherits(split, "ukb_ml_survival_split")) {
    stop("'split' must be a ukb_ml_survival_split object.", call. = FALSE)
  }
  train_dev <- if (!is.null(split$validation)) {
    rbind(as.data.frame(split$train), as.data.frame(split$validation))
  } else {
    as.data.frame(split$train)
  }
  fit <- .sml_fit_core(formula, train_dev, model, params = best_params, seed = seed, verbose = verbose)
  fit$feature_spec <- feature_spec
  fit$selected_features <- if (!is.null(feature_spec)) feature_spec$selected_features else fit$predictors
  fit$fit_data_role <- if (!is.null(split$validation)) "train_plus_validation" else "train"
  fit
}

#' Evaluate Survival ML Once on the Frozen Test Set
#'
#' @description
#' Computes final survival ML metrics on the frozen test set. The primary metric
#' is Harrell's C-index. Time-specific Brier scores use inverse probability of
#' censoring weighting, with the censoring distribution estimated by
#' Kaplan-Meier in the evaluation set.
#'
#' @param object A \code{ukb_ml_survival_final} object.
#' @param split A \code{ukb_ml_survival_split} object.
#' @param times Time points for survival probability prediction.
#' @param verbose Logical.
#' @param ... Additional arguments.
#'
#' @return A \code{ukb_ml_survival_test_eval} object.
#'
#' @export
ukb_ml_survival_evaluate_test <- function(object,
                                          split,
                                          times = c(1, 3, 5, 10),
                                          verbose = TRUE,
                                          ...) {
  if (!inherits(object, "ukb_ml_survival_final")) {
    stop("'object' must be a ukb_ml_survival_final object.", call. = FALSE)
  }
  if (!inherits(split, "ukb_ml_survival_split")) {
    stop("'split' must be a ukb_ml_survival_split object.", call. = FALSE)
  }
  test <- as.data.frame(split$test)
  risk <- .sml_predict_core(object, test, type = "risk")
  survival_prob <- .sml_predict_core(object, test, times = times, type = "survival")
  metrics <- c(
    c_index = .sml_c_index(test[[object$time_var]], test[[object$event_var]], risk),
    .sml_brier_at_times(test[[object$time_var]], test[[object$event_var]], survival_prob, times)
  )
  pred_df <- data.frame(
    time = test[[object$time_var]],
    event = test[[object$event_var]],
    risk = as.numeric(risk),
    survival_prob,
    stringsAsFactors = FALSE
  )

  out <- list(
    metrics = metrics,
    predictions = pred_df,
    times = times,
    evaluated_at = Sys.time(),
    test_rows = nrow(test)
  )
  class(out) <- "ukb_ml_survival_test_eval"
  if (isTRUE(verbose)) {
    message("[ukb_ml_survival_evaluate_test] Final frozen-test metrics:")
    for (nm in names(metrics)) message(sprintf("  %s: %.4f", nm, metrics[[nm]]))
  }
  out
}

#' Run a Frozen-Test Survival ML Workflow
#'
#' @description
#' High-level survival ML workflow for time-to-event prediction. The test set is
#' frozen before feature selection, hyperparameter tuning, final refit, and final
#' evaluation.
#'
#' @param formula Survival formula, for example \code{Surv(time, event) ~ x1 + x2}.
#' @param data Optional full dataset. Required when \code{split} is NULL.
#' @param split Optional \code{ukb_ml_survival_split} object.
#' @param model \code{"cox"}, \code{"rsf"}, \code{"gbm_surv"}, or
#'   \code{"coxnet"}.
#' @param split_params List passed to \code{\link{ukb_ml_survival_split_data}}.
#' @param feature_select \code{"none"}, \code{"filter"}, or \code{"glmnet"}.
#' @param feature_params List passed to
#'   \code{\link{ukb_ml_survival_feature_select}}.
#' @param tune Logical. Run hyperparameter tuning.
#' @param tune_params List passed to \code{\link{ukb_ml_survival_tune}}.
#' @param evaluation_times Time points for survival probability prediction.
#' @param fit_final Logical. Refit final model.
#' @param evaluate_test Logical. Evaluate once on frozen test set.
#' @param seed Optional random seed.
#' @param verbose Logical.
#' @param ... Additional arguments.
#'
#' @return A \code{ukb_ml_survival_workflow} object.
#'
#' @export
ukb_ml_survival_workflow <- function(formula,
                                     data = NULL,
                                     split = NULL,
                                     model = c("cox", "rsf", "gbm_surv", "coxnet"),
                                     split_params = list(),
                                     feature_select = c("none", "filter", "glmnet"),
                                     feature_params = list(),
                                     tune = TRUE,
                                     tune_params = list(),
                                     evaluation_times = c(1, 3, 5, 10),
                                     fit_final = TRUE,
                                     evaluate_test = TRUE,
                                     seed = NULL,
                                     verbose = TRUE,
                                     ...) {
  model <- .sml_check_model(model[[1]])
  feature_select <- match.arg(feature_select)

  data_for_parse <- if (!is.null(split)) split$train else data
  if (is.null(data_for_parse)) {
    stop("Either 'data' or 'split' must be supplied.", call. = FALSE)
  }
  parsed <- .sml_parse_formula(formula, data_for_parse)

  if (is.null(split)) {
    split <- do.call(ukb_ml_survival_split_data, c(
      list(df = data, time = parsed$time_var, event = parsed$event_var, seed = seed, verbose = verbose),
      split_params
    ))
  } else if (!inherits(split, "ukb_ml_survival_split")) {
    stop("'split' must be a ukb_ml_survival_split object.", call. = FALSE)
  }

  feature_result <- do.call(ukb_ml_survival_feature_select, c(
    list(split = split, formula = formula, method = feature_select, seed = seed, verbose = verbose),
    feature_params
  ))

  tune_result <- NULL
  best_params <- list()
  if (isTRUE(tune)) {
    tune_result <- do.call(ukb_ml_survival_tune, c(
      list(split = split, formula = feature_result$formula, model = model, seed = seed, verbose = verbose),
      tune_params
    ))
    best_params <- tune_result$best_params
  }

  final_model <- NULL
  test_eval <- NULL
  if (isTRUE(fit_final)) {
    final_model <- ukb_ml_survival_fit_final(
      split = split,
      formula = feature_result$formula,
      model = model,
      best_params = best_params,
      feature_spec = feature_result,
      seed = seed,
      verbose = verbose
    )
  }
  if (isTRUE(evaluate_test)) {
    if (is.null(final_model)) {
      stop("evaluate_test = TRUE requires fit_final = TRUE.", call. = FALSE)
    }
    test_eval <- ukb_ml_survival_evaluate_test(final_model, split, times = evaluation_times, verbose = verbose)
  }

  out <- list(
    split = split,
    feature_result = feature_result,
    tune_result = tune_result,
    final_model = final_model,
    final_test_metrics = if (!is.null(test_eval)) test_eval$metrics else NULL,
    final_test_predictions = if (!is.null(test_eval)) test_eval$predictions else NULL,
    evaluation_times = evaluation_times,
    model = model,
    formula = formula,
    time_var = parsed$time_var,
    event_var = parsed$event_var,
    seed = seed,
    call = match.call()
  )
  class(out) <- "ukb_ml_survival_workflow"
  out
}

# Variable Importance for Survival

#' Get Variable Importance for Survival Model
#'
#' @param object A ukb_ml_surv object
#' @param ... Additional arguments
#'
#' @return Data frame with variable importance
#'
#' @export
ukb_ml_survival_importance <- function(object, ...) {
  if (inherits(object, "ukb_ml_survival_workflow")) {
    if (is.null(object$final_model)) {
      stop("The workflow does not contain a final model.", call. = FALSE)
    }
    object <- object$final_model
  }
  
  if (inherits(object, "ukb_ml_survival_final")) {
    model <- object$fitted_model
    model_type <- object$model
  } else {
    model <- object$model
    model_type <- object$model_type
  }
  
  importance <- switch(model_type,
    cox = {
      coefs <- stats::coef(model)
      data.frame(
        variable = names(coefs),
        importance = abs(as.numeric(coefs)),
        stringsAsFactors = FALSE
      )
    },
    rsf = {
      imp <- model$importance
      data.frame(
        variable = names(imp),
        importance = as.numeric(imp),
        stringsAsFactors = FALSE
      )
    },
    gbm_surv = {
      .check_ml_package("gbm")
      imp <- summary(model, plotit = FALSE)
      data.frame(
        variable = imp$var,
        importance = imp$rel.inf,
        stringsAsFactors = FALSE
      )
    },
    coxnet = {
      coefs <- as.matrix(coef(model, s = "lambda.min"))
      data.frame(
        variable = rownames(coefs),
        importance = abs(coefs[, 1]),
        stringsAsFactors = FALSE
      )
    }
  )
  
  importance <- importance[order(importance$importance, decreasing = TRUE), ]
  rownames(importance) <- NULL
  
  importance
}

# SHAP for Survival Models

#' SHAP Values for Survival Models
#'
#' @description
#' Compute SHAP values for survival ML models at a specific time point.
#'
#' @param object A ukb_ml_surv object
#' @param data Data for SHAP computation
#' @param time_point Time point for SHAP calculation
#' @param nsim Number of Monte Carlo samples
#' @param sample_n Subsample size
#' @param seed Random seed
#' @param verbose Print progress
#' @param method Model-agnostic SHAP backend.
#'   \code{"permutation"} is retained as an alias for \code{"fastshap"}.
#' @param background_data Optional reference data. Training data are used by
#'   default.
#' @param background_n Maximum number of background observations.
#' @param ... Additional arguments
#'
#' @return A ukb_shap object
#'
#' @export
ukb_ml_survival_shap <- function(object,
                                 data = NULL,
                                 time_point = 5,
                                 nsim = 50,
                                 sample_n = NULL,
                                 seed = NULL,
                                 verbose = TRUE,
                                 method = c("auto", "kernelshap", "fastshap", "permutation", "kernalshap"),
                                 background_data = NULL,
                                 background_n = 200,
                                 ...) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  if (inherits(object, "ukb_ml_survival_workflow")) {
    if (is.null(object$final_model)) {
      stop("The workflow does not contain a final model.", call. = FALSE)
    }
    if (is.null(data)) data <- object$split$test
    if (is.null(background_data)) background_data <- object$split$train
    object <- object$final_model
  }

  if (inherits(object, "ukb_ml_survival_final")) {
    if (is.null(data)) data <- object$train_data
    if (is.null(background_data)) background_data <- object$train_data
    data <- as.data.frame(data)
    features <- object$selected_features %||% object$predictors
    X <- data[, features, drop = FALSE]

    if (!is.null(sample_n) && sample_n < nrow(X)) {
      idx <- sample(nrow(X), sample_n)
      X <- X[idx, , drop = FALSE]
      if (verbose) message(sprintf("Using %d sampled observations for SHAP", sample_n))
    }

    if (verbose) message(sprintf("Computing SHAP values at time = %g...", time_point))

    pred_wrapper <- function(model, newdata) {
      .sml_predict_core(model, newdata = newdata, times = time_point, type = "survival")[, 1]
    }

    background_X <- as.data.frame(background_data)[, features, drop = FALSE]
    background_X <- .ukb_shap_sample_background(background_X, background_n)
    resolved_method <- .ukb_shap_resolve_method(method, object$model)
    backend <- .ukb_shap_model_agnostic(
      method = resolved_method,
      object = object,
      X = X,
      background_X = background_X,
      pred_wrapper = pred_wrapper,
      nsim = nsim,
      verbose = verbose,
      ...
    )

    result <- list(
      shap_values = backend$shap_values,
      baseline = backend$baseline,
      feature_names = colnames(X),
      feature_values = X,
      model_type = object$model,
      task = "survival",
      time_point = time_point,
      method = resolved_method,
      output_scale = "survival_probability",
      background_n = nrow(background_X),
      local_accuracy_error = backend$local_accuracy_error
    )
    class(result) <- "ukb_shap"
    if (verbose) message("SHAP computation complete")
    return(result)
  }
  
  if (is.null(data)) {
    data <- object$test_data
  }
  if (is.null(background_data)) {
    background_data <- object$train_data
  }
  
  X <- data[, object$predictors, drop = FALSE]
  
  if (!is.null(sample_n) && sample_n < nrow(X)) {
    idx <- sample(nrow(X), sample_n)
    X <- X[idx, , drop = FALSE]
    data <- data[idx, , drop = FALSE]
    if (verbose) message(sprintf("Using %d sampled observations for SHAP", sample_n))
  }
  
  if (verbose) message(sprintf("Computing SHAP values at time = %g...", time_point))
  
  pred_wrapper <- function(model, newdata) {
    prediction_object <- object
    prediction_object$model <- model
    ukb_ml_survival_predict(
      prediction_object,
      newdata = newdata,
      times = time_point,
      type = "survival"
    )[, 1]
  }

  background_X <- as.data.frame(background_data)[, object$predictors, drop = FALSE]
  background_X <- .ukb_shap_sample_background(background_X, background_n)
  resolved_method <- .ukb_shap_resolve_method(method, object$model_type)
  backend <- .ukb_shap_model_agnostic(
    method = resolved_method,
    object = object$model,
    X = X,
    background_X = background_X,
    pred_wrapper = pred_wrapper,
    nsim = nsim,
    verbose = verbose,
    ...
  )

  result <- list(
    shap_values = backend$shap_values,
    baseline = backend$baseline,
    feature_names = object$predictors,
    feature_values = X,
    model_type = object$model_type,
    task = "survival",
    time_point = time_point,
    method = resolved_method,
    output_scale = "survival_probability",
    background_n = nrow(background_X),
    local_accuracy_error = backend$local_accuracy_error
  )
  
  class(result) <- "ukb_shap"
  
  if (verbose) message("SHAP computation complete")
  
  result
}

# S3 Methods

#' @export
print.ukb_ml_survival_split <- function(x, ...) {
  cat("\n")
  cat("UKB Survival ML Split\n")
  cat("\n")
  cat(sprintf("Train: %d\n", nrow(x$train)))
  cat(sprintf("Validation: %d\n", if (!is.null(x$validation)) nrow(x$validation) else 0L))
  cat(sprintf("Test: %d\n", nrow(x$test)))
  cat(sprintf("Time variable: %s\n", x$time_var))
  cat(sprintf("Event variable: %s\n", x$event_var))
  invisible(x)
}

#' @export
print.ukb_ml_survival_workflow <- function(x, ...) {
  cat("\n")
  cat("UKB Survival ML Workflow\n")
  cat("\n")
  cat(sprintf("Model: %s\n", x$model))
  cat(sprintf("Time variable: %s\n", x$time_var))
  cat(sprintf("Event variable: %s\n", x$event_var))
  if (!is.null(x$split)) {
    cat(sprintf(
      "Split: train=%d, validation=%d, test=%d\n",
      nrow(x$split$train),
      if (!is.null(x$split$validation)) nrow(x$split$validation) else 0L,
      nrow(x$split$test)
    ))
  }
  if (!is.null(x$final_test_metrics)) {
    cat("\nFinal frozen-test metrics:\n")
    for (nm in names(x$final_test_metrics)) {
      cat(sprintf("  %s: %.3f\n", nm, x$final_test_metrics[[nm]]))
    }
  }
  invisible(x)
}

#' @export
predict.ukb_ml_survival_final <- function(object,
                                          newdata = NULL,
                                          times = c(1, 3, 5),
                                          type = c("survival", "risk", "chf"),
                                          ...) {
  ukb_ml_survival_predict(object, newdata = newdata, times = times, type = type, ...)
}

#' @export
predict.ukb_ml_survival_workflow <- function(object,
                                             newdata = NULL,
                                             times = c(1, 3, 5),
                                             type = c("survival", "risk", "chf"),
                                             ...) {
  ukb_ml_survival_predict(object, newdata = newdata, times = times, type = type, ...)
}

#' @export
print.ukb_ml_surv <- function(x, ...) {
  cat("\n")
  cat("UKB Survival Machine Learning Model\n")
  cat("\n")
  cat(sprintf("Model: %s\n", x$model_type))
  cat(sprintf("Time variable: %s\n", x$time_var))
  cat(sprintf("Event variable: %s\n", x$event_var))
  cat(sprintf("Predictors: %d variables\n", length(x$predictors)))
  cat(sprintf("Train size: %d\n", nrow(x$train_data)))
  cat(sprintf("Test size: %d\n", nrow(x$test_data)))
  cat(sprintf("\nTest C-index: %.3f\n", x$c_index))
  invisible(x)
}

#' @export
summary.ukb_ml_surv <- function(object, ...) {
  print(object)
  
  cat("\nVariable Importance (Top 10):\n")
  imp <- ukb_ml_survival_importance(object)
  print(head(imp, 10))
  
  invisible(object)
}

#' @export
predict.ukb_ml_surv <- function(object, newdata = NULL, times = c(1, 3, 5), type = "survival", ...) {
  ukb_ml_survival_predict(object, newdata = newdata, times = times, type = type, ...)
}
