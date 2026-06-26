#' @keywords internal
.mlw_infer_outcome_type <- function(y) {
  y_nonmissing <- y[!is.na(y)]
  if (length(y_nonmissing) == 0L) {
    stop("Outcome contains only missing values.", call. = FALSE)
  }

  if (is.factor(y_nonmissing) || is.character(y_nonmissing)) {
    n_levels <- length(unique(as.character(y_nonmissing)))
    if (n_levels == 2L) return("binary")
    if (n_levels > 2L) return("multiclass")
  }

  if (is.numeric(y_nonmissing) || is.integer(y_nonmissing) || is.logical(y_nonmissing)) {
    vals <- sort(unique(y_nonmissing))
    if (length(vals) == 2L && all(vals %in% c(0, 1, FALSE, TRUE))) {
      return("binary")
    }
    return("continuous")
  }

  stop("Could not infer outcome_type; set it explicitly.", call. = FALSE)
}

#' @keywords internal
.mlw_resolve_outcome_type <- function(y, outcome_type = c("auto", "binary", "multiclass", "continuous")) {
  outcome_type <- match.arg(outcome_type)
  if (outcome_type == "auto") {
    outcome_type <- .mlw_infer_outcome_type(y)
  }

  y_nonmissing <- y[!is.na(y)]
  if (outcome_type == "binary") {
    n_levels <- length(unique(as.character(y_nonmissing)))
    if (n_levels != 2L) {
      stop("Binary outcome_type requires exactly 2 non-missing outcome values.", call. = FALSE)
    }
  } else if (outcome_type == "multiclass") {
    n_levels <- length(unique(as.character(y_nonmissing)))
    if (n_levels < 3L) {
      stop("Multiclass outcome_type requires at least 3 non-missing outcome values.", call. = FALSE)
    }
  } else if (outcome_type == "continuous" && !(is.numeric(y_nonmissing) || is.integer(y_nonmissing))) {
    stop("Continuous outcome_type requires a numeric outcome.", call. = FALSE)
  }

  outcome_type
}

#' @keywords internal
.mlw_outcome_classes <- function(y) {
  y_nonmissing <- y[!is.na(y)]
  if (is.factor(y_nonmissing)) {
    lv <- levels(droplevels(y_nonmissing))
    return(lv[!is.na(lv) & nzchar(lv)])
  }
  sort(unique(as.character(y_nonmissing)))
}

#' @keywords internal
.mlw_parse_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a model formula.", call. = FALSE)
  }

  response <- all.vars(formula)[1]
  if (!response %in% names(data)) {
    stop(sprintf("Outcome column not found: %s", response), call. = FALSE)
  }

  terms_obj <- terms(formula, data = data)
  predictors <- attr(terms_obj, "term.labels")
  if (length(predictors) == 1L && predictors == ".") {
    predictors <- setdiff(names(data), response)
  }
  if (length(predictors) == 0L) {
    stop("Formula must include at least one predictor.", call. = FALSE)
  }

  missing_predictors <- setdiff(all.vars(stats::delete.response(terms_obj)), names(data))
  if (length(missing_predictors) > 0L) {
    stop("Predictor column(s) not found: ", paste(missing_predictors, collapse = ", "), call. = FALSE)
  }

  list(response = response, predictors = predictors, terms = terms_obj)
}

#' @keywords internal
.mlw_formula_from_features <- function(outcome, features) {
  if (length(features) == 0L) {
    stop("No features available for model fitting.", call. = FALSE)
  }
  stats::as.formula(paste(outcome, "~", paste(features, collapse = " + ")))
}

#' @keywords internal
.mlw_model_frame <- function(formula, data, outcome_type, classes = NULL) {
  data <- as.data.frame(data)
  parsed <- .mlw_parse_formula(formula, data)
  vars <- unique(c(parsed$response, all.vars(stats::delete.response(parsed$terms))))
  complete_idx <- stats::complete.cases(data[, vars, drop = FALSE])
  mf <- data[complete_idx, vars, drop = FALSE]

  if (nrow(mf) == 0L) {
    stop("No complete-case rows available for ML fitting/evaluation.", call. = FALSE)
  }

  y <- mf[[parsed$response]]
  if (outcome_type %in% c("binary", "multiclass")) {
    if (is.null(classes)) {
      classes <- .mlw_outcome_classes(y)
    }
    y <- factor(as.character(y), levels = classes)
    if (any(is.na(y))) {
      stop("Data contains outcome class(es) not present in training classes.", call. = FALSE)
    }
    mf[[parsed$response]] <- y
  } else {
    mf[[parsed$response]] <- as.numeric(y)
  }

  list(
    data = mf,
    outcome = parsed$response,
    predictors = parsed$predictors,
    complete_idx = which(complete_idx),
    classes = classes
  )
}

#' @keywords internal
.mlw_model_matrix <- function(formula, data, outcome_type, classes = NULL) {
  mf <- .mlw_model_frame(formula, data, outcome_type, classes = classes)
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
    y = mf$data[[mf$outcome]],
    model_frame = mf$data,
    outcome = mf$outcome,
    predictors = mf$predictors,
    complete_idx = mf$complete_idx,
    x_terms = x_terms,
    contrasts = contrasts,
    xlevels = stats::.getXlevels(x_terms, mf$data),
    feature_map = feature_map,
    classes = mf$classes
  )
}

#' @keywords internal
.mlw_check_model <- function(model, outcome_type) {
  model <- match.arg(model, c(
    "logistic", "linear", "lm", "rf", "xgboost", "glmnet",
    "svm", "nnet", "rpart", "naive_bayes"
  ))
  if (model == "lm") model <- "linear"

  if (model == "logistic" && outcome_type != "binary") {
    stop("model = 'logistic' is only available for binary outcomes.", call. = FALSE)
  }
  if (model == "linear" && outcome_type != "continuous") {
    stop("model = 'linear' is only available for continuous outcomes.", call. = FALSE)
  }
  if (model == "naive_bayes" && outcome_type == "continuous") {
    stop("model = 'naive_bayes' is only available for binary or multiclass outcomes.", call. = FALSE)
  }

  model
}

#' List Supported Machine Learning Models
#'
#' @description
#' Returns the machine-learning algorithms supported by the UKBAnalytica ML
#' workflow, including eligible outcome types, required R package, and default
#' tuning parameters.
#'
#' @param outcome_type Optional outcome type filter: \code{"all"},
#'   \code{"binary"}, \code{"multiclass"}, or \code{"continuous"}.
#'
#' @return A data.frame describing supported models.
#'
#' @examples
#' ukb_ml_supported_models("binary")
#'
#' @export
ukb_ml_supported_models <- function(outcome_type = c("all", "binary", "multiclass", "continuous")) {
  outcome_type <- match.arg(outcome_type)
  out <- data.frame(
    model = c("logistic", "linear", "rf", "xgboost", "glmnet", "svm", "nnet", "rpart", "naive_bayes"),
    label = c(
      "Logistic regression",
      "Linear regression",
      "Random forest",
      "XGBoost",
      "Regularized regression",
      "Support vector machine",
      "Neural network",
      "Decision tree",
      "Naive Bayes"
    ),
    outcome_types = c(
      "binary",
      "continuous",
      "binary,multiclass,continuous",
      "binary,multiclass,continuous",
      "binary,multiclass,continuous",
      "binary,multiclass,continuous",
      "binary,multiclass,continuous",
      "binary,multiclass,continuous",
      "binary,multiclass"
    ),
    package = c("stats", "stats", "ranger", "xgboost", "glmnet", "e1071", "nnet", "rpart", "e1071"),
    default_tuning = c(
      "none",
      "none",
      "num.trees, mtry, min.node.size",
      "nrounds, max_depth, eta",
      "alpha",
      "cost, gamma",
      "size, decay",
      "cp, maxdepth, minsplit",
      "laplace"
    ),
    stringsAsFactors = FALSE
  )

  if (outcome_type != "all") {
    keep <- vapply(strsplit(out$outcome_types, ","), function(x) outcome_type %in% x, logical(1))
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' @keywords internal
.mlw_default_metric <- function(outcome_type) {
  switch(
    outcome_type,
    binary = "auc",
    multiclass = "macro_f1",
    continuous = "rmse"
  )
}

#' @keywords internal
.mlw_metric_direction <- function(metric) {
  !(metric %in% c("rmse", "mae", "logloss", "brier"))
}

#' @keywords internal
.mlw_auc_binary <- function(truth01, prob) {
  keep <- !is.na(truth01) & !is.na(prob)
  truth01 <- as.integer(truth01[keep])
  prob <- as.numeric(prob[keep])
  if (length(unique(truth01)) != 2L) return(NA_real_)

  n_pos <- sum(truth01 == 1L)
  n_neg <- sum(truth01 == 0L)
  ranks <- rank(prob, ties.method = "average")
  (sum(ranks[truth01 == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' @keywords internal
.mlw_binary_metrics <- function(truth, prob, threshold = 0.5, positive_class = NULL) {
  truth <- droplevels(as.factor(truth))
  if (nlevels(truth) != 2L) {
    stop("Binary metrics require exactly 2 truth classes.", call. = FALSE)
  }

  if (is.null(positive_class)) {
    positive_class <- levels(truth)[2]
  }
  if (!positive_class %in% levels(truth)) {
    stop("'positive_class' must be one of the truth levels.", call. = FALSE)
  }

  truth01 <- as.integer(truth == positive_class)
  pred01 <- as.integer(prob >= threshold)

  tp <- sum(pred01 == 1L & truth01 == 1L, na.rm = TRUE)
  tn <- sum(pred01 == 0L & truth01 == 0L, na.rm = TRUE)
  fp <- sum(pred01 == 1L & truth01 == 0L, na.rm = TRUE)
  fn <- sum(pred01 == 0L & truth01 == 1L, na.rm = TRUE)
  n <- tp + tn + fp + fn

  c(
    auc = .mlw_auc_binary(truth01, prob),
    accuracy = if (n > 0) (tp + tn) / n else NA_real_,
    sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
    specificity = if (tn + fp > 0) tn / (tn + fp) else NA_real_,
    ppv = if (tp + fp > 0) tp / (tp + fp) else NA_real_,
    npv = if (tn + fn > 0) tn / (tn + fn) else NA_real_,
    f1 = if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_,
    brier = mean((prob - truth01)^2, na.rm = TRUE)
  )
}

#' @keywords internal
.mlw_multiclass_metrics <- function(truth, prob) {
  truth <- droplevels(as.factor(truth))
  classes <- levels(truth)
  if (is.null(colnames(prob))) {
    stop("Multiclass probabilities must have class column names.", call. = FALSE)
  }
  prob <- prob[, classes, drop = FALSE]
  pred <- factor(classes[max.col(prob, ties.method = "first")], levels = classes)

  accuracy <- mean(pred == truth, na.rm = TRUE)
  per_class <- lapply(classes, function(cls) {
    tp <- sum(pred == cls & truth == cls, na.rm = TRUE)
    fp <- sum(pred == cls & truth != cls, na.rm = TRUE)
    fn <- sum(pred != cls & truth == cls, na.rm = TRUE)
    support <- sum(truth == cls, na.rm = TRUE)
    precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    recall <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    f1 <- if (!is.na(precision) && !is.na(recall) && precision + recall > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA_real_
    }
    c(precision = precision, recall = recall, f1 = f1, support = support)
  })
  per_class <- do.call(rbind, per_class)

  truth_index <- match(as.character(truth), classes)
  p_true <- prob[cbind(seq_along(truth_index), truth_index)]
  p_true <- pmax(pmin(p_true, 1 - 1e-15), 1e-15)

  c(
    accuracy = accuracy,
    balanced_accuracy = mean(per_class[, "recall"], na.rm = TRUE),
    macro_f1 = mean(per_class[, "f1"], na.rm = TRUE),
    weighted_f1 = stats::weighted.mean(per_class[, "f1"], per_class[, "support"], na.rm = TRUE),
    macro_precision = mean(per_class[, "precision"], na.rm = TRUE),
    macro_recall = mean(per_class[, "recall"], na.rm = TRUE),
    logloss = -mean(log(p_true), na.rm = TRUE)
  )
}

#' @keywords internal
.mlw_continuous_metrics <- function(truth, pred) {
  truth <- as.numeric(truth)
  pred <- as.numeric(pred)
  keep <- !is.na(truth) & !is.na(pred)
  truth <- truth[keep]
  pred <- pred[keep]
  c(
    rmse = sqrt(mean((truth - pred)^2)),
    mae = mean(abs(truth - pred)),
    rsquared = 1 - sum((truth - pred)^2) / sum((truth - mean(truth))^2)
  )
}

#' @keywords internal
.mlw_metrics <- function(truth,
                         prediction,
                         outcome_type,
                         threshold = 0.5,
                         positive_class = NULL) {
  switch(
    outcome_type,
    binary = .mlw_binary_metrics(truth, prediction, threshold, positive_class),
    multiclass = .mlw_multiclass_metrics(truth, prediction),
    continuous = .mlw_continuous_metrics(truth, prediction)
  )
}

#' @keywords internal
.ukb_ml_fit_core <- function(formula,
                             data,
                             model,
                             outcome_type = c("binary", "multiclass", "continuous"),
                             params = list(),
                             classes = NULL,
                             seed = NULL,
                             verbose = FALSE,
                             ...) {
  outcome_type <- match.arg(outcome_type)
  model <- .mlw_check_model(model, outcome_type)
  if (!is.null(seed)) set.seed(seed)

  prep <- .mlw_model_matrix(formula, data, outcome_type, classes = classes)
  mf <- prep$model_frame
  y <- prep$y
  classes <- prep$classes

  fitted <- switch(
    model,
    logistic = {
      stats::glm(formula, data = mf, family = stats::binomial())
    },
    linear = {
      stats::lm(formula, data = mf)
    },
    rf = {
      .check_ml_package("ranger")
      default_params <- list(
        num.trees = 300,
        mtry = NULL,
        min.node.size = if (outcome_type == "continuous") 5 else 1,
        importance = "permutation"
      )
      fit_params <- utils::modifyList(default_params, params)
      if (outcome_type != "continuous") {
        fit_params$probability <- TRUE
      }
      do.call(ranger::ranger, c(list(formula = formula, data = mf), fit_params))
    },
    xgboost = {
      .check_ml_package("xgboost")
      default_params <- list(
        nrounds = 80,
        max_depth = 3,
        eta = 0.1,
        subsample = 1,
        colsample_bytree = 1,
        verbosity = 0
      )
      fit_params <- utils::modifyList(default_params, params)
      nrounds <- fit_params$nrounds
      fit_params$nrounds <- NULL
      if (outcome_type == "binary") {
        label <- as.integer(y == classes[2])
        fit_params$objective <- fit_params$objective %||% "binary:logistic"
        fit_params$eval_metric <- fit_params$eval_metric %||% "auc"
      } else if (outcome_type == "multiclass") {
        label <- as.integer(y) - 1L
        fit_params$objective <- fit_params$objective %||% "multi:softprob"
        fit_params$eval_metric <- fit_params$eval_metric %||% "mlogloss"
        fit_params$num_class <- length(classes)
      } else {
        label <- as.numeric(y)
        fit_params$objective <- fit_params$objective %||% "reg:squarederror"
        fit_params$eval_metric <- fit_params$eval_metric %||% "rmse"
      }
      dtrain <- xgboost::xgb.DMatrix(data = prep$X, label = label)
      xgboost::xgb.train(params = fit_params, data = dtrain, nrounds = nrounds)
    },
    glmnet = {
      .check_ml_package("glmnet")
      default_params <- list(
        alpha = 1,
        family = switch(
          outcome_type,
          binary = "binomial",
          multiclass = "multinomial",
          continuous = "gaussian"
        )
      )
      fit_params <- utils::modifyList(default_params, params)
      do.call(glmnet::cv.glmnet, c(list(x = prep$X, y = y), fit_params))
    },
    svm = {
      .check_ml_package("e1071")
      default_params <- list(kernel = "radial", scale = TRUE)
      fit_params <- utils::modifyList(default_params, params)
      if (outcome_type != "continuous") {
        fit_params$probability <- TRUE
      }
      do.call(e1071::svm, c(list(formula = formula, data = mf), fit_params))
    },
    rpart = {
      .check_ml_package("rpart")
      default_params <- list(cp = 0.01, minsplit = 20, maxdepth = 30)
      fit_params <- utils::modifyList(default_params, params)
      control_args <- fit_params[names(fit_params) %in% c("cp", "minsplit", "maxdepth")]
      fit_params[names(control_args)] <- NULL
      fit_params$control <- fit_params$control %||% do.call(rpart::rpart.control, control_args)
      fit_params$method <- fit_params$method %||% if (outcome_type == "continuous") "anova" else "class"
      do.call(rpart::rpart, c(list(formula = formula, data = mf), fit_params))
    },
    naive_bayes = {
      .check_ml_package("e1071")
      default_params <- list(laplace = 0)
      fit_params <- utils::modifyList(default_params, params)
      do.call(e1071::naiveBayes, c(list(formula = formula, data = mf), fit_params))
    },
    nnet = {
      .check_ml_package("nnet")
      if (outcome_type %in% c("binary", "multiclass")) {
        default_params <- list(trace = FALSE, maxit = 200)
        fit_params <- utils::modifyList(default_params, params)
        do.call(nnet::multinom, c(list(formula = formula, data = mf), fit_params))
      } else {
        default_params <- list(size = 5, decay = 0.01, maxit = 200, trace = FALSE, linout = TRUE)
        fit_params <- utils::modifyList(default_params, params)
        do.call(nnet::nnet, c(list(formula = formula, data = mf), fit_params))
      }
    }
  )

  result <- list(
    fitted_model = fitted,
    model = model,
    outcome_type = outcome_type,
    formula = formula,
    outcome = prep$outcome,
    predictors = prep$predictors,
    classes = classes,
    positive_class = if (outcome_type == "binary") classes[2] else NULL,
    params = params,
    x_terms = prep$x_terms,
    contrasts = prep$contrasts,
    xlevels = prep$xlevels,
    feature_names = colnames(prep$X),
    feature_map = prep$feature_map,
    train_rows = nrow(mf),
    call = match.call()
  )
  class(result) <- c("ukb_ml_core", "ukb_ml_final")
  result
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
.ukb_ml_predict_core <- function(object,
                                 newdata,
                                 type = c("response", "prob", "class")) {
  type <- match.arg(type)
  outcome_type <- object$outcome_type
  model <- object$model

  if (outcome_type == "continuous") {
    pred <- switch(
      model,
      linear = stats::predict(object$fitted_model, newdata = newdata),
      rf = predict(object$fitted_model, data = newdata)$predictions,
      xgboost = {
        X <- stats::model.matrix(object$x_terms, data = newdata, contrasts.arg = object$contrasts)
        if ("(Intercept)" %in% colnames(X)) X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
        predict(object$fitted_model, xgboost::xgb.DMatrix(X))
      },
      glmnet = {
        X <- stats::model.matrix(object$x_terms, data = newdata, contrasts.arg = object$contrasts)
        if ("(Intercept)" %in% colnames(X)) X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
        as.numeric(stats::predict(object$fitted_model, newx = X, s = "lambda.min"))
      },
      svm = as.numeric(stats::predict(object$fitted_model, newdata = newdata)),
      rpart = as.numeric(stats::predict(object$fitted_model, newdata = newdata)),
      nnet = as.numeric(stats::predict(object$fitted_model, newdata = newdata))
    )
    return(as.numeric(pred))
  }

  classes <- object$classes
  prob <- switch(
    model,
    logistic = {
      p <- as.numeric(stats::predict(object$fitted_model, newdata = newdata, type = "response"))
      out <- cbind(1 - p, p)
      colnames(out) <- classes
      out
    },
    rf = {
      pred <- predict(object$fitted_model, data = newdata)$predictions
      pred[, classes, drop = FALSE]
    },
    xgboost = {
      .check_ml_package("xgboost")
      X <- stats::model.matrix(object$x_terms, data = newdata, contrasts.arg = object$contrasts)
      if ("(Intercept)" %in% colnames(X)) X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
      raw <- predict(object$fitted_model, xgboost::xgb.DMatrix(X))
      if (outcome_type == "binary") {
        out <- cbind(1 - raw, raw)
        colnames(out) <- classes
        out
      } else {
        out <- matrix(raw, ncol = length(classes), byrow = TRUE)
        colnames(out) <- classes
        out
      }
    },
    glmnet = {
      X <- stats::model.matrix(object$x_terms, data = newdata, contrasts.arg = object$contrasts)
      if ("(Intercept)" %in% colnames(X)) X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
      raw <- stats::predict(object$fitted_model, newx = X, s = "lambda.min", type = "response")
      if (outcome_type == "binary") {
        p <- as.numeric(raw)
        out <- cbind(1 - p, p)
        colnames(out) <- classes
        out
      } else {
        arr <- raw[, , 1]
        arr[, classes, drop = FALSE]
      }
    },
    svm = {
      raw <- stats::predict(object$fitted_model, newdata = newdata, probability = TRUE)
      probs <- attr(raw, "probabilities")
      probs[, classes, drop = FALSE]
    },
    rpart = {
      probs <- stats::predict(object$fitted_model, newdata = newdata, type = "prob")
      probs[, classes, drop = FALSE]
    },
    naive_bayes = {
      probs <- stats::predict(object$fitted_model, newdata = newdata, type = "raw")
      probs[, classes, drop = FALSE]
    },
    nnet = {
      probs <- stats::predict(object$fitted_model, newdata = newdata, type = "probs")
      if (is.null(dim(probs))) {
        out <- cbind(1 - probs, probs)
        colnames(out) <- classes
        out
      } else {
        probs[, classes, drop = FALSE]
      }
    }
  )

  if (outcome_type == "binary" && type %in% c("response", "prob")) {
    return(as.numeric(prob[, object$positive_class]))
  }
  if (type == "prob") {
    return(prob)
  }
  factor(classes[max.col(prob, ties.method = "first")], levels = classes)
}

#' Standardize Manual ML Train/Test Splits
#'
#' @description
#' Converts user-provided train/test (and optional validation) data frames into a
#' standardized \code{ukb_ml_split} object. This object is consumed by the high
#' level ML workflow and keeps the test set frozen until final evaluation.
#'
#' @param train_data Training/development data.
#' @param test_data Frozen test data.
#' @param validation_data Optional validation data.
#' @param id_col Optional participant ID column used to check overlap.
#' @param check_overlap Logical. Check duplicated and overlapping IDs.
#' @param outcome Outcome column name.
#' @param outcome_type One of \code{"auto"}, \code{"binary"},
#'   \code{"multiclass"}, or \code{"continuous"}.
#'
#' @return A \code{ukb_ml_split} object.
#'
#' @export
ukb_ml_as_split <- function(train_data,
                            test_data,
                            validation_data = NULL,
                            id_col = NULL,
                            check_overlap = TRUE,
                            outcome = NULL,
                            outcome_type = c("auto", "binary", "multiclass", "continuous")) {
  if (!is.data.frame(train_data) || !is.data.frame(test_data)) {
    stop("'train_data' and 'test_data' must be data.frames or data.tables.", call. = FALSE)
  }
  if (nrow(train_data) == 0L || nrow(test_data) == 0L) {
    stop("'train_data' and 'test_data' must both contain rows.", call. = FALSE)
  }
  if (!is.null(validation_data) && (!is.data.frame(validation_data) || nrow(validation_data) == 0L)) {
    stop("'validation_data' must be NULL or a non-empty data.frame/data.table.", call. = FALSE)
  }

  train_data <- data.table::as.data.table(data.table::copy(train_data))
  test_data <- data.table::as.data.table(data.table::copy(test_data))
  validation_data <- if (!is.null(validation_data)) data.table::as.data.table(data.table::copy(validation_data)) else NULL

  if (!is.null(outcome)) {
    for (nm in c("train_data", "test_data", if (!is.null(validation_data)) "validation_data")) {
      obj <- get(nm)
      if (!outcome %in% names(obj)) {
        stop(sprintf("Outcome column '%s' not found in %s.", outcome, nm), call. = FALSE)
      }
    }
    outcome_type <- .mlw_resolve_outcome_type(train_data[[outcome]], outcome_type)

    if (outcome_type %in% c("binary", "multiclass")) {
      train_classes <- .mlw_outcome_classes(train_data[[outcome]])
      check_classes <- function(x, label) {
        classes <- .mlw_outcome_classes(x[[outcome]])
        missing <- setdiff(classes, train_classes)
        if (length(missing) > 0L) {
          stop(sprintf("%s contains outcome class(es) absent from train_data: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
        }
      }
      check_classes(test_data, "test_data")
      if (!is.null(validation_data)) check_classes(validation_data, "validation_data")
    }
  } else {
    outcome_type <- match.arg(outcome_type)
    if (outcome_type == "auto") outcome_type <- NA_character_
  }

  if (isTRUE(check_overlap) && !is.null(id_col)) {
    all_sets <- list(train = train_data, test = test_data)
    if (!is.null(validation_data)) all_sets$validation <- validation_data
    for (nm in names(all_sets)) {
      if (!id_col %in% names(all_sets[[nm]])) {
        stop(sprintf("id_col '%s' not found in %s data.", id_col, nm), call. = FALSE)
      }
      if (any(duplicated(all_sets[[nm]][[id_col]]))) {
        stop(sprintf("Duplicated IDs found in %s data.", nm), call. = FALSE)
      }
    }
    if (length(intersect(train_data[[id_col]], test_data[[id_col]])) > 0L) {
      stop("Overlapping IDs found between train_data and test_data.", call. = FALSE)
    }
    if (!is.null(validation_data)) {
      if (length(intersect(train_data[[id_col]], validation_data[[id_col]])) > 0L ||
          length(intersect(test_data[[id_col]], validation_data[[id_col]])) > 0L) {
        stop("Overlapping IDs found across train/test/validation data.", call. = FALSE)
      }
    }
  }

  out <- list(
    train = train_data,
    validation = validation_data,
    internal_validation = validation_data,
    test = test_data,
    id_col = id_col,
    outcome = outcome,
    outcome_type = outcome_type,
    split_method = "manual",
    split_info = list(
      train_n = nrow(train_data),
      validation_n = if (!is.null(validation_data)) nrow(validation_data) else 0L,
      test_n = nrow(test_data),
      check_overlap = check_overlap
    )
  )
  class(out) <- "ukb_ml_split"
  out
}

#' @keywords internal
.mlw_stratified_sample <- function(idx, strata, prop) {
  groups <- split(idx, strata[idx], drop = TRUE)
  selected <- integer(0)
  for (g in groups) {
    n_g <- length(g)
    n_take <- as.integer(round(n_g * prop))
    if (n_g > 1L) {
      n_take <- max(1L, min(n_g - 1L, n_take))
    } else {
      n_take <- if (prop >= 0.5) 1L else 0L
    }
    if (n_take > 0L) {
      selected <- c(selected, sample(g, n_take))
    }
  }
  sort(unique(selected))
}

#' @keywords internal
.mlw_make_strata <- function(df, outcome, outcome_type, stratify_by, stratify_col, regression_bins) {
  if (identical(stratify_by, "none")) {
    return(rep("all", nrow(df)))
  }
  if (identical(stratify_by, "custom")) {
    if (is.null(stratify_col) || !stratify_col %in% names(df)) {
      stop("stratify_by = 'custom' requires a valid 'stratify_col'.", call. = FALSE)
    }
    strata <- as.character(df[[stratify_col]])
  } else if (identical(stratify_by, "outcome") || identical(stratify_by, "auto")) {
    if (is.null(outcome) || !outcome %in% names(df)) {
      return(rep("all", nrow(df)))
    }
    if (outcome_type == "continuous") {
      probs <- seq(0, 1, length.out = regression_bins + 1L)
      qs <- unique(stats::quantile(df[[outcome]], probs = probs, na.rm = TRUE, type = 7))
      if (length(qs) < 3L) {
        strata <- rep("all", nrow(df))
      } else {
        strata <- as.character(cut(df[[outcome]], breaks = qs, include.lowest = TRUE))
      }
    } else {
      strata <- as.character(df[[outcome]])
    }
  } else {
    stop("Invalid 'stratify_by'.", call. = FALSE)
  }
  strata[is.na(strata) | strata == ""] <- "<NA>"
  strata
}

#' Split Data into Frozen ML Train/Test Sets
#'
#' @description
#' Creates a standardized \code{ukb_ml_split} object for the high-level ML
#' workflow. Supports train/test and train/validation/test splits. The older
#' \code{split_ratio}/\code{stratify_by = <column>} calling style is still
#' accepted for compatibility.
#'
#' @param df A data.frame or data.table.
#' @param outcome Outcome column name. If NULL, a legacy random split is
#'   returned with \code{internal_validation} populated.
#' @param outcome_type One of \code{"auto"}, \code{"binary"},
#'   \code{"multiclass"}, or \code{"continuous"}.
#' @param split Either \code{"train_test"} or \code{"train_valid_test"}.
#' @param train_ratio Training proportion.
#' @param validation_ratio Validation proportion for train/validation/test.
#' @param test_ratio Test proportion.
#' @param split_ratio Deprecated compatibility alias for \code{train_ratio}.
#' @param stratify_by \code{"auto"}, \code{"outcome"}, \code{"custom"},
#'   \code{"none"}, or an older-style column name.
#' @param stratify_col Column used when \code{stratify_by = "custom"}.
#' @param regression_bins Number of quantile bins for continuous outcome
#'   stratification.
#' @param seed Optional random seed.
#' @param verbose Logical. Print split summary.
#'
#' @return A \code{ukb_ml_split} object.
#'
#' @export
ukb_ml_split_data <- function(df,
                              outcome = NULL,
                              outcome_type = c("auto", "binary", "multiclass", "continuous"),
                              split = c("train_test", "train_valid_test"),
                              train_ratio = 0.7,
                              validation_ratio = 0.1,
                              test_ratio = 0.2,
                              split_ratio = NULL,
                              stratify_by = c("auto", "outcome", "custom", "none"),
                              stratify_col = NULL,
                              regression_bins = 5,
                              seed = NULL,
                              verbose = TRUE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame or data.table", call. = FALSE)
  }
  if (nrow(df) < 2L) {
    stop("`df` must contain at least 2 rows", call. = FALSE)
  }
  legacy_split_ratio <- !is.null(split_ratio)
  if (legacy_split_ratio) {
    train_ratio <- split_ratio
    test_ratio <- 1 - split_ratio
    split <- "train_test"
  }

  split <- match.arg(split)
  stratify_by_raw <- stratify_by
  stratify_by <- stratify_by[[1]]
  if (!stratify_by %in% c("auto", "outcome", "custom", "none")) {
    stratify_col <- stratify_by
    stratify_by <- "custom"
  }

  if (!is.null(seed)) set.seed(seed)
  df <- data.table::as.data.table(data.table::copy(df))

  if (!is.null(outcome)) {
    if (!outcome %in% names(df)) {
      stop(sprintf("Outcome column not found: %s", outcome), call. = FALSE)
    }
    outcome_type <- .mlw_resolve_outcome_type(df[[outcome]], outcome_type)
  } else {
    outcome_type <- match.arg(outcome_type)
    if (outcome_type == "auto") outcome_type <- NA_character_
    stratify_by <- if (!is.null(stratify_col)) "custom" else "none"
  }

  n <- nrow(df)
  idx_all <- seq_len(n)
  strata <- .mlw_make_strata(df, outcome, outcome_type, stratify_by, stratify_col, regression_bins)

  if (split == "train_test") {
    if (train_ratio <= 0 || train_ratio >= 1) {
      stop("'train_ratio' must be between 0 and 1.", call. = FALSE)
    }
    train_idx <- if (stratify_by == "none") {
      sort(sample(idx_all, max(1L, min(n - 1L, round(n * train_ratio)))))
    } else {
      .mlw_stratified_sample(idx_all, strata, train_ratio)
    }
    test_idx <- setdiff(idx_all, train_idx)
    validation_idx <- integer(0)
  } else {
    ratio_sum <- train_ratio + validation_ratio + test_ratio
    if (!isTRUE(all.equal(ratio_sum, 1, tolerance = 1e-8))) {
      stop("train_ratio + validation_ratio + test_ratio must equal 1.", call. = FALSE)
    }
    train_idx <- if (stratify_by == "none") {
      sort(sample(idx_all, max(1L, min(n - 2L, round(n * train_ratio)))))
    } else {
      .mlw_stratified_sample(idx_all, strata, train_ratio)
    }
    remaining <- setdiff(idx_all, train_idx)
    remaining_strata <- strata
    valid_prop <- validation_ratio / (validation_ratio + test_ratio)
    validation_idx <- if (stratify_by == "none") {
      sort(sample(remaining, max(1L, min(length(remaining) - 1L, round(length(remaining) * valid_prop)))))
    } else {
      .mlw_stratified_sample(remaining, remaining_strata, valid_prop)
    }
    test_idx <- setdiff(remaining, validation_idx)
  }

  split_obj <- ukb_ml_as_split(
    train_data = df[train_idx],
    validation_data = if (length(validation_idx) > 0L) df[validation_idx] else NULL,
    test_data = df[test_idx],
    outcome = outcome,
    outcome_type = if (is.na(outcome_type)) "auto" else outcome_type,
    check_overlap = FALSE
  )
  if (legacy_split_ratio && is.null(split_obj$validation)) {
    split_obj$internal_validation <- split_obj$test
  }
  split_obj$split_method <- "package"
  split_obj$split_info <- utils::modifyList(split_obj$split_info, list(
    split = split,
    train_ratio = train_ratio,
    validation_ratio = if (split == "train_valid_test") validation_ratio else 0,
    test_ratio = if (split == "train_valid_test") test_ratio else 1 - train_ratio,
    stratify_by = stratify_by,
    stratify_col = stratify_col,
    seed = seed,
    legacy_stratify_by = paste(stratify_by_raw, collapse = ",")
  ))

  if (isTRUE(verbose)) {
    message(sprintf(
      "[ukb_ml_split_data] train=%d, validation=%d, test=%d, outcome_type=%s",
      nrow(split_obj$train),
      if (!is.null(split_obj$validation)) nrow(split_obj$validation) else 0L,
      nrow(split_obj$test),
      ifelse(is.na(split_obj$outcome_type), "unspecified", split_obj$outcome_type)
    ))
  }

  split_obj
}

#' Select Features for UKB ML Workflows
#'
#' @description
#' Performs optional feature selection using only the training/development data
#' in a \code{ukb_ml_split}. The test set is never used for feature selection.
#'
#' @param split A \code{ukb_ml_split} object.
#' @param formula Model formula.
#' @param method \code{"none"}, \code{"boruta"}, \code{"filter"}, or
#'   \code{"glmnet"}.
#' @param outcome_type Outcome type. Defaults to the split outcome type.
#' @param max_features Optional maximum number of selected features.
#' @param boruta_params Parameters passed to \code{Boruta::Boruta()}.
#' @param keep_tentative Logical. Keep Boruta tentative features.
#' @param seed Optional random seed.
#' @param verbose Logical.
#'
#' @return A \code{ukb_ml_feature} object.
#'
#' @export
ukb_ml_feature_select <- function(split,
                                  formula,
                                  method = c("none", "boruta", "filter", "glmnet"),
                                  outcome_type = c("auto", "binary", "multiclass", "continuous"),
                                  max_features = NULL,
                                  boruta_params = list(),
                                  keep_tentative = TRUE,
                                  seed = NULL,
                                  verbose = TRUE) {
  if (!inherits(split, "ukb_ml_split")) {
    stop("'split' must be a ukb_ml_split object.", call. = FALSE)
  }
  method <- match.arg(method)
  parsed <- .mlw_parse_formula(formula, split$train)
  outcome_type <- if (match.arg(outcome_type) == "auto") split$outcome_type else match.arg(outcome_type)
  features <- parsed$predictors

  if (!is.null(seed)) set.seed(seed)

  selected <- features
  status <- data.frame(feature = features, status = "Selected", stringsAsFactors = FALSE)
  info <- list()

  if (method == "boruta") {
    .check_ml_package("Boruta")
    fit_data <- .mlw_model_frame(formula, split$train, outcome_type)$data
    boruta_call <- c(
      list(formula = formula, data = fit_data, doTrace = if (isTRUE(verbose)) 1 else 0),
      boruta_params
    )
    boruta_fit <- do.call(Boruta::Boruta, boruta_call)
    decision <- Boruta::attStats(boruta_fit)
    status <- data.frame(
      feature = rownames(decision),
      status = decision$decision,
      mean_importance = decision$meanImp,
      stringsAsFactors = FALSE
    )
    keep_status <- if (isTRUE(keep_tentative)) c("Confirmed", "Tentative") else "Confirmed"
    selected <- status$feature[status$status %in% keep_status]
    info$boruta <- boruta_fit
  } else if (method == "filter") {
    mf <- .mlw_model_frame(formula, split$train, outcome_type)$data
    y <- mf[[parsed$response]]
    scores <- vapply(features, function(v) {
      x <- mf[[v]]
      if (is.numeric(x) || is.integer(x)) {
        if (outcome_type == "continuous") {
          abs(stats::cor(x, y, use = "complete.obs"))
        } else {
          abs(stats::cor(x, as.numeric(as.factor(y)), use = "complete.obs"))
        }
      } else {
        tbl <- table(x, y)
        suppressWarnings(as.numeric(stats::chisq.test(tbl)$statistic))
      }
    }, numeric(1))
    scores[!is.finite(scores)] <- 0
    status <- data.frame(feature = names(scores), score = as.numeric(scores), status = "Ranked", stringsAsFactors = FALSE)
    status <- status[order(status$score, decreasing = TRUE), , drop = FALSE]
    selected <- status$feature
  } else if (method == "glmnet") {
    .check_ml_package("glmnet")
    fit <- .ukb_ml_fit_core(formula, split$train, model = "glmnet", outcome_type = outcome_type, seed = seed)
    coef_obj <- stats::coef(fit$fitted_model, s = "lambda.min")
    coef_list <- if (is.list(coef_obj)) coef_obj else list(coef_obj)
    nonzero <- unique(unlist(lapply(coef_list, function(z) {
      z <- as.matrix(z)
      rownames(z)[z[, 1] != 0]
    })))
    nonzero <- setdiff(nonzero, "(Intercept)")
    selected <- unique(fit$feature_map$feature[fit$feature_map$matrix_column %in% nonzero])
    if (length(selected) == 0L) {
      selected <- unique(sub("([^:]+).*", "\\1", nonzero))
    }
    selected <- intersect(features, selected)
    status <- data.frame(feature = features, status = ifelse(features %in% selected, "Selected", "Rejected"), stringsAsFactors = FALSE)
    info$glmnet <- fit$fitted_model
  }

  if (!is.null(max_features) && length(selected) > max_features) {
    selected <- selected[seq_len(max_features)]
  }
  if (length(selected) == 0L) {
    stop("Feature selection returned zero features.", call. = FALSE)
  }

  result <- list(
    method = method,
    selected_features = selected,
    selected_status = status,
    formula = .mlw_formula_from_features(parsed$response, selected),
    outcome = parsed$response,
    outcome_type = outcome_type,
    feature_info = info
  )
  class(result) <- "ukb_ml_feature"
  if (isTRUE(verbose)) {
    message(sprintf("[ukb_ml_feature_select] method=%s, selected=%d feature(s)", method, length(selected)))
  }
  result
}

#' @keywords internal
.mlw_param_grid <- function(model, outcome_type, param_grid = NULL, param_space = NULL, search = "grid", n_iter = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(param_grid) && is.null(param_space)) {
    param_grid <- switch(
      model,
      logistic = list(.none = 1),
      linear = list(.none = 1),
      rf = list(num.trees = c(100, 300), mtry = c(2, NA), min.node.size = c(1, 5)),
      xgboost = list(nrounds = c(50, 100), max_depth = c(2, 4), eta = c(0.05, 0.1)),
      glmnet = list(alpha = c(0, 0.5, 1)),
      svm = list(cost = c(1, 2), gamma = c(0.01, 0.1)),
      rpart = list(cp = c(0.001, 0.01), maxdepth = c(5, 10), minsplit = c(10, 20)),
      naive_bayes = list(laplace = c(0, 1)),
      nnet = if (outcome_type == "continuous") {
        list(size = c(3, 5), decay = c(0, 0.01))
      } else {
        list(decay = c(0, 0.01))
      }
    )
  }

  if (search == "grid") {
    grid_source <- if (!is.null(param_grid)) param_grid else param_space
    if (is.data.frame(grid_source)) {
      grid <- grid_source
    } else {
      grid <- do.call(expand.grid, c(grid_source, stringsAsFactors = FALSE))
    }
  } else {
    space <- if (!is.null(param_space)) param_space else param_grid
    if (is.null(space)) stop("'param_space' or 'param_grid' is required for random search.", call. = FALSE)
    n_iter <- if (is.null(n_iter)) 10L else as.integer(n_iter)
    sampled <- lapply(space, function(vals) sample(vals, n_iter, replace = TRUE))
    grid <- as.data.frame(sampled, stringsAsFactors = FALSE)
    grid <- unique(grid)
  }

  if (".none" %in% names(grid)) {
    grid$.none <- NULL
  }
  if (nrow(grid) == 0L) {
    grid <- data.frame(.row = 1L)
    grid$.row <- NULL
  }
  grid
}

#' @keywords internal
.mlw_row_to_params <- function(row_df) {
  params <- as.list(row_df)
  params <- params[!names(params) %in% c(".score", ".fold")]
  params <- lapply(params, function(x) {
    x <- x[[1]]
    if (is.na(x)) NULL else x
  })
  params[!vapply(params, is.null, logical(1))]
}

#' @keywords internal
.mlw_eval_metric <- function(metrics, metric) {
  if (!metric %in% names(metrics)) {
    return(NA_real_)
  }
  as.numeric(metrics[[metric]])
}

#' @keywords internal
.mlw_create_folds <- function(y, k, outcome_type) {
  n <- length(y)
  k <- max(2L, min(as.integer(k), n))
  if (outcome_type %in% c("binary", "multiclass")) {
    y <- as.factor(y)
    class_counts <- table(y)
    if (min(class_counts) < 2L) {
      stop("Cross-validation requires at least 2 rows per outcome class; use a validation split or fewer classes.", call. = FALSE)
    }
    k <- min(k, min(class_counts))
    folds <- vector("list", k)
    for (lv in levels(y)) {
      idx <- which(y == lv)
      idx <- sample(idx)
      fold_ids <- rep(seq_len(k), length.out = length(idx))
      for (i in seq_len(k)) {
        folds[[i]] <- c(folds[[i]], idx[fold_ids == i])
      }
    }
    lapply(folds, sort)
  } else {
    idx <- sample(seq_len(n))
    split(idx, cut(seq_along(idx), k, labels = FALSE))
  }
}

#' Tune ML Hyperparameters Without Touching the Test Set
#'
#' @description
#' Searches model hyperparameters using only the training/development portion of
#' a \code{ukb_ml_split}. The frozen test set is never used.
#'
#' @param split A \code{ukb_ml_split} object.
#' @param formula Model formula.
#' @param model Model type.
#' @param outcome_type Outcome type.
#' @param search \code{"grid"}, \code{"random"}, or \code{"bayes"}. Bayesian
#'   search currently requires \code{rBayesianOptimization}; if unavailable it
#'   falls back to random search with the same parameter space.
#' @param param_grid List or data.frame of candidate parameters.
#' @param param_space Parameter space for random or Bayesian search.
#' @param n_iter Number of random/Bayesian iterations.
#' @param resampling \code{"cv"} or \code{"validation"}.
#' @param folds Number of CV folds.
#' @param metric Metric to optimize.
#' @param maximize Logical. Whether higher metric values are better.
#' @param seed Optional random seed.
#' @param verbose Logical.
#' @param ... Reserved for future extensions.
#'
#' @return A \code{ukb_ml_tune} object.
#'
#' @export
ukb_ml_tune <- function(split,
                        formula,
                        model,
                        outcome_type = c("auto", "binary", "multiclass", "continuous"),
                        search = c("grid", "random", "bayes"),
                        param_grid = NULL,
                        param_space = NULL,
                        n_iter = NULL,
                        resampling = c("cv", "validation"),
                        folds = 5,
                        metric = NULL,
                        maximize = NULL,
                        seed = NULL,
                        verbose = TRUE,
                        ...) {
  if (!inherits(split, "ukb_ml_split")) {
    stop("'split' must be a ukb_ml_split object.", call. = FALSE)
  }
  search <- match.arg(search)
  if (search == "bayes" && !requireNamespace("rBayesianOptimization", quietly = TRUE)) {
    if (isTRUE(verbose)) {
      message("[ukb_ml_tune] rBayesianOptimization unavailable; using random search fallback.")
    }
    search <- "random"
  }
  resampling <- match.arg(resampling)
  if (resampling == "validation" && is.null(split$validation)) {
    if (isTRUE(verbose)) message("[ukb_ml_tune] No validation set; using CV instead.")
    resampling <- "cv"
  }

  parsed <- .mlw_parse_formula(formula, split$train)
  outcome_type <- if (match.arg(outcome_type) == "auto") split$outcome_type else match.arg(outcome_type)
  model <- .mlw_check_model(model, outcome_type)
  metric <- metric %||% .mlw_default_metric(outcome_type)
  maximize <- maximize %||% .mlw_metric_direction(metric)
  grid <- .mlw_param_grid(model, outcome_type, param_grid, param_space, search, n_iter, seed)

  if (!is.null(seed)) set.seed(seed)
  classes <- if (outcome_type %in% c("binary", "multiclass")) {
    .mlw_outcome_classes(split$train[[parsed$response]])
  } else {
    NULL
  }

  results <- vector("list", nrow(grid))
  oof_truth <- NULL
  oof_pred <- NULL
  best_score_seen <- if (isTRUE(maximize)) -Inf else Inf
  best_oof <- NULL

  for (i in seq_len(nrow(grid))) {
    params <- .mlw_row_to_params(grid[i, , drop = FALSE])
    scores <- numeric(0)
    candidate_oof_truth <- NULL
    candidate_oof_pred <- NULL

    if (resampling == "validation") {
      fit <- .ukb_ml_fit_core(formula, split$train, model, outcome_type, params, classes = classes, seed = seed, verbose = FALSE)
      pred <- .ukb_ml_predict_core(fit, split$validation, type = if (outcome_type == "continuous") "response" else "prob")
      truth <- split$validation[[parsed$response]]
      metrics <- .mlw_metrics(truth, pred, outcome_type, threshold = 0.5, positive_class = fit$positive_class)
      scores <- .mlw_eval_metric(metrics, metric)
    } else {
      prep_train <- .mlw_model_frame(formula, split$train, outcome_type, classes = classes)
      fold_idx <- .mlw_create_folds(prep_train$data[[parsed$response]], folds, outcome_type)
      candidate_oof_truth <- prep_train$data[[parsed$response]]
      if (outcome_type == "multiclass") {
        candidate_oof_pred <- matrix(NA_real_, nrow = nrow(prep_train$data), ncol = length(classes))
        colnames(candidate_oof_pred) <- classes
      } else {
        candidate_oof_pred <- rep(NA_real_, nrow(prep_train$data))
      }

      for (fold in seq_along(fold_idx)) {
        valid_idx <- fold_idx[[fold]]
        train_idx <- setdiff(seq_len(nrow(prep_train$data)), valid_idx)
        fit <- .ukb_ml_fit_core(formula, prep_train$data[train_idx, , drop = FALSE], model, outcome_type, params, classes = classes, seed = seed, verbose = FALSE)
        pred <- .ukb_ml_predict_core(fit, prep_train$data[valid_idx, , drop = FALSE], type = if (outcome_type == "continuous") "response" else "prob")
        truth <- prep_train$data[[parsed$response]][valid_idx]
        metrics <- .mlw_metrics(truth, pred, outcome_type, threshold = 0.5, positive_class = fit$positive_class)
        scores <- c(scores, .mlw_eval_metric(metrics, metric))
        if (outcome_type == "multiclass") {
          candidate_oof_pred[valid_idx, ] <- pred
        } else {
          candidate_oof_pred[valid_idx] <- pred
        }
      }
    }

    score <- mean(scores, na.rm = TRUE)
    results[[i]] <- cbind(grid[i, , drop = FALSE], data.frame(.score = score, stringsAsFactors = FALSE))

    is_better <- if (isTRUE(maximize)) score > best_score_seen else score < best_score_seen
    if (isTRUE(is_better)) {
      best_score_seen <- score
      if (!is.null(candidate_oof_truth)) {
        best_oof <- list(truth = candidate_oof_truth, prediction = candidate_oof_pred)
      }
    }
    if (isTRUE(verbose)) {
      message(sprintf("[ukb_ml_tune] %s %d/%d: %s = %.4f", search, i, nrow(grid), metric, score))
    }
  }

  results_df <- do.call(rbind, results)
  best_idx <- if (isTRUE(maximize)) which.max(results_df$.score) else which.min(results_df$.score)
  best_params <- .mlw_row_to_params(results_df[best_idx, , drop = FALSE])

  result <- list(
    model = model,
    outcome_type = outcome_type,
    search = search,
    metric = metric,
    maximize = maximize,
    results = results_df,
    best_params = best_params,
    best_score = results_df$.score[[best_idx]],
    oof = best_oof,
    split_info = split$split_info,
    tuning_info = list(resampling = resampling, folds = folds, n_candidates = nrow(grid))
  )
  class(result) <- "ukb_ml_tune"
  result
}

#' Learn a Binary Classification Threshold
#'
#' @description
#' Selects a binary classification threshold using a fixed value or Youden index
#' on training-development predictions. The test set should never be supplied to
#' this function.
#'
#' @param truth True binary outcome values.
#' @param prob Predicted probability for the positive class.
#' @param method \code{"fixed"} or \code{"youden"}.
#' @param fixed_threshold Threshold used when \code{method = "fixed"}.
#' @param positive_class Optional positive class label.
#'
#' @return A \code{ukb_ml_threshold} object.
#'
#' @export
ukb_ml_threshold <- function(truth,
                             prob,
                             method = c("fixed", "youden"),
                             fixed_threshold = 0.5,
                             positive_class = NULL) {
  method <- match.arg(method)
  truth <- droplevels(as.factor(truth))
  if (nlevels(truth) != 2L) {
    stop("ukb_ml_threshold() only supports binary outcomes.", call. = FALSE)
  }
  if (is.null(positive_class)) {
    positive_class <- levels(truth)[2]
  }
  truth01 <- as.integer(truth == positive_class)
  prob <- as.numeric(prob)
  keep <- !is.na(truth01) & !is.na(prob)
  truth01 <- truth01[keep]
  prob <- prob[keep]

  if (method == "fixed") {
    threshold <- fixed_threshold
  } else {
    candidates <- sort(unique(prob))
    stats <- vapply(candidates, function(th) {
      pred <- as.integer(prob >= th)
      tp <- sum(pred == 1L & truth01 == 1L)
      tn <- sum(pred == 0L & truth01 == 0L)
      fp <- sum(pred == 1L & truth01 == 0L)
      fn <- sum(pred == 0L & truth01 == 1L)
      sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
      spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
      sens + spec - 1
    }, numeric(1))
    threshold <- candidates[which.max(stats)]
  }

  met <- .mlw_binary_metrics(factor(ifelse(truth01 == 1L, positive_class, setdiff(levels(truth), positive_class)[1]), levels = levels(truth)), prob, threshold, positive_class)
  result <- list(
    method = method,
    threshold = as.numeric(threshold),
    sensitivity = unname(met["sensitivity"]),
    specificity = unname(met["specificity"]),
    positive_class = positive_class
  )
  class(result) <- "ukb_ml_threshold"
  result
}

#' Refit the Final ML Model on Training Development Data
#'
#' @description
#' Fits the final model with selected features and tuned parameters using train
#' or train plus validation data. The frozen test set is not used.
#'
#' @param split A \code{ukb_ml_split} object.
#' @param formula Model formula.
#' @param model Model type.
#' @param best_params Best hyperparameters.
#' @param outcome_type Outcome type.
#' @param feature_spec Optional \code{ukb_ml_feature} object.
#' @param threshold Optional \code{ukb_ml_threshold} object.
#' @param use_validation_in_refit Logical. If TRUE, refit on train + validation.
#' @param seed Optional random seed.
#' @param verbose Logical.
#' @param ... Additional arguments.
#'
#' @return A \code{ukb_ml_final} object.
#'
#' @export
ukb_ml_fit_final <- function(split,
                             formula,
                             model,
                             best_params = list(),
                             outcome_type = c("auto", "binary", "multiclass", "continuous"),
                             feature_spec = NULL,
                             threshold = NULL,
                             use_validation_in_refit = TRUE,
                             seed = NULL,
                             verbose = TRUE,
                             ...) {
  if (!inherits(split, "ukb_ml_split")) {
    stop("'split' must be a ukb_ml_split object.", call. = FALSE)
  }
  parsed <- .mlw_parse_formula(formula, split$train)
  outcome_type <- if (match.arg(outcome_type) == "auto") split$outcome_type else match.arg(outcome_type)
  fit_formula <- if (!is.null(feature_spec)) feature_spec$formula else formula

  refit_data <- split$train
  refit_label <- "train"
  if (!is.null(split$validation) && isTRUE(use_validation_in_refit)) {
    refit_data <- data.table::rbindlist(list(split$train, split$validation), use.names = TRUE, fill = TRUE)
    refit_label <- "train_plus_validation"
  }

  classes <- if (outcome_type %in% c("binary", "multiclass")) {
    .mlw_outcome_classes(split$train[[parsed$response]])
  } else {
    NULL
  }
  fit <- .ukb_ml_fit_core(fit_formula, refit_data, model, outcome_type, best_params, classes = classes, seed = seed, verbose = verbose, ...)

  fit$best_params <- best_params
  fit$selected_features <- if (!is.null(feature_spec)) feature_spec$selected_features else .mlw_parse_formula(fit_formula, refit_data)$predictors
  fit$threshold <- threshold
  fit$train_rows <- nrow(split$train)
  fit$validation_rows <- if (!is.null(split$validation)) nrow(split$validation) else 0L
  fit$test_rows <- nrow(split$test)
  fit$refit_data <- refit_label
  fit$split_info <- split$split_info
  class(fit) <- c("ukb_ml_final", "ukb_ml_core")

  if (isTRUE(verbose)) {
    message(sprintf("[ukb_ml_fit_final] model=%s, refit_data=%s, features=%d", fit$model, refit_label, length(fit$selected_features)))
  }
  fit
}

#' Evaluate the Final Model Once on the Frozen Test Set
#'
#' @description
#' Applies the final model, selected features, tuned hyperparameters, and fixed
#' threshold to the frozen test set exactly once.
#'
#' @param object A \code{ukb_ml_final} object.
#' @param split A \code{ukb_ml_split} object.
#' @param metrics Optional metric names to return.
#' @param threshold Optional threshold override for binary classification.
#' @param positive_class Optional positive class label.
#' @param verbose Logical.
#'
#' @return A \code{ukb_ml_test_eval} object.
#'
#' @export
ukb_ml_evaluate_test <- function(object,
                                 split,
                                 metrics = NULL,
                                 threshold = NULL,
                                 positive_class = NULL,
                                 verbose = TRUE) {
  if (!inherits(object, "ukb_ml_final")) {
    stop("'object' must be a ukb_ml_final object.", call. = FALSE)
  }
  if (!inherits(split, "ukb_ml_split")) {
    stop("'split' must be a ukb_ml_split object.", call. = FALSE)
  }

  outcome_type <- object$outcome_type
  truth <- split$test[[object$outcome]]
  pred <- .ukb_ml_predict_core(object, split$test, type = if (outcome_type == "continuous") "response" else "prob")

  if (outcome_type == "binary") {
    th <- threshold %||% if (!is.null(object$threshold)) object$threshold$threshold else 0.5
    positive_class <- positive_class %||% object$positive_class
    all_metrics <- .mlw_metrics(truth, pred, outcome_type, threshold = th, positive_class = positive_class)
  } else {
    th <- NULL
    all_metrics <- .mlw_metrics(truth, pred, outcome_type)
  }

  if (!is.null(metrics)) {
    all_metrics <- all_metrics[intersect(metrics, names(all_metrics))]
  }

  pred_out <- data.frame(
    row_id = seq_len(nrow(split$test)),
    truth = truth,
    stringsAsFactors = FALSE
  )
  if (outcome_type == "binary") {
    pred_out$prob <- pred
    pred_out$pred_class <- factor(ifelse(pred >= th, positive_class, setdiff(object$classes, positive_class)[1]), levels = object$classes)
  } else if (outcome_type == "multiclass") {
    pred_out$pred_class <- factor(colnames(pred)[max.col(pred, ties.method = "first")], levels = object$classes)
    pred_out <- cbind(pred_out, as.data.frame(pred, check.names = FALSE))
  } else {
    pred_out$prediction <- pred
  }

  result <- list(
    metrics = all_metrics,
    predictions = pred_out,
    threshold = th,
    evaluated_at = Sys.time(),
    test_rows = nrow(split$test),
    outcome_type = outcome_type
  )
  class(result) <- "ukb_ml_test_eval"

  if (isTRUE(verbose)) {
    message("[ukb_ml_evaluate_test] Final frozen-test metrics:")
    for (nm in names(all_metrics)) {
      message(sprintf("  %s: %.4f", nm, all_metrics[[nm]]))
    }
  }
  result
}

.mlw_default_xgboost_grid <- function(y,
                                      positive_class = NULL,
                                      nthread = max(1L, min(4L, parallel::detectCores(logical = FALSE)))) {
  classes <- .mlw_outcome_classes(y)
  if (length(classes) != 2L) {
    stop("Default XGBoost grid with scale_pos_weight requires a binary outcome.", call. = FALSE)
  }
  positive_class <- positive_class %||% classes[[2]]
  if (!positive_class %in% classes) {
    stop("`positive_class` must be one of: ", paste(classes, collapse = ", "), call. = FALSE)
  }
  y_chr <- as.character(y)
  pos <- sum(y_chr == positive_class, na.rm = TRUE)
  neg <- sum(!is.na(y_chr) & y_chr != positive_class)
  scale_pos_weight <- if (pos > 0L) neg / pos else 1

  expand.grid(
    nrounds = c(100, 200),
    max_depth = c(2, 3),
    eta = c(0.03, 0.10),
    subsample = 0.80,
    colsample_bytree = c(0.60, 0.80),
    min_child_weight = 1,
    gamma = 0,
    lambda = 1,
    alpha = 0,
    scale_pos_weight = scale_pos_weight,
    nthread = nthread,
    stringsAsFactors = FALSE
  )
}

.mlw_feature_set_grid <- function(param_grid, model_id) {
  if (is.null(param_grid)) {
    return(NULL)
  }
  if (!is.data.frame(param_grid) && is.list(param_grid) && model_id %in% names(param_grid)) {
    return(param_grid[[model_id]])
  }
  param_grid
}

#' Run a Complete Single-Model UKB ML Flow
#'
#' @description
#' High-level single-model interface for common UK Biobank machine-learning
#' analyses. The function can create or consume a frozen train/test split, tune
#' model hyperparameters, learn a binary threshold, fit the final model, evaluate
#' the frozen test set, prepare ROC data, and optionally compute SHAP values.
#'
#' @param formula Model formula. Required unless both \code{outcome} and
#'   \code{features} are supplied.
#' @param data Optional full dataset. Used to create a split when \code{split} is
#'   \code{NULL} and \code{train_data}/\code{test_data} are not supplied.
#' @param split Optional \code{ukb_ml_split} object.
#' @param train_data,test_data,validation_data Optional pre-split datasets used
#'   when \code{split} is \code{NULL}.
#' @param id_col Optional participant ID column for overlap checks and output
#'   predictions.
#' @param outcome Optional outcome column. Defaults to the response in
#'   \code{formula} or \code{split$outcome}.
#' @param features Optional feature names. Used when \code{formula} is
#'   \code{NULL}.
#' @param model Model type passed to \code{\link{ukb_ml_tune}} and
#'   \code{\link{ukb_ml_fit_final}}.
#' @param model_id,model_label Optional model identifier and display label.
#' @param outcome_type Outcome type.
#' @param split_params List passed to \code{\link{ukb_ml_split_data}} when
#'   splitting a full \code{data} object.
#' @param param_grid Optional hyperparameter grid.
#' @param tune Logical. Run hyperparameter tuning.
#' @param tune_params Additional arguments passed to \code{\link{ukb_ml_tune}}.
#' @param best_params Optional final model parameters when \code{tune = FALSE}.
#' @param threshold_method \code{"none"}, \code{"fixed"}, or \code{"youden"}.
#' @param threshold_params Additional arguments passed to
#'   \code{\link{ukb_ml_threshold}}.
#' @param metrics Optional metric names passed to \code{\link{ukb_ml_evaluate_test}}.
#' @param positive_class Optional positive class label for binary outcomes.
#' @param use_validation_in_refit Logical passed to \code{\link{ukb_ml_fit_final}}.
#' @param compute_shap Logical. Compute SHAP values for the final model.
#' @param shap_data Optional data used for SHAP. Defaults to the frozen test set.
#' @param shap_params Additional arguments passed to \code{\link{ukb_shap}}.
#' @param seed Optional random seed.
#' @param verbose Logical.
#'
#' @return A \code{ukb_ml_flow} object with standardized components:
#'   \code{split}, \code{formula}, \code{features}, \code{tune},
#'   \code{threshold}, \code{final_model}, \code{test_eval},
#'   \code{metrics}, \code{predictions}, \code{roc}, and optional \code{shap}.
#'
#' @export
ukb_ml_flow <- function(formula = NULL,
                        data = NULL,
                        split = NULL,
                        train_data = NULL,
                        test_data = NULL,
                        validation_data = NULL,
                        id_col = NULL,
                        outcome = NULL,
                        features = NULL,
                        model = "xgboost",
                        model_id = "model",
                        model_label = NULL,
                        outcome_type = c("auto", "binary", "multiclass", "continuous"),
                        split_params = list(),
                        param_grid = NULL,
                        tune = TRUE,
                        tune_params = list(),
                        best_params = NULL,
                        threshold_method = c("none", "fixed", "youden"),
                        threshold_params = list(),
                        metrics = NULL,
                        positive_class = NULL,
                        use_validation_in_refit = FALSE,
                        compute_shap = FALSE,
                        shap_data = NULL,
                        shap_params = list(),
                        seed = NULL,
                        verbose = TRUE) {
  outcome_type_arg <- match.arg(outcome_type)
  threshold_method <- match.arg(threshold_method)
  model_label <- model_label %||% model_id

  if (!is.null(formula)) {
    if (!inherits(formula, "formula")) {
      stop("`formula` must be a model formula.", call. = FALSE)
    }
    response <- all.vars(formula)[1]
    outcome <- outcome %||% response
  } else {
    if (is.null(outcome) || is.null(features)) {
      stop("Supply `formula`, or supply both `outcome` and `features`.", call. = FALSE)
    }
    formula <- .mlw_formula_from_features(outcome, features)
  }

  if (is.null(split)) {
    if (!is.null(train_data) || !is.null(test_data)) {
      if (is.null(train_data) || is.null(test_data)) {
        stop("Both `train_data` and `test_data` are required when using pre-split data.", call. = FALSE)
      }
      split <- ukb_ml_as_split(
        train_data = train_data,
        test_data = test_data,
        validation_data = validation_data,
        id_col = id_col,
        check_overlap = !is.null(id_col),
        outcome = outcome,
        outcome_type = outcome_type_arg
      )
    } else if (!is.null(data)) {
      split_call <- utils::modifyList(
        list(
          df = data,
          outcome = outcome,
          outcome_type = outcome_type_arg,
          seed = seed,
          verbose = verbose
        ),
        split_params
      )
      split <- do.call(ukb_ml_split_data, split_call)
    } else {
      stop("Supply `split`, `train_data` + `test_data`, or full `data`.", call. = FALSE)
    }
  } else if (!inherits(split, "ukb_ml_split")) {
    stop("`split` must be a ukb_ml_split object.", call. = FALSE)
  }

  outcome <- outcome %||% split$outcome
  if (is.null(outcome)) {
    stop("Could not determine outcome. Supply `outcome` or use a split with `split$outcome`.", call. = FALSE)
  }
  outcome_type_final <- if (outcome_type_arg == "auto") {
    split$outcome_type %||% .mlw_resolve_outcome_type(split$train[[outcome]], "auto")
  } else {
    outcome_type_arg
  }

  parsed <- .mlw_parse_formula(formula, split$train)
  features <- parsed$predictors
  available_features <- intersect(features, intersect(names(split$train), names(split$test)))
  if (length(available_features) != length(features)) {
    missing_features <- setdiff(features, available_features)
    stop("Feature column(s) missing from train/test split: ", paste(missing_features, collapse = ", "), call. = FALSE)
  }

  positive_class <- if (outcome_type_final == "binary") {
    classes <- .mlw_outcome_classes(split$train[[outcome]])
    as.character(positive_class %||% classes[[2]])
  } else {
    NULL
  }
  if (outcome_type_final != "binary" && threshold_method != "none") {
    stop("Threshold learning is only available for binary outcomes.", call. = FALSE)
  }

  grid <- param_grid
  if (is.null(grid) && identical(model, "xgboost") && outcome_type_final == "binary") {
    grid <- .mlw_default_xgboost_grid(split$train[[outcome]], positive_class = positive_class)
  }

  tune_res <- NULL
  if (isTRUE(tune)) {
    tune_call <- utils::modifyList(
      list(
        split = split,
        formula = formula,
        model = model,
        outcome_type = outcome_type_final,
        search = "grid",
        param_grid = grid,
        resampling = "cv",
        folds = 10,
        metric = if (outcome_type_final == "continuous") "rmse" else "auc",
        maximize = if (outcome_type_final == "continuous") FALSE else TRUE,
        seed = seed,
        verbose = verbose
      ),
      tune_params
    )
    tune_res <- do.call(ukb_ml_tune, tune_call)
    best_params <- tune_res$best_params
  } else {
    best_params <- best_params %||% list()
  }

  threshold_res <- NULL
  if (outcome_type_final == "binary" && threshold_method != "none") {
    threshold_source <- NULL
    if (!is.null(tune_res) && !is.null(tune_res$oof)) {
      threshold_truth <- tune_res$oof$truth
      threshold_prob <- as.numeric(tune_res$oof$prediction)
      threshold_source <- "training_oof"
    } else {
      temp_fit <- ukb_ml_fit_final(
        split = split,
        formula = formula,
        model = model,
        best_params = best_params,
        outcome_type = outcome_type_final,
        use_validation_in_refit = FALSE,
        seed = seed,
        verbose = FALSE
      )
      if (!is.null(split$validation)) {
        threshold_truth <- split$validation[[outcome]]
        threshold_prob <- as.numeric(.ukb_ml_predict_core(temp_fit, split$validation, type = "prob"))
        threshold_source <- "validation"
      } else {
        threshold_truth <- split$train[[outcome]]
        threshold_prob <- as.numeric(.ukb_ml_predict_core(temp_fit, split$train, type = "prob"))
        threshold_source <- "train"
      }
    }
    threshold_call <- utils::modifyList(
      list(
        truth = threshold_truth,
        prob = threshold_prob,
        method = threshold_method,
        positive_class = positive_class
      ),
      threshold_params
    )
    threshold_res <- do.call(ukb_ml_threshold, threshold_call)
    threshold_res$source <- threshold_source
  }

  final_model <- ukb_ml_fit_final(
    split = split,
    formula = formula,
    model = model,
    best_params = best_params,
    outcome_type = outcome_type_final,
    threshold = threshold_res,
    use_validation_in_refit = use_validation_in_refit,
    seed = seed,
    verbose = verbose
  )

  test_eval <- ukb_ml_evaluate_test(
    object = final_model,
    split = split,
    metrics = metrics,
    threshold = if (!is.null(threshold_res)) threshold_res$threshold else NULL,
    positive_class = positive_class,
    verbose = verbose
  )

  predictions <- test_eval$predictions
  predictions$model_id <- model_id
  predictions$model_label <- model_label
  if (!is.null(split$id_col) && split$id_col %in% names(split$test)) {
    predictions[[split$id_col]] <- split$test[[split$id_col]]
  }
  first_cols <- c("model_id", "model_label", split$id_col)
  predictions <- predictions[, c(intersect(first_cols, names(predictions)), setdiff(names(predictions), first_cols)), drop = FALSE]

  metrics_df <- data.frame(
    model_id = model_id,
    model_label = model_label,
    metric = names(test_eval$metrics),
    value = as.numeric(test_eval$metrics),
    stringsAsFactors = FALSE
  )
  if (!is.null(threshold_res)) {
    metrics_df <- rbind(
      metrics_df,
      data.frame(
        model_id = model_id,
        model_label = model_label,
        metric = paste0("threshold_", threshold_res$method),
        value = threshold_res$threshold,
        stringsAsFactors = FALSE
      )
    )
  }

  roc_df <- NULL
  if (outcome_type_final == "binary" && "prob" %in% names(predictions)) {
    roc_df <- ukb_ml_roc_data(
      truth = predictions$truth,
      prob = predictions$prob,
      model_id = model_id,
      model_label = model_label,
      positive_class = positive_class
    )
  }

  shap_res <- NULL
  if (isTRUE(compute_shap)) {
    shap_data <- shap_data %||% split$test
    shap_call <- utils::modifyList(
      list(
        object = final_model,
        data = shap_data,
        seed = seed,
        verbose = verbose
      ),
      shap_params
    )
    shap_res <- do.call(ukb_shap, shap_call)
  }

  out <- list(
    split = split,
    formula = formula,
    outcome = outcome,
    features = features,
    model_id = model_id,
    model_label = model_label,
    model = model,
    outcome_type = outcome_type_final,
    tune = tune_res,
    threshold = threshold_res,
    final_model = final_model,
    test_eval = test_eval,
    metrics = metrics_df,
    predictions = predictions,
    roc = roc_df,
    shap = shap_res,
    call = match.call()
  )
  class(out) <- "ukb_ml_flow"
  out
}

#' Plot a UKB ML Flow Object
#'
#' @param x A \code{ukb_ml_flow} object.
#' @param type Plot type: \code{"roc"} or \code{"shap_beeswarm"}.
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return A ggplot2 object.
#'
#' @export
plot.ukb_ml_flow <- function(x, type = c("roc", "shap_beeswarm"), ...) {
  type <- match.arg(type)
  if (!inherits(x, "ukb_ml_flow")) {
    stop("`x` must be a ukb_ml_flow object.", call. = FALSE)
  }
  if (type == "roc") {
    if (is.null(x$roc)) {
      stop("This flow object does not contain ROC data.", call. = FALSE)
    }
    return(plot_ml_roc_compare(x$roc, ...))
  }
  if (is.null(x$shap)) {
    stop("This flow object does not contain SHAP values. Run ukb_ml_flow(..., compute_shap = TRUE).", call. = FALSE)
  }
  plot_shap_beeswarm(x$shap, ...)
}

.mlw_named_arg <- function(x, primary, secondary = NULL, default = NULL) {
  if (is.null(x)) {
    return(default)
  }
  if (!is.null(names(x))) {
    if (!is.null(primary) && primary %in% names(x)) {
      return(x[[primary]])
    }
    if (!is.null(secondary) && secondary %in% names(x)) {
      return(x[[secondary]])
    }
  }
  x
}

.mlw_select_compare_arg <- function(x, combo_id, model, feature_set_id) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x) && !is.null(names(x))) {
    if (combo_id %in% names(x)) return(x[[combo_id]])
    if (model %in% names(x)) return(x[[model]])
    if (feature_set_id %in% names(x)) return(x[[feature_set_id]])
  }
  x
}

#' Compare Multiple Feature Sets and/or Models
#'
#' @description
#' Batch-runs \code{\link{ukb_ml_flow}} across feature-set and model
#' combinations. The same frozen train/test split is reused for every
#' combination, making the output suitable for comparing different feature
#' groups, different machine-learning algorithms, or the full
#' feature-set-by-model grid.
#'
#' @param formula Optional base formula. The response is used as the outcome.
#'   Predictors are used as the default feature set when \code{feature_sets} is
#'   \code{NULL}.
#' @param data,split,train_data,test_data,validation_data,id_col Passed to
#'   \code{\link{ukb_ml_flow}}.
#' @param outcome Optional outcome column. Required when \code{formula} is
#'   \code{NULL}.
#' @param feature_sets Optional named list of feature vectors. If \code{NULL},
#'   one feature set is derived from \code{formula} or \code{features}.
#' @param features Optional feature names used when \code{formula} and
#'   \code{feature_sets} are \code{NULL}.
#' @param models Character vector of models supported by
#'   \code{\link{ukb_ml_supported_models}}.
#' @param compare Comparison mode: \code{"auto"}, \code{"feature_sets"},
#'   \code{"models"}, or \code{"both"}. In \code{"auto"} mode, all supplied
#'   feature-set and model combinations are evaluated.
#' @param outcome_type Outcome type passed to \code{\link{ukb_ml_flow}}.
#' @param feature_set_labels Optional labels for feature sets.
#' @param model_labels Optional labels for models.
#' @param param_grid Optional hyperparameter grid. Can be a single grid shared
#'   by all combinations, a named list keyed by model, feature set, or
#'   \code{"feature_set__model"}.
#' @param tune_params Optional list passed to \code{\link{ukb_ml_tune}}. Can also
#'   be keyed by model, feature set, or combination.
#' @param threshold_params Optional list passed to \code{\link{ukb_ml_threshold}}.
#'   Can also be keyed by model, feature set, or combination.
#' @param ... Additional arguments passed to \code{\link{ukb_ml_flow}}, including
#'   \code{outcome_type}, \code{split_params}, \code{threshold_method},
#'   \code{metrics}, \code{positive_class}, \code{use_validation_in_refit},
#'   \code{compute_shap}, \code{shap_params}, \code{seed}, and \code{verbose}.
#'
#' @return A \code{ukb_ml_flow_compare} object containing \code{flows},
#'   \code{metrics}, \code{comparison}, \code{predictions}, \code{roc}, and
#'   \code{thresholds}.
#'
#' @export
ukb_ml_compare_flows <- function(formula = NULL,
                                 data = NULL,
                                 split = NULL,
                                 train_data = NULL,
                                 test_data = NULL,
                                 validation_data = NULL,
                                 id_col = NULL,
                                 outcome = NULL,
                                 feature_sets = NULL,
                                 features = NULL,
                                 models = "xgboost",
                                 compare = c("auto", "feature_sets", "models", "both"),
                                 outcome_type = c("auto", "binary", "multiclass", "continuous"),
                                 feature_set_labels = NULL,
                                 model_labels = NULL,
                                 param_grid = NULL,
                                 tune_params = list(),
                                 threshold_params = list(),
                                 ...) {
  compare <- match.arg(compare)
  if (!is.character(models) || length(models) == 0L || anyNA(models)) {
    stop("`models` must be a non-empty character vector.", call. = FALSE)
  }
  outcome_type_arg <- match.arg(outcome_type)

  if (is.null(outcome) && !is.null(formula)) {
    ref_data <- split$train %||% train_data %||% data
    if (!is.null(ref_data)) {
      parsed_formula <- .mlw_parse_formula(formula, ref_data)
      outcome <- parsed_formula$response
    }
  }
  if (is.null(outcome) && !is.null(split)) {
    outcome <- split$outcome
  }

  outcome_type_for_model_check <- outcome_type_arg
  if (outcome_type_for_model_check == "auto") {
    y_ref <- NULL
    if (!is.null(split) && !is.null(outcome) && outcome %in% names(split$train)) {
      y_ref <- split$train[[outcome]]
    } else if (!is.null(train_data) && !is.null(outcome) && outcome %in% names(train_data)) {
      y_ref <- train_data[[outcome]]
    } else if (!is.null(data) && !is.null(outcome) && outcome %in% names(data)) {
      y_ref <- data[[outcome]]
    }
    outcome_type_for_model_check <- if (!is.null(y_ref)) .mlw_resolve_outcome_type(y_ref, "auto") else "binary"
  }
  models <- vapply(models, .mlw_check_model, character(1), outcome_type = outcome_type_for_model_check)

  if (is.null(feature_sets)) {
    if (!is.null(formula)) {
      ref_data <- split$train %||% train_data %||% data
      if (is.null(ref_data)) {
        stop("A data source is required to parse `formula` when `feature_sets` is NULL.", call. = FALSE)
      }
      parsed <- .mlw_parse_formula(formula, ref_data)
      outcome <- outcome %||% parsed$response
      feature_sets <- list(features = parsed$predictors)
    } else {
      if (is.null(outcome) || is.null(features)) {
        stop("Supply `feature_sets`, or supply `formula`, or supply both `outcome` and `features`.", call. = FALSE)
      }
      feature_sets <- list(features = features)
    }
  }
  if (!is.list(feature_sets) || length(feature_sets) == 0L) {
    stop("`feature_sets` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(feature_sets)) || any(!nzchar(names(feature_sets)))) {
    names(feature_sets) <- paste0("feature_set_", seq_along(feature_sets))
  }

  if (compare == "feature_sets" && length(models) != 1L) {
    stop("compare = 'feature_sets' requires exactly one model.", call. = FALSE)
  }
  if (compare == "models" && length(feature_sets) != 1L) {
    stop("compare = 'models' requires exactly one feature set.", call. = FALSE)
  }

  if (is.null(feature_set_labels)) {
    feature_set_labels <- stats::setNames(names(feature_sets), names(feature_sets))
  } else if (is.null(names(feature_set_labels))) {
    feature_set_labels <- stats::setNames(as.character(feature_set_labels), names(feature_sets))
  }
  if (is.null(model_labels)) {
    model_meta <- ukb_ml_supported_models()
    model_labels <- stats::setNames(model_meta$label, model_meta$model)
  } else if (is.null(names(model_labels))) {
    model_labels <- stats::setNames(as.character(model_labels), models)
  }

  flows <- list()
  for (feature_set_id in names(feature_sets)) {
    feature_set_label <- unname(feature_set_labels[[feature_set_id]] %||% feature_set_id)
    for (model in models) {
      combo_id <- paste(feature_set_id, model, sep = "__")
      model_label <- unname(model_labels[[model]] %||% model)
      combo_label <- if (length(feature_sets) > 1L && length(models) > 1L) {
        paste(feature_set_label, model_label, sep = " - ")
      } else if (length(feature_sets) > 1L) {
        feature_set_label
      } else {
        model_label
      }

      flow <- ukb_ml_flow(
        formula = .mlw_formula_from_features(outcome, feature_sets[[feature_set_id]]),
        data = data,
        split = split,
        train_data = train_data,
        test_data = test_data,
        validation_data = validation_data,
        id_col = id_col,
        outcome = outcome,
        model = model,
        model_id = combo_id,
        model_label = combo_label,
        outcome_type = outcome_type_arg,
        param_grid = .mlw_select_compare_arg(param_grid, combo_id, model, feature_set_id),
        tune_params = .mlw_select_compare_arg(tune_params, combo_id, model, feature_set_id) %||% list(),
        threshold_params = .mlw_select_compare_arg(threshold_params, combo_id, model, feature_set_id) %||% list(),
        ...
      )
      flow$feature_set_id <- feature_set_id
      flow$feature_set_label <- feature_set_label
      flows[[combo_id]] <- flow
    }
  }

  add_meta <- function(df, flow) {
    if (is.null(df)) {
      return(NULL)
    }
    df$feature_set_id <- flow$feature_set_id
    df$feature_set_label <- flow$feature_set_label
    df$model <- flow$model
    df
  }

  metrics_all <- do.call(rbind, lapply(flows, function(x) add_meta(x$metrics, x)))
  predictions_all <- do.call(rbind, lapply(flows, function(x) add_meta(x$predictions, x)))
  roc_all <- do.call(rbind, lapply(flows, function(x) add_meta(x$roc, x)))

  thresholds_all <- do.call(rbind, lapply(flows, function(x) {
    if (is.null(x$threshold)) {
      return(NULL)
    }
    out <- data.frame(
      model_id = x$model_id,
      model_label = x$model_label,
      feature_set_id = x$feature_set_id,
      feature_set_label = x$feature_set_label,
      model = x$model,
      method = x$threshold$method,
      source = x$threshold$source,
      threshold = x$threshold$threshold,
      sensitivity = x$threshold$sensitivity,
      specificity = x$threshold$specificity,
      positive_class = x$threshold$positive_class,
      stringsAsFactors = FALSE
    )
    out
  }))

  comparison <- stats::reshape(
    metrics_all,
    idvar = c("model_id", "model_label", "feature_set_id", "feature_set_label", "model"),
    timevar = "metric",
    direction = "wide"
  )
  names(comparison) <- sub("^value\\.", "", names(comparison))
  if ("auc" %in% names(comparison)) {
    comparison <- comparison[order(comparison$auc, decreasing = TRUE), , drop = FALSE]
  }
  rownames(comparison) <- NULL

  out <- list(
    flows = flows,
    metrics = metrics_all,
    comparison = comparison,
    predictions = predictions_all,
    roc = roc_all,
    thresholds = thresholds_all,
    design = expand.grid(
      feature_set_id = names(feature_sets),
      model = models,
      stringsAsFactors = FALSE
    ),
    call = match.call()
  )
  class(out) <- "ukb_ml_flow_compare"
  out
}

#' Plot a UKB ML Flow Comparison Object
#'
#' @param x A \code{ukb_ml_flow_compare} object.
#' @param type Plot type. Currently \code{"roc"}.
#' @param ... Additional arguments passed to \code{\link{plot_ml_roc_compare}}.
#'
#' @return A ggplot2 object.
#'
#' @export
plot.ukb_ml_flow_compare <- function(x, type = c("roc"), ...) {
  type <- match.arg(type)
  if (!inherits(x, "ukb_ml_flow_compare")) {
    stop("`x` must be a ukb_ml_flow_compare object.", call. = FALSE)
  }
  if (type == "roc") {
    if (is.null(x$roc)) {
      stop("This comparison object does not contain ROC data.", call. = FALSE)
    }
    return(plot_ml_roc_compare(x$roc, ...))
  }
}

#' Compare Multiple Feature Sets with a Frozen-Test ML Workflow
#'
#' @description
#' Runs the same machine-learning workflow across multiple feature sets using a
#' shared \code{ukb_ml_split}. For binary outcomes, the function can tune models
#' by cross-validation, learn a threshold on training-development predictions,
#' refit the final model, evaluate the frozen test set, and return unified
#' metrics, prediction, threshold, and ROC tables.
#'
#' @param split A \code{ukb_ml_split} object.
#' @param feature_sets Named list of character vectors. Each vector contains the
#'   feature names used by one model.
#' @param outcome Optional outcome column. Defaults to \code{split$outcome}.
#' @param model Model type passed to \code{\link{ukb_ml_tune}} and
#'   \code{\link{ukb_ml_fit_final}}.
#' @param outcome_type Outcome type. Currently this helper is intended for
#'   binary classification.
#' @param model_labels Optional labels for feature sets. Can be a named vector or
#'   a vector in the same order as \code{feature_sets}.
#' @param param_grid Optional parameter grid. Can be a single grid shared by all
#'   models or a named list keyed by feature-set name.
#' @param tune_params Additional arguments passed to \code{\link{ukb_ml_tune}}.
#' @param threshold_method \code{"none"}, \code{"fixed"}, or \code{"youden"}.
#' @param threshold_params Additional arguments passed to
#'   \code{\link{ukb_ml_threshold}}.
#' @param metrics Optional metric names passed to \code{\link{ukb_ml_evaluate_test}}.
#' @param positive_class Optional positive class label for binary outcomes.
#' @param use_validation_in_refit Logical passed to \code{\link{ukb_ml_fit_final}}.
#' @param seed Optional random seed.
#' @param verbose Logical.
#'
#' @return A \code{ukb_ml_feature_set_compare} object containing per-feature-set
#'   models and unified result tables.
#'
#' @export
ukb_ml_compare_feature_sets <- function(split,
                                        feature_sets,
                                        outcome = NULL,
                                        model = "xgboost",
                                        outcome_type = c("auto", "binary"),
                                        model_labels = NULL,
                                        param_grid = NULL,
                                        tune_params = list(),
                                        threshold_method = c("none", "fixed", "youden"),
                                        threshold_params = list(),
                                        metrics = c("auc", "accuracy", "sensitivity", "specificity", "ppv", "npv", "f1", "brier"),
                                        positive_class = NULL,
                                        use_validation_in_refit = FALSE,
                                        seed = NULL,
                                        verbose = TRUE) {
  if (!inherits(split, "ukb_ml_split")) {
    stop("`split` must be a ukb_ml_split object.", call. = FALSE)
  }
  if (!is.list(feature_sets) || length(feature_sets) == 0L) {
    stop("`feature_sets` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(feature_sets)) || any(!nzchar(names(feature_sets)))) {
    names(feature_sets) <- paste0("model_", seq_along(feature_sets))
  }

  outcome <- outcome %||% split$outcome
  if (is.null(outcome) || !outcome %in% names(split$train)) {
    stop("`outcome` must be supplied or stored in `split$outcome`.", call. = FALSE)
  }
  outcome_type <- match.arg(outcome_type)
  outcome_type <- if (outcome_type == "auto") split$outcome_type %||% .mlw_resolve_outcome_type(split$train[[outcome]], "auto") else outcome_type
  if (!identical(outcome_type, "binary")) {
    stop("ukb_ml_compare_feature_sets() currently supports binary outcomes.", call. = FALSE)
  }

  classes <- .mlw_outcome_classes(split$train[[outcome]])
  positive_class <- positive_class %||% classes[[2]]
  threshold_method <- match.arg(threshold_method)

  if (is.null(model_labels)) {
    model_labels <- names(feature_sets)
  } else if (is.null(names(model_labels))) {
    model_labels <- stats::setNames(as.character(model_labels), names(feature_sets))
  }

  model_results <- vector("list", length(feature_sets))
  names(model_results) <- names(feature_sets)

  for (model_id in names(feature_sets)) {
    model_label <- unname(model_labels[[model_id]] %||% model_id)
    features <- unique(as.character(feature_sets[[model_id]]))
    features <- intersect(features, intersect(names(split$train), names(split$test)))
    if (length(features) == 0L) {
      stop("No available features for feature set: ", model_id, call. = FALSE)
    }
    if (isTRUE(verbose)) {
      message(sprintf("[ukb_ml_compare_feature_sets] %s: %d feature(s)", model_label, length(features)))
    }

    ml_formula <- .mlw_formula_from_features(outcome, features)
    grid_i <- .mlw_feature_set_grid(param_grid, model_id)
    if (is.null(grid_i) && identical(model, "xgboost")) {
      grid_i <- .mlw_default_xgboost_grid(split$train[[outcome]], positive_class = positive_class)
    }

    tune_call <- utils::modifyList(
      list(
        split = split,
        formula = ml_formula,
        model = model,
        outcome_type = outcome_type,
        search = "grid",
        param_grid = grid_i,
        resampling = "cv",
        folds = 10,
        metric = "auc",
        maximize = TRUE,
        seed = seed,
        verbose = verbose
      ),
      tune_params
    )
    tune_res <- do.call(ukb_ml_tune, tune_call)

    threshold_res <- NULL
    if (threshold_method != "none") {
      threshold_source <- NULL
      if (!is.null(tune_res$oof)) {
        threshold_truth <- tune_res$oof$truth
        threshold_prob <- as.numeric(tune_res$oof$prediction)
        threshold_source <- "training_oof"
      } else {
        temp_fit <- ukb_ml_fit_final(
          split = split,
          formula = ml_formula,
          model = model,
          best_params = tune_res$best_params,
          outcome_type = outcome_type,
          use_validation_in_refit = FALSE,
          seed = seed,
          verbose = FALSE
        )
        if (!is.null(split$validation)) {
          threshold_truth <- split$validation[[outcome]]
          threshold_prob <- as.numeric(.ukb_ml_predict_core(temp_fit, split$validation, type = "prob"))
          threshold_source <- "validation"
        } else {
          threshold_truth <- split$train[[outcome]]
          threshold_prob <- as.numeric(.ukb_ml_predict_core(temp_fit, split$train, type = "prob"))
          threshold_source <- "train"
        }
      }
      threshold_call <- utils::modifyList(
        list(
          truth = threshold_truth,
          prob = threshold_prob,
          method = threshold_method,
          positive_class = positive_class
        ),
        threshold_params
      )
      threshold_res <- do.call(ukb_ml_threshold, threshold_call)
      threshold_res$source <- threshold_source
    }

    final_model <- ukb_ml_fit_final(
      split = split,
      formula = ml_formula,
      model = model,
      best_params = tune_res$best_params,
      outcome_type = outcome_type,
      threshold = threshold_res,
      use_validation_in_refit = use_validation_in_refit,
      seed = seed,
      verbose = verbose
    )

    eval_res <- ukb_ml_evaluate_test(
      object = final_model,
      split = split,
      metrics = metrics,
      threshold = if (!is.null(threshold_res)) threshold_res$threshold else NULL,
      positive_class = positive_class,
      verbose = verbose
    )

    pred <- eval_res$predictions
    pred$model_id <- model_id
    pred$model_label <- model_label
    if (!is.null(split$id_col) && split$id_col %in% names(split$test)) {
      pred[[split$id_col]] <- split$test[[split$id_col]]
    }
    first_cols <- c("model_id", "model_label", split$id_col)
    pred <- pred[, c(intersect(first_cols, names(pred)), setdiff(names(pred), first_cols)), drop = FALSE]

    metrics_df <- data.frame(
      model_id = model_id,
      model_label = model_label,
      metric = names(eval_res$metrics),
      value = as.numeric(eval_res$metrics),
      stringsAsFactors = FALSE
    )
    if (!is.null(threshold_res)) {
      metrics_df <- rbind(
        metrics_df,
        data.frame(
          model_id = model_id,
          model_label = model_label,
          metric = paste0("threshold_", threshold_res$method),
          value = threshold_res$threshold,
          stringsAsFactors = FALSE
        )
      )
    }

    roc_df <- ukb_ml_roc_data(
      truth = pred$truth,
      prob = pred$prob,
      model_id = model_id,
      model_label = model_label,
      positive_class = positive_class
    )

    threshold_df <- if (!is.null(threshold_res)) {
      data.frame(
        model_id = model_id,
        model_label = model_label,
        method = threshold_res$method,
        source = threshold_res$source,
        threshold = threshold_res$threshold,
        sensitivity = threshold_res$sensitivity,
        specificity = threshold_res$specificity,
        positive_class = threshold_res$positive_class,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }

    model_results[[model_id]] <- list(
      model_id = model_id,
      model_label = model_label,
      features = features,
      formula = ml_formula,
      tune = tune_res,
      threshold = threshold_res,
      final_model = final_model,
      test_eval = eval_res,
      predictions = pred,
      metrics = metrics_df,
      roc = roc_df,
      threshold_table = threshold_df
    )
  }

  metrics_all <- do.call(rbind, lapply(model_results, `[[`, "metrics"))
  predictions_all <- do.call(rbind, lapply(model_results, `[[`, "predictions"))
  roc_all <- do.call(rbind, lapply(model_results, `[[`, "roc"))
  threshold_all <- do.call(rbind, lapply(model_results, `[[`, "threshold_table"))

  comparison <- stats::reshape(
    metrics_all,
    idvar = c("model_id", "model_label"),
    timevar = "metric",
    direction = "wide"
  )
  names(comparison) <- sub("^value\\.", "", names(comparison))
  if ("auc" %in% names(comparison)) {
    comparison <- comparison[order(comparison$auc, decreasing = TRUE), , drop = FALSE]
  }
  rownames(comparison) <- NULL

  out <- list(
    models = model_results,
    metrics = metrics_all,
    comparison = comparison,
    predictions = predictions_all,
    roc = roc_all,
    thresholds = threshold_all,
    call = match.call()
  )
  class(out) <- "ukb_ml_feature_set_compare"
  out
}

#' Run a Frozen-Test UKB ML Workflow
#'
#' @description
#' High-level, publication-oriented ML workflow for binary, multiclass, and
#' continuous outcomes. The test set is frozen before feature selection,
#' hyperparameter tuning, threshold learning, and final refit.
#'
#' @param formula Model formula.
#' @param data Optional full dataset. Required when \code{split} is NULL.
#' @param split Optional \code{ukb_ml_split} object.
#' @param model Model type: \code{"logistic"}, \code{"linear"}, \code{"rf"},
#'   \code{"xgboost"}, \code{"glmnet"}, \code{"svm"}, \code{"nnet"},
#'   \code{"rpart"}, or \code{"naive_bayes"}.
#' @param outcome_type \code{"auto"}, \code{"binary"}, \code{"multiclass"}, or
#'   \code{"continuous"}.
#' @param split_params List passed to \code{\link{ukb_ml_split_data}} when
#'   \code{split} is NULL.
#' @param feature_select \code{"none"}, \code{"boruta"}, \code{"filter"}, or
#'   \code{"glmnet"}.
#' @param feature_params List passed to \code{\link{ukb_ml_feature_select}}.
#' @param tune Logical. Run hyperparameter tuning.
#' @param tune_params List passed to \code{\link{ukb_ml_tune}}.
#' @param threshold_method \code{"none"}, \code{"fixed"}, or \code{"youden"}.
#' @param threshold_params List passed to \code{\link{ukb_ml_threshold}}.
#' @param fit_final Logical. Refit final model.
#' @param evaluate_test Logical. Evaluate once on frozen test set.
#' @param seed Optional random seed.
#' @param verbose Logical.
#' @param ... Additional arguments.
#'
#' @return A \code{ukb_ml_workflow} object.
#'
#' @export
ukb_ml_workflow <- function(formula,
                            data = NULL,
                            split = NULL,
                            model,
                            outcome_type = c("auto", "binary", "multiclass", "continuous"),
                            split_params = list(),
                            feature_select = c("none", "boruta", "filter", "glmnet"),
                            feature_params = list(),
                            tune = TRUE,
                            tune_params = list(),
                            threshold_method = c("none", "fixed", "youden"),
                            threshold_params = list(),
                            fit_final = TRUE,
                            evaluate_test = TRUE,
                            seed = NULL,
                            verbose = TRUE,
                            ...) {
  if (is.null(split) && is.null(data)) {
    stop("Either 'data' or 'split' must be supplied.", call. = FALSE)
  }
  if (!is.null(split) && !inherits(split, "ukb_ml_split")) {
    stop("'split' must be a ukb_ml_split object.", call. = FALSE)
  }

  parsed_data <- if (!is.null(split)) split$train else data
  parsed <- .mlw_parse_formula(formula, parsed_data)
  outcome_type_arg <- match.arg(outcome_type)

  if (is.null(split)) {
    split_call <- c(
      list(df = data, outcome = parsed$response, outcome_type = outcome_type_arg, seed = seed, verbose = verbose),
      split_params
    )
    split <- do.call(ukb_ml_split_data, split_call)
  }

  outcome_type_final <- if (outcome_type_arg == "auto") split$outcome_type else outcome_type_arg
  model <- .mlw_check_model(model, outcome_type_final)

  feature_select <- match.arg(feature_select)
  feature_result <- do.call(
    ukb_ml_feature_select,
    c(
      list(split = split, formula = formula, method = feature_select, outcome_type = outcome_type_final, seed = seed, verbose = verbose),
      feature_params
    )
  )
  fit_formula <- feature_result$formula

  if (isTRUE(tune)) {
    tune_result <- do.call(
      ukb_ml_tune,
      c(
        list(split = split, formula = fit_formula, model = model, outcome_type = outcome_type_final, seed = seed, verbose = verbose),
        tune_params
      )
    )
    best_params <- tune_result$best_params
  } else {
    tune_result <- NULL
    best_params <- tune_params$best_params %||% list()
  }

  threshold_method <- match.arg(threshold_method)
  threshold_result <- NULL
  if (outcome_type_final == "binary" && threshold_method != "none") {
    if (!is.null(split$validation)) {
      temp_fit <- .ukb_ml_fit_core(fit_formula, split$train, model, outcome_type_final, best_params, seed = seed, verbose = FALSE)
      prob <- .ukb_ml_predict_core(temp_fit, split$validation, type = "prob")
      truth <- split$validation[[temp_fit$outcome]]
      threshold_result <- do.call(
        ukb_ml_threshold,
        c(list(truth = truth, prob = prob, method = threshold_method, positive_class = temp_fit$positive_class), threshold_params)
      )
      threshold_result$source <- "validation"
    } else if (!is.null(tune_result$oof)) {
      threshold_result <- do.call(
        ukb_ml_threshold,
        c(list(truth = tune_result$oof$truth, prob = tune_result$oof$prediction, method = threshold_method), threshold_params)
      )
      threshold_result$source <- "oof"
    } else {
      temp_fit <- .ukb_ml_fit_core(fit_formula, split$train, model, outcome_type_final, best_params, seed = seed, verbose = FALSE)
      prob <- .ukb_ml_predict_core(temp_fit, split$train, type = "prob")
      threshold_result <- do.call(
        ukb_ml_threshold,
        c(list(truth = split$train[[temp_fit$outcome]], prob = prob, method = threshold_method, positive_class = temp_fit$positive_class), threshold_params)
      )
      threshold_result$source <- "train"
    }
    if (isTRUE(verbose)) {
      message(sprintf("[ukb_ml_workflow] threshold=%.4f (%s)", threshold_result$threshold, threshold_result$source))
    }
  } else if (outcome_type_final != "binary" && threshold_method != "none") {
    stop("Threshold learning is only available for binary outcomes.", call. = FALSE)
  }

  final_model <- NULL
  final_eval <- NULL
  if (isTRUE(fit_final)) {
    final_model <- ukb_ml_fit_final(
      split = split,
      formula = fit_formula,
      model = model,
      best_params = best_params,
      outcome_type = outcome_type_final,
      feature_spec = feature_result,
      threshold = threshold_result,
      seed = seed,
      verbose = verbose
    )
  }

  if (isTRUE(evaluate_test)) {
    if (is.null(final_model)) {
      stop("evaluate_test = TRUE requires fit_final = TRUE.", call. = FALSE)
    }
    final_eval <- ukb_ml_evaluate_test(final_model, split, verbose = verbose)
  }

  result <- list(
    split = split,
    feature_result = feature_result,
    tune_result = tune_result,
    threshold_result = threshold_result,
    final_model = final_model,
    final_test_metrics = if (!is.null(final_eval)) final_eval$metrics else NULL,
    final_test_predictions = if (!is.null(final_eval)) final_eval$predictions else NULL,
    workflow_info = list(
      model = model,
      outcome_type = outcome_type_final,
      test_set_used_only_for_final_eval = isTRUE(evaluate_test),
      seed = seed,
      created_at = Sys.time()
    ),
    call = match.call()
  )
  class(result) <- "ukb_ml_workflow"
  result
}

#' @export
print.ukb_ml_split <- function(x, ...) {
  cat("\nUKB ML Split\n")
  cat(sprintf("Train: %d\n", nrow(x$train)))
  cat(sprintf("Validation: %d\n", if (!is.null(x$validation)) nrow(x$validation) else 0L))
  cat(sprintf("Test: %d\n", nrow(x$test)))
  cat(sprintf("Outcome: %s\n", x$outcome %||% "<unspecified>"))
  cat(sprintf("Outcome type: %s\n", x$outcome_type %||% "<unspecified>"))
  invisible(x)
}

#' @export
print.ukb_ml_flow <- function(x, ...) {
  cat("\nUKB ML Flow\n")
  cat(sprintf("Model: %s\n", x$model %||% "<unspecified>"))
  cat(sprintf("Label: %s\n", x$model_label %||% x$model_id %||% "<unspecified>"))
  cat(sprintf("Outcome: %s\n", x$outcome %||% "<unspecified>"))
  cat(sprintf("Outcome type: %s\n", x$outcome_type %||% "<unspecified>"))
  cat(sprintf("Features: %d\n", length(x$features)))
  if (!is.null(x$tune)) {
    cat(sprintf("Best %s: %.4f\n", x$tune$metric, x$tune$best_score))
  }
  if (!is.null(x$threshold)) {
    cat(sprintf("Threshold: %.4f (%s)\n", x$threshold$threshold, x$threshold$source %||% x$threshold$method))
  }
  if (!is.null(x$test_eval) && length(x$test_eval$metrics) > 0L) {
    cat("Frozen test metrics:\n")
    for (nm in names(x$test_eval$metrics)) {
      cat(sprintf("  %s: %.4f\n", nm, x$test_eval$metrics[[nm]]))
    }
  }
  invisible(x)
}

#' @export
print.ukb_ml_workflow <- function(x, ...) {
  cat("\nUKB ML Workflow\n")
  cat(sprintf("Model: %s\n", x$workflow_info$model))
  cat(sprintf("Outcome type: %s\n", x$workflow_info$outcome_type))
  if (!is.null(x$feature_result)) {
    cat(sprintf("Features: %d selected by %s\n", length(x$feature_result$selected_features), x$feature_result$method))
  }
  if (!is.null(x$tune_result)) {
    cat(sprintf("Best %s: %.4f\n", x$tune_result$metric, x$tune_result$best_score))
  }
  if (!is.null(x$threshold_result)) {
    cat(sprintf("Threshold: %.4f (%s)\n", x$threshold_result$threshold, x$threshold_result$source %||% x$threshold_result$method))
  }
  if (!is.null(x$final_test_metrics)) {
    cat("Final test metrics:\n")
    for (nm in names(x$final_test_metrics)) {
      cat(sprintf("  %s: %.4f\n", nm, x$final_test_metrics[[nm]]))
    }
  }
  invisible(x)
}

#' @export
predict.ukb_ml_final <- function(object, newdata, type = c("response", "prob", "class"), ...) {
  .ukb_ml_predict_core(object, newdata, type = type)
}
