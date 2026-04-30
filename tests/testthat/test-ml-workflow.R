if (!exists("ukb_ml_workflow", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "ml_model.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ml_workflow.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ml_shap.R"), local = FALSE)
}

test_that("ukb_ml_split_data creates frozen splits and keeps legacy alias", {
  set.seed(1)
  df <- data.frame(
    eid = seq_len(80),
    age = rnorm(80, 60, 7),
    bmi = rnorm(80, 27, 4),
    case = factor(rep(c("control", "case"), each = 40), levels = c("control", "case"))
  )

  split_obj <- ukb_ml_split_data(
    df,
    outcome = "case",
    split = "train_valid_test",
    train_ratio = 0.6,
    validation_ratio = 0.2,
    test_ratio = 0.2,
    seed = 42,
    verbose = FALSE
  )

  expect_s3_class(split_obj, "ukb_ml_split")
  expect_equal(split_obj$outcome_type, "binary")
  expect_equal(nrow(split_obj$train) + nrow(split_obj$validation) + nrow(split_obj$test), nrow(df))
  expect_false(any(split_obj$test$eid %in% split_obj$train$eid))

  legacy <- ukb_ml_split_data(df, split_ratio = 0.75, stratify_by = "case", seed = 42, verbose = FALSE)
  expect_s3_class(legacy, "ukb_ml_split")
  expect_equal(nrow(legacy$internal_validation), nrow(legacy$test))
})

test_that("ukb_ml_as_split validates manual split overlap", {
  df <- data.frame(
    eid = seq_len(12),
    x = rnorm(12),
    y = factor(rep(c("control", "case"), each = 6), levels = c("control", "case"))
  )

  expect_error(
    ukb_ml_as_split(
      train_data = df[1:8, ],
      test_data = df[8:12, ],
      id_col = "eid",
      outcome = "y"
    ),
    "Overlapping IDs"
  )
})

test_that("ukb_ml_workflow runs binary classification with threshold learning", {
  set.seed(2)
  n <- 140
  df <- data.frame(
    age = rnorm(n, 60, 8),
    bmi = rnorm(n, 27, 4),
    sex = factor(sample(c("F", "M"), n, TRUE))
  )
  eta <- -7 + 0.07 * df$age + 0.11 * df$bmi + 0.35 * (df$sex == "M")
  df$case <- factor(ifelse(runif(n) < stats::plogis(eta), "case", "control"), levels = c("control", "case"))

  wf <- ukb_ml_workflow(
    case ~ .,
    data = df,
    model = "logistic",
    split_params = list(split = "train_valid_test", train_ratio = 0.65, validation_ratio = 0.15, test_ratio = 0.20),
    feature_select = "filter",
    feature_params = list(max_features = 3),
    tune = TRUE,
    tune_params = list(resampling = "validation", metric = "auc"),
    threshold_method = "youden",
    seed = 123,
    verbose = FALSE
  )

  expect_s3_class(wf, "ukb_ml_workflow")
  expect_true(is.finite(wf$final_test_metrics[["auc"]]))
  expect_true(wf$threshold_result$threshold >= 0 && wf$threshold_result$threshold <= 1)
  expect_equal(nrow(wf$final_test_predictions), nrow(wf$split$test))
})

test_that("ukb_ml_workflow runs continuous regression", {
  set.seed(3)
  n <- 120
  df <- data.frame(
    age = rnorm(n, 58, 7),
    bmi = rnorm(n, 26, 4),
    sex = factor(sample(c("F", "M"), n, TRUE))
  )
  df$y <- 0.4 * df$age + 1.2 * df$bmi + ifelse(df$sex == "M", 3, 0) + rnorm(n, 0, 4)

  wf <- ukb_ml_workflow(
    y ~ age + bmi + sex,
    data = df,
    model = "linear",
    outcome_type = "continuous",
    split_params = list(split = "train_test", train_ratio = 0.75),
    feature_select = "none",
    tune = TRUE,
    tune_params = list(folds = 3, metric = "rmse"),
    seed = 99,
    verbose = FALSE
  )

  expect_s3_class(wf, "ukb_ml_workflow")
  expect_true(is.finite(wf$final_test_metrics[["rmse"]]))
  expect_equal(nrow(wf$final_test_predictions), nrow(wf$split$test))
})

test_that("glmnet feature selection maps dummy variables back to source features", {
  testthat::skip_if_not_installed("glmnet")

  set.seed(5)
  df <- data.frame(
    x = rnorm(80),
    sex = factor(rep(c("F", "M"), 40))
  )
  df$y <- factor(
    ifelse(df$x + 1.2 * (df$sex == "M") + rnorm(80) > 0, "case", "control"),
    levels = c("control", "case")
  )
  split_obj <- ukb_ml_split_data(df, outcome = "y", split = "train_test", seed = 1, verbose = FALSE)

  selected <- ukb_ml_feature_select(split_obj, y ~ x + sex, method = "glmnet", verbose = FALSE)

  expect_true(all(selected$selected_features %in% c("x", "sex")))
  expect_true("sex" %in% selected$selected_features)
})

test_that("new lightweight workflow models run only when selected", {
  set.seed(6)
  n <- 120
  df <- data.frame(
    age = rnorm(n, 60, 8),
    bmi = rnorm(n, 27, 4),
    sex = factor(sample(c("F", "M"), n, TRUE))
  )
  eta <- -7 + 0.08 * df$age + 0.10 * df$bmi + 0.4 * (df$sex == "M")
  df$case <- factor(ifelse(runif(n) < stats::plogis(eta), "case", "control"), levels = c("control", "case"))

  tree <- ukb_ml_workflow(
    case ~ age + bmi + sex,
    data = df,
    model = "rpart",
    split_params = list(split = "train_test", train_ratio = 0.75),
    tune = TRUE,
    tune_params = list(resampling = "cv", folds = 3, param_grid = data.frame(cp = c(0.001, 0.01), maxdepth = c(3, 5), minsplit = c(10, 10))),
    threshold_method = "youden",
    seed = 42,
    verbose = FALSE
  )
  expect_s3_class(tree, "ukb_ml_workflow")
  expect_true(is.finite(tree$final_test_metrics[["auc"]]))

  nb <- ukb_ml_workflow(
    case ~ age + bmi + sex,
    data = df,
    model = "naive_bayes",
    split_params = list(split = "train_valid_test"),
    tune = TRUE,
    tune_params = list(resampling = "validation", param_grid = data.frame(laplace = c(0, 1))),
    threshold_method = "youden",
    seed = 42,
    verbose = FALSE
  )
  expect_s3_class(nb, "ukb_ml_workflow")
  expect_true(is.finite(nb$final_test_metrics[["auc"]]))
})

test_that("ukb_shap supports ukb_ml_workflow objects", {
  testthat::skip_if_not_installed("fastshap")

  set.seed(7)
  n <- 60
  df <- data.frame(
    age = rnorm(n, 60, 8),
    bmi = rnorm(n, 27, 4)
  )
  df$case <- factor(
    ifelse(df$age + df$bmi + rnorm(n, 0, 5) > 88, "case", "control"),
    levels = c("control", "case")
  )

  wf <- ukb_ml_workflow(
    case ~ age + bmi,
    data = df,
    model = "logistic",
    split_params = list(split = "train_test", train_ratio = 0.75),
    tune = FALSE,
    threshold_method = "fixed",
    seed = 1,
    verbose = FALSE
  )

  shap <- ukb_shap(wf, nsim = 2, sample_n = 5, seed = 1, verbose = FALSE)

  expect_s3_class(shap, "ukb_shap")
  expect_equal(colnames(shap$shap_values), c("age", "bmi"))
  expect_equal(nrow(shap$shap_values), 5)
})

test_that("ukb_ml_workflow runs multiclass classification with manual split", {
  testthat::skip_if_not_installed("nnet")

  set.seed(4)
  n_each <- 36
  df <- data.frame(
    eid = seq_len(n_each * 3),
    x1 = c(rnorm(n_each, 1), rnorm(n_each, -1), rnorm(n_each, 0)),
    x2 = c(rnorm(n_each, 0), rnorm(n_each, 1), rnorm(n_each, -1)),
    x3 = factor(sample(c("low", "high"), n_each * 3, TRUE)),
    outcome = factor(rep(c("A", "B", "C"), each = n_each))
  )

  train_idx <- unlist(lapply(split(seq_len(nrow(df)), df$outcome), function(idx) idx[1:24]))
  valid_idx <- unlist(lapply(split(seq_len(nrow(df)), df$outcome), function(idx) idx[25:30]))
  test_idx <- unlist(lapply(split(seq_len(nrow(df)), df$outcome), function(idx) idx[31:36]))
  split_obj <- ukb_ml_as_split(
    train_data = df[train_idx, ],
    validation_data = df[valid_idx, ],
    test_data = df[test_idx, ],
    id_col = "eid",
    outcome = "outcome",
    outcome_type = "multiclass"
  )

  wf <- ukb_ml_workflow(
    outcome ~ x1 + x2 + x3,
    split = split_obj,
    model = "nnet",
    outcome_type = "multiclass",
    feature_select = "none",
    tune = TRUE,
    tune_params = list(resampling = "validation", param_grid = data.frame(decay = c(0, 0.01))),
    threshold_method = "none",
    seed = 101,
    verbose = FALSE
  )

  expect_s3_class(wf, "ukb_ml_workflow")
  expect_true(is.finite(wf$final_test_metrics[["macro_f1"]]))
  expect_equal(nrow(wf$final_test_predictions), length(test_idx))
})

test_that("ukb_ml_threshold is binary-only", {
  expect_error(
    ukb_ml_threshold(
      truth = factor(c("A", "B", "C")),
      prob = c(0.1, 0.2, 0.3),
      method = "youden"
    ),
    "binary"
  )
})
