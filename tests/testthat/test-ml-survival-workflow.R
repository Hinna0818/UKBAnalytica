if (!exists("ukb_ml_survival_workflow", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "ml_model.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ml_workflow.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ml_survival.R"), local = FALSE)
}

simulate_survival_ml_data <- function(n = 160, seed = 1) {
  set.seed(seed)
  df <- data.frame(
    eid = seq_len(n),
    age = rnorm(n, 60, 8),
    bmi = rnorm(n, 27, 4),
    sex = factor(sample(c("F", "M"), n, TRUE))
  )
  lp <- 0.04 * (df$age - 60) + 0.08 * (df$bmi - 27) + 0.3 * (df$sex == "M")
  event_time <- stats::rexp(n, rate = 0.08 * exp(lp))
  censor_time <- stats::rexp(n, rate = 0.04)
  df$time <- pmax(0.01, pmin(event_time, censor_time))
  df$event <- as.integer(event_time <= censor_time)
  df
}

test_that("ukb_ml_survival_split_data creates frozen survival splits", {
  df <- simulate_survival_ml_data()

  split <- ukb_ml_survival_split_data(
    df,
    time = "time",
    event = "event",
    split = "train_valid_test",
    train_ratio = 0.6,
    validation_ratio = 0.2,
    test_ratio = 0.2,
    seed = 42,
    verbose = FALSE
  )

  expect_s3_class(split, "ukb_ml_survival_split")
  expect_equal(split$time_var, "time")
  expect_equal(split$event_var, "event")
  expect_equal(nrow(split$train) + nrow(split$validation) + nrow(split$test), nrow(df))
  expect_true(all(c(0, 1) %in% unique(split$train$event)))
})

test_that("ukb_ml_survival_as_split validates manual split overlap", {
  df <- simulate_survival_ml_data(n = 30)

  expect_error(
    ukb_ml_survival_as_split(
      train_data = df[1:20, ],
      test_data = df[20:30, ],
      time = "time",
      event = "event",
      id_col = "eid"
    ),
    "Overlapping IDs"
  )
})

test_that("ukb_ml_survival_workflow runs Cox frozen-test workflow", {
  df <- simulate_survival_ml_data()

  wf <- ukb_ml_survival_workflow(
    survival::Surv(time, event) ~ age + bmi + sex,
    data = df,
    model = "cox",
    split_params = list(
      split = "train_valid_test",
      train_ratio = 0.65,
      validation_ratio = 0.15,
      test_ratio = 0.20
    ),
    feature_select = "filter",
    feature_params = list(max_features = 3),
    tune = TRUE,
    tune_params = list(resampling = "validation"),
    evaluation_times = c(2, 5),
    seed = 42,
    verbose = FALSE
  )

  expect_s3_class(wf, "ukb_ml_survival_workflow")
  expect_s3_class(wf$final_model, "ukb_ml_survival_final")
  expect_true(is.finite(wf$final_test_metrics[["c_index"]]))
  expect_true(all(wf$final_test_predictions$t2 >= 0 & wf$final_test_predictions$t2 <= 1))
  expect_equal(nrow(wf$final_test_predictions), nrow(wf$split$test))

  pred <- predict(wf, times = c(2, 5), type = "survival")
  expect_equal(dim(pred), c(nrow(wf$split$test), 2L))
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("ukb_ml_survival_workflow supports user-supplied split", {
  df <- simulate_survival_ml_data()
  train_idx <- 1:100
  valid_idx <- 101:125
  test_idx <- 126:160

  split <- ukb_ml_survival_as_split(
    train_data = df[train_idx, ],
    validation_data = df[valid_idx, ],
    test_data = df[test_idx, ],
    time = "time",
    event = "event",
    id_col = "eid"
  )

  wf <- ukb_ml_survival_workflow(
    survival::Surv(time, event) ~ age + bmi + sex,
    split = split,
    model = "cox",
    feature_select = "none",
    tune = FALSE,
    evaluation_times = c(3),
    seed = 42,
    verbose = FALSE
  )

  expect_s3_class(wf, "ukb_ml_survival_workflow")
  expect_true(is.finite(wf$final_test_metrics[["c_index"]]))
  expect_equal(nrow(wf$final_test_predictions), length(test_idx))
})
