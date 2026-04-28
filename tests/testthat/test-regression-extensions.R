if (!exists("runmulti_cox_zph", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "regression.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "sensitivity.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "regression_extensions.R"), local = FALSE)
}

skip_if_not_installed("survival")

make_regression_extension_test_data <- function() {
  set.seed(42)
  n <- 200
  dt <- data.frame(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.55),
    age = rnorm(n, 60, 8),
    sex = factor(sample(c("F", "M"), n, TRUE)),
    ph_ecog = sample(0:2, n, TRUE),
    stringsAsFactors = FALSE
  )

  quartiles <- stats::quantile(dt$age, probs = c(0, 0.25, 0.5, 0.75, 1))
  breaks <- unique(quartiles)
  dt$age_q <- cut(
    dt$age,
    breaks = breaks,
    include.lowest = TRUE,
    labels = paste0("Q", seq_len(length(breaks) - 1))
  )

  dt$fg_status <- sample(
    c(0, 1, 2),
    nrow(dt),
    replace = TRUE,
    prob = c(0.55, 0.3, 0.15)
  )

  early_event_idx <- which(dt$status == 1 & dt$time <= 1)
  if (length(early_event_idx) == 0) {
    dt$status[[1]] <- 1
    dt$time[[1]] <- 0.5
  }

  dt
}

test_that("runmulti_cox_zph returns effect and PH diagnostics", {
  dt <- make_regression_extension_test_data()

  res <- runmulti_cox_zph(
    data = dt,
    main_var = c("age", "sex"),
    covariates = "ph_ecog",
    endpoint = c("time", "status")
  )

  expect_s3_class(res, "data.frame")
  expect_equal(sort(unique(res$variable)), c("age", "sex"))
  expect_true(all(c("HR", "zph_pvalue", "zph_global_pvalue") %in% names(res)))
  expect_true(all(is.logical(res$ph_violation)))
})

test_that("runmulti_trend returns grouped estimates with p for trend", {
  dt <- make_regression_extension_test_data()

  res <- runmulti_trend(
    data = dt,
    main_var = "age_q",
    covariates = "ph_ecog",
    model_type = "cox",
    endpoint = c("time", "status")
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("estimate", "p_trend", "trend_estimate", "score_method") %in% names(res)))
  expect_equal(unique(res$variable), "age_q")
  expect_true(length(unique(res$p_trend)) == 1)
})

test_that("runmulti_competing returns subdistribution hazard ratios", {
  dt <- make_regression_extension_test_data()

  res <- runmulti_competing(
    data = dt,
    main_var = "age_q",
    covariates = "ph_ecog",
    time_col = "time",
    event_col = "fg_status"
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("SHR", "n_compete") %in% names(res)))
  expect_equal(unique(res$variable), "age_q")
  expect_true(all(res$n_compete >= 0))
})

test_that("runmulti_cox_lag returns lag-specific HR tables", {
  dt <- make_regression_extension_test_data()

  res <- runmulti_cox_lag(
    data = dt,
    main_var = c("age", "sex"),
    covariates = "ph_ecog",
    endpoint = c("time", "status"),
    lag_years = c(0, 1),
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("lag_years", "n_input", "n_removed", "HR") %in% names(res)))
  expect_equal(sort(unique(res$lag_years)), c(0, 1))
  expect_true(all(res$n_input == nrow(dt)))
  expect_true(any(res$n_removed > 0))
})
