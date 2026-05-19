if (!exists("run_rcs", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "rcs.R"), local = FALSE)
}

skip_if_not_installed("survival")
library(survival)

make_rcs_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  data.frame(
    time    = rexp(n, 0.05),
    status  = rbinom(n, 1, 0.4),
    bin_out = rbinom(n, 1, 0.35),
    cont_out = rnorm(n),
    x       = rnorm(n, mean = 40, sd = 10),
    age     = rnorm(n, mean = 55, sd = 8),
    sex     = rbinom(n, 1, 0.5),
    stringsAsFactors = FALSE
  )
}

# ── Input validation ──────────────────────────────────────────────────────────

test_that("run_rcs errors on missing exposure column", {
  dt <- make_rcs_data()
  expect_error(
    run_rcs(dt, exposure = "no_such_col", model_type = "cox",
            endpoint = c("time", "status"), backend = "ns"),
    "not found"
  )
})

test_that("run_rcs errors on missing covariate column", {
  dt <- make_rcs_data()
  expect_error(
    run_rcs(dt, exposure = "x", covariates = "missing_cov",
            model_type = "cox", endpoint = c("time", "status"), backend = "ns"),
    "not found"
  )
})

test_that("run_rcs errors when outcome missing for logistic", {
  dt <- make_rcs_data()
  expect_error(
    run_rcs(dt, exposure = "x", model_type = "logistic", backend = "ns"),
    "'outcome' must be specified"
  )
})

test_that("run_rcs errors when outcome missing for linear", {
  dt <- make_rcs_data()
  expect_error(
    run_rcs(dt, exposure = "x", model_type = "linear", backend = "ns"),
    "'outcome' must be specified"
  )
})

test_that("run_rcs errors when cox endpoint missing", {
  dt <- make_rcs_data()
  expect_error(
    run_rcs(dt, exposure = "x", model_type = "cox",
            endpoint = NULL, backend = "ns"),
    "endpoint"
  )
})

# ── ns backend: Cox ───────────────────────────────────────────────────────────

test_that("ns backend cox returns ukb_rcs class", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_s3_class(fit, "ukb_rcs")
  expect_s3_class(fit, "list")
})

test_that("ns backend cox prediction has correct structure", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_true(is.data.frame(fit$prediction))
  expect_true(all(c("x", "estimate", "lower95", "upper95") %in% names(fit$prediction)))
  expect_equal(nrow(fit$prediction), 200L)
})

test_that("ns backend cox estimate is HR-scale (positive)", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_true(all(fit$prediction$estimate > 0))
  expect_true(all(fit$prediction$lower95  > 0))
  expect_true(all(fit$prediction$upper95  > 0))
})

test_that("ns backend cox CI is correctly ordered", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_true(all(fit$prediction$lower95 <= fit$prediction$estimate + 1e-8))
  expect_true(all(fit$prediction$upper95 >= fit$prediction$estimate - 1e-8))
})

test_that("ns backend cox reference point has estimate ≈ 1", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3, ref_quantile = 0.5)
  ref_row <- fit$prediction[which.min(abs(fit$prediction$x - fit$ref)), ]
  expect_equal(ref_row$estimate, 1, tolerance = 1e-6)
})

test_that("ns backend cox with covariates works", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", covariates = c("age", "sex"),
                 model_type = "cox", endpoint = c("time", "status"),
                 backend = "ns", knots = 3)
  expect_s3_class(fit, "ukb_rcs")
  expect_equal(fit$covariates, c("age", "sex"))
})

test_that("ns backend cox returns n and n_event", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_true(is.numeric(fit$n)       && fit$n > 0)
  expect_true(is.numeric(fit$n_event) && fit$n_event > 0)
})

# ── ns backend: logistic ──────────────────────────────────────────────────────

test_that("ns backend logistic returns OR-scale predictions", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "logistic",
                 outcome = "bin_out", backend = "ns", knots = 3)
  expect_s3_class(fit, "ukb_rcs")
  expect_true(all(fit$prediction$estimate > 0))
  ref_row <- fit$prediction[which.min(abs(fit$prediction$x - fit$ref)), ]
  expect_equal(ref_row$estimate, 1, tolerance = 1e-6)
})

# ── ns backend: linear ────────────────────────────────────────────────────────

test_that("ns backend linear returns beta-difference predictions", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "linear",
                 outcome = "cont_out", backend = "ns", knots = 3)
  expect_s3_class(fit, "ukb_rcs")
  # Beta differences: reference row should be ≈ 0
  ref_row <- fit$prediction[which.min(abs(fit$prediction$x - fit$ref)), ]
  expect_equal(ref_row$estimate, 0, tolerance = 1e-6)
  # CI can be negative (differences around zero)
  expect_true(all(fit$prediction$lower95 <= fit$prediction$upper95 + 1e-8))
})

# ── AIC table & knot selection ────────────────────────────────────────────────

test_that("ns backend returns aic_table when knots auto-selected", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knot_range = 3:5)
  expect_true(is.data.frame(fit$aic_table))
  expect_true(all(c("knots", "AIC") %in% names(fit$aic_table)))
  expect_equal(nrow(fit$aic_table), 3L)
  expect_true(fit$knots %in% 3:5)
})

test_that("ns backend respects user-specified knots", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 4)
  expect_equal(fit$knots, 4L)
})

# ── Return fields ─────────────────────────────────────────────────────────────

test_that("run_rcs returns all required fields", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  required <- c("model", "model_type", "backend", "exposure", "covariates",
                "endpoint", "outcome", "knots", "ref", "n", "n_event",
                "p_overall", "p_nonlinear", "prediction", "distribution",
                "aic_table")
  expect_true(all(required %in% names(fit)))
})

test_that("run_rcs distribution field contains exposure values", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_true(is.data.frame(fit$distribution))
  expect_true("x" %in% names(fit$distribution))
  expect_true(nrow(fit$distribution) > 0)
})

# ── plot_rcs ──────────────────────────────────────────────────────────────────

test_that("plot_rcs returns a ggplot object (ns cox)", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  p <- plot_rcs(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_rcs works with distribution = 'histogram'", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_s3_class(plot_rcs(fit, distribution = "histogram"), "ggplot")
})

test_that("plot_rcs works with distribution = 'density'", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_s3_class(plot_rcs(fit, distribution = "density"), "ggplot")
})

test_that("plot_rcs works with distribution = 'rug'", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_s3_class(plot_rcs(fit, distribution = "rug"), "ggplot")
})

test_that("plot_rcs works with show_distribution = FALSE", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_s3_class(plot_rcs(fit, show_distribution = FALSE), "ggplot")
})

test_that("plot_rcs works for logistic model", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "logistic",
                 outcome = "bin_out", backend = "ns", knots = 3)
  expect_s3_class(plot_rcs(fit), "ggplot")
})

test_that("plot_rcs works for linear model", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "linear",
                 outcome = "cont_out", backend = "ns", knots = 3)
  expect_s3_class(plot_rcs(fit), "ggplot")
})

test_that("plot.ukb_rcs S3 method dispatches correctly", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "ns",
                 knots = 3)
  expect_s3_class(plot(fit), "ggplot")
})

# ── rms backend (skip if not installed) ───────────────────────────────────────

skip_if_not_installed("rms")

test_that("rms backend cox returns ukb_rcs class", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "rms",
                 knots = 3)
  expect_s3_class(fit, "ukb_rcs")
  expect_equal(fit$backend, "rms")
})

test_that("rms backend cox prediction has correct structure", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "rms",
                 knots = 3)
  expect_true(all(c("x", "estimate", "lower95", "upper95") %in% names(fit$prediction)))
  expect_true(all(fit$prediction$estimate > 0))
})

test_that("rms backend aic_table covers knot_range", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "rms",
                 knot_range = 3:5)
  expect_equal(nrow(fit$aic_table), 3L)
  expect_true(fit$knots %in% 3:5)
})

test_that("rms backend plot_rcs returns ggplot", {
  dt  <- make_rcs_data()
  fit <- run_rcs(dt, exposure = "x", model_type = "cox",
                 endpoint = c("time", "status"), backend = "rms",
                 knots = 3)
  expect_s3_class(plot_rcs(fit), "ggplot")
})
