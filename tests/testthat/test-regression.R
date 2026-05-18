if (!exists("run_regression", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "regression.R"), local = FALSE)
}

skip_if_not_installed("survival")
library(survival)

make_reg_test_data <- function() {
  set.seed(123)
  n <- 150
  data.frame(
    time    = rexp(n, 0.1),
    status  = rbinom(n, 1, 0.5),
    outcome = rnorm(n),
    bin_out = rbinom(n, 1, 0.4),
    count   = rpois(n, lambda = 5),
    pos_out = abs(rnorm(n)) + 0.1,
    x1      = rnorm(n),
    x2      = rnorm(n),
    cov1    = rnorm(n),
    stringsAsFactors = FALSE
  )
}

# ── run_regression (cox) ─────────────────────────────────────────────────────

test_that("run_regression cox returns correct structure", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = c("x1", "x2"),
                        type = "cox", endpoint = c("time", "status"))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expected_cols <- c("variable", "HR", "lower95", "upper95", "pvalue",
                     "n", "n_event")
  expect_true(all(expected_cols %in% names(res)))
  expect_setequal(res$variable, c("x1", "x2"))
  expect_true(all(res$HR > 0))
})

test_that("run_regression cox with covariates works", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = "x1", covariates = "cov1",
                        type = "cox", endpoint = c("time", "status"))

  expect_equal(nrow(res), 1)
  expect_equal(res$variable, "x1")
})

# ── run_regression (lm) ──────────────────────────────────────────────────────

test_that("run_regression lm returns correct structure", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = c("x1", "x2"),
                        type = "lm", outcome = "outcome")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("variable", "beta", "lower95", "upper95", "pvalue")
                  %in% names(res)))
  expect_setequal(res$variable, c("x1", "x2"))
})

test_that("run_regression lm with covariates works", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = "x1", covariates = "cov1",
                        type = "lm", outcome = "outcome")

  expect_equal(nrow(res), 1)
})

test_that("run_regression lm errors when outcome missing", {
  dt <- make_reg_test_data()
  expect_error(run_regression(dt, main_var = "x1", type = "lm"),
               "'outcome' must be specified")
})

# ── run_regression (logit) ───────────────────────────────────────────────────

test_that("run_regression logit returns correct structure", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = c("x1", "x2"),
                        type = "logit", outcome = "bin_out")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("variable", "OR", "lower95", "upper95", "pvalue")
                  %in% names(res)))
  expect_true(all(res$OR > 0))
})

test_that("run_regression logit with covariates works", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = "x1", covariates = "cov1",
                        type = "logit", outcome = "bin_out")

  expect_equal(nrow(res), 1)
})

test_that("run_regression logit errors when outcome missing", {
  dt <- make_reg_test_data()
  expect_error(run_regression(dt, main_var = "x1", type = "logit"),
               "'outcome' must be specified")
})

# ── runmulti_glm ─────────────────────────────────────────────────────────────

test_that("runmulti_glm poisson returns correct structure", {
  dt  <- make_reg_test_data()
  res <- runmulti_glm(dt, main_var = c("x1", "x2"),
                      family = "poisson", outcome = "count")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("variable", "family", "link", "beta",
                    "lower95", "upper95", "pvalue", "n") %in% names(res)))
  expect_setequal(res$variable, c("x1", "x2"))
  expect_true(all(res$family == "poisson"))
  expect_true(all(res$link == "log"))
  expect_equal(res$n[[1]], nrow(dt))
})

test_that("runmulti_glm poisson with covariates works", {
  dt  <- make_reg_test_data()
  res <- runmulti_glm(dt, main_var = "x1", covariates = "cov1",
                      family = "poisson", outcome = "count")

  expect_equal(nrow(res), 1)
  expect_equal(res$variable, "x1")
})

test_that("runmulti_glm quasipoisson uses Wald CI", {
  dt  <- make_reg_test_data()
  res <- runmulti_glm(dt, main_var = "x1",
                      family = "quasipoisson", outcome = "count")

  expect_s3_class(res, "data.frame")
  expect_true(all(c("beta", "lower95", "upper95") %in% names(res)))
  expect_true(all(is.finite(res$lower95)))
  expect_true(all(is.finite(res$upper95)))
  expect_true(all(res$lower95 < res$beta))
  expect_true(all(res$upper95 > res$beta))
})

test_that("runmulti_glm Gamma with log link works", {
  dt  <- make_reg_test_data()
  res <- runmulti_glm(dt, main_var = c("x1", "x2"),
                      family = stats::Gamma(link = "log"),
                      outcome = "pos_out")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(res$link == "log"))
})

test_that("runmulti_glm accepts function form of family", {
  dt  <- make_reg_test_data()
  res <- runmulti_glm(dt, main_var = "x1",
                      family = stats::poisson, outcome = "count")

  expect_equal(res$family[[1]], "poisson")
})

test_that("runmulti_glm gaussian matches runmulti_lm betas", {
  dt  <- make_reg_test_data()
  res_glm <- runmulti_glm(dt, main_var = "x1",
                          family = "gaussian", outcome = "outcome")
  res_lm  <- runmulti_lm(dt, main_var = "x1", outcome = "outcome")

  expect_equal(round(res_glm$beta, 3), res_lm$beta)
})

test_that("runmulti_glm errors on unknown family string", {
  dt <- make_reg_test_data()
  expect_error(
    runmulti_glm(dt, main_var = "x1", family = "notafamily",
                 outcome = "count"),
    "Unknown GLM family"
  )
})

# ── run_regression (glm) ─────────────────────────────────────────────────────

test_that("run_regression glm dispatches to runmulti_glm", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = c("x1", "x2"),
                        type = "glm", family = "poisson",
                        outcome = "count")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true("family" %in% names(res))
})

test_that("run_regression glm errors when outcome missing", {
  dt <- make_reg_test_data()
  expect_error(
    run_regression(dt, main_var = "x1", type = "glm"),
    "'outcome' must be specified"
  )
})

# ── type validation ───────────────────────────────────────────────────────────

test_that("run_regression errors on invalid type", {
  dt <- make_reg_test_data()
  expect_error(run_regression(dt, main_var = "x1", type = "tobit"),
               "should be one of")
})

# ── runmulti_negbin ───────────────────────────────────────────────────────────

test_that("runmulti_negbin returns correct structure", {
  dt  <- make_reg_test_data()
  res <- runmulti_negbin(dt, main_var = c("x1", "x2"), outcome = "count")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("variable", "IRR", "lower95", "upper95",
                    "pvalue", "theta", "n") %in% names(res)))
  expect_setequal(res$variable, c("x1", "x2"))
  expect_true(all(res$IRR > 0))
  expect_true(all(res$theta > 0))
  expect_equal(res$n[[1]], nrow(dt))
})

test_that("runmulti_negbin with covariates works", {
  dt  <- make_reg_test_data()
  res <- runmulti_negbin(dt, main_var = "x1",
                         covariates = "cov1", outcome = "count")

  expect_equal(nrow(res), 1)
  expect_equal(res$variable, "x1")
})

test_that("runmulti_negbin IRR CIs are on exp scale", {
  dt  <- make_reg_test_data()
  res <- runmulti_negbin(dt, main_var = "x1", outcome = "count")

  expect_true(res$lower95 < res$IRR)
  expect_true(res$upper95 > res$IRR)
})

test_that("run_regression negbin dispatches correctly", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = c("x1", "x2"),
                        type = "negbin", outcome = "count")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true("IRR" %in% names(res))
})

test_that("run_regression negbin errors when outcome missing", {
  dt <- make_reg_test_data()
  expect_error(run_regression(dt, main_var = "x1", type = "negbin"),
               "'outcome' must be specified")
})

# ── runmulti_gam ─────────────────────────────────────────────────────────────

skip_if_not_installed("mgcv")

test_that("runmulti_gam smooth returns correct structure", {
  dt  <- make_reg_test_data()
  res <- runmulti_gam(dt, main_var = c("x1", "x2"), outcome = "outcome")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("variable", "edf", "ref_df", "stat", "stat_type",
                    "pvalue", "family", "link", "n") %in% names(res)))
  expect_setequal(res$variable, c("x1", "x2"))
  expect_true(all(res$edf >= 1))
  expect_equal(res$n[[1]], nrow(dt))
})

test_that("runmulti_gam smooth with covariates works", {
  dt  <- make_reg_test_data()
  res <- runmulti_gam(dt, main_var = "x1", covariates = "cov1",
                      outcome = "outcome")

  expect_equal(nrow(res), 1)
  expect_equal(res$variable, "x1")
})

test_that("runmulti_gam linear mode returns beta structure", {
  dt  <- make_reg_test_data()
  res <- runmulti_gam(dt, main_var = c("x1", "x2"),
                      outcome = "outcome", smooth = FALSE)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("variable", "beta", "lower95", "upper95",
                    "pvalue", "family", "link") %in% names(res)))
  expect_true(all(is.finite(res$beta)))
})

test_that("runmulti_gam poisson family works", {
  dt  <- make_reg_test_data()
  res <- runmulti_gam(dt, main_var = "x1", outcome = "count",
                      family = "poisson")

  expect_equal(res$family[[1]], "poisson")
  expect_equal(res$link[[1]], "log")
})

test_that("run_regression gam smooth dispatches correctly", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = c("x1", "x2"),
                        type = "gam", outcome = "outcome")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true("edf" %in% names(res))
})

test_that("run_regression gam linear dispatches correctly", {
  dt  <- make_reg_test_data()
  res <- run_regression(dt, main_var = "x1",
                        type = "gam", outcome = "outcome", smooth = FALSE)

  expect_true("beta" %in% names(res))
})

test_that("run_regression gam errors when outcome missing", {
  dt <- make_reg_test_data()
  expect_error(run_regression(dt, main_var = "x1", type = "gam"),
               "'outcome' must be specified")
})
