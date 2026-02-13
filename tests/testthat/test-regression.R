test_that("runmulti_lm works for univariate models", {
  data <- mtcars
  res <- runmulti_lm(data, main_var = c("hp", "wt"), outcome = "mpg")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_equal(res$variable, c("hp", "wt"))
  expect_true(all(c("variable", "beta", "lower95", "upper95", "pvalue") %in% names(res)))
  # hp should have negative beta (more hp -> less mpg)
  expect_true(res$beta[res$variable == "hp"] < 0)
  expect_true(res$beta[res$variable == "wt"] < 0)
})

test_that("runmulti_lm works for multivariate models", {
  data <- mtcars
  res <- runmulti_lm(data,
                     main_var = c("hp", "wt"),
                     covariates = c("cyl"),
                     outcome = "mpg")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(res$pvalue > 0 & res$pvalue <= 1))
  expect_true(all(res$lower95 < res$upper95))
})

test_that("runmulti_logit works for univariate models", {
  data <- mtcars
  data$vs <- as.integer(data$vs)

  res <- runmulti_logit(data, main_var = c("hp", "wt"), outcome = "vs")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_equal(res$variable, c("hp", "wt"))
  expect_true(all(c("variable", "OR", "lower95", "upper95", "pvalue") %in% names(res)))
  # OR should be positive
  expect_true(all(res$OR > 0))
  expect_true(all(res$lower95 < res$upper95))
})

test_that("runmulti_logit works for multivariate models", {
  data <- mtcars
  data$vs <- as.integer(data$vs)

  res <- runmulti_logit(data,
                        main_var = c("hp", "wt"),
                        covariates = c("cyl"),
                        outcome = "vs")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(res$OR > 0))
})

test_that("runmulti_logit rejects non-binary outcome", {
  data <- mtcars
  expect_error(
    runmulti_logit(data, main_var = c("hp"), outcome = "mpg"),
    "must be binary"
  )
})

test_that("runmulti_cox works for univariate models", {
  skip_if_not_installed("survival")
  library(survival)

  data <- lung
  data <- data[complete.cases(data[, c("time", "status", "age", "sex")]), ]

  res <- runmulti_cox(data, main_var = c("age", "sex"), endpoint = c("time", "status"))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_equal(res$variable, c("age", "sex"))
  expect_true(all(c("variable", "HR", "lower95", "upper95", "pvalue") %in% names(res)))
  expect_true(all(res$HR > 0))
  expect_true(all(res$lower95 < res$upper95))
})

test_that("runmulti_cox works for multivariate models", {
  skip_if_not_installed("survival")
  library(survival)

  data <- lung
  data <- data[complete.cases(data[, c("time", "status", "age", "sex", "ph.ecog")]), ]

  res <- runmulti_cox(data,
                      main_var = c("age", "sex"),
                      covariates = c("ph.ecog"),
                      endpoint = c("time", "status"))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(res$HR > 0))
})

test_that("all regression functions reject invalid inputs", {
  data <- mtcars

  # missing variable
  expect_error(
    runmulti_lm(data, main_var = c("nonexistent"), outcome = "mpg"),
    "not found in data"
  )

  # non-data.frame
  expect_error(
    runmulti_lm("not_a_df", main_var = c("hp"), outcome = "mpg"),
    "must be a data.frame"
  )

  # empty main_var
  expect_error(
    runmulti_logit(data, main_var = character(0), outcome = "vs"),
    "non-empty"
  )
})
