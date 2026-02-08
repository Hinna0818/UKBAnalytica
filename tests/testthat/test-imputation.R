test_that("impute_mice_merge returns m completed datasets", {
  skip_if_not_installed("mice")

  set.seed(1)
  df <- data.frame(
    eid = 1:50,
    outcome = rnorm(50),
    x = rnorm(50),
    y = rnorm(50),
    stringsAsFactors = FALSE
  )
  df$x[sample.int(50, 5)] <- NA
  df$y[sample.int(50, 5)] <- NA

  res <- run_imputation(
    data = df,
    id_col = "eid",
    vars = c("x", "y"),
    m = 2,
    maxit = 2,
    seed = 123,
    print = FALSE
  )

  expect_true(is.list(res))
  expect_true(all(c("imp", "data_list") %in% names(res)))
  expect_length(res$data_list, 2)

  for (i in seq_along(res$data_list)) {
    dfi <- res$data_list[[i]]
    expect_true(all(c("eid", "outcome", "x", "y") %in% names(dfi)))
    expect_equal(nrow(dfi), nrow(df))
    expect_setequal(dfi$eid, df$eid)
  }
})


test_that("impute_mice_merge can treat factor variables", {
  skip_if_not_installed("mice")

  set.seed(1)
  df <- data.frame(
    eid = 1:30,
    sex = sample(c("M", "F", NA_character_), 30, replace = TRUE),
    age = sample(c(40:60, NA_integer_), 30, replace = TRUE),
    stringsAsFactors = FALSE
  )

  res <- run_imputation(
    data = df,
    id_col = "eid",
    vars = c("sex", "age"),
    factor_vars = "sex",
    m = 1,
    maxit = 1,
    seed = 1,
    print = FALSE
  )

  d1 <- res$data_list[[1]]
  expect_true(is.factor(d1$sex) || is.character(d1$sex))
  expect_true(is.numeric(d1$age))
})
