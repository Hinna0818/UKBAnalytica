test_that("run_correlation works with valid input", {
  skip_if_not_installed("data.table")
  
  set.seed(123)
  df <- data.frame(
    x = rnorm(50),
    y = rnorm(50),
    z = rnorm(50)
  )
  
  corr_mat <- run_correlation(df, vars = c("x", "y", "z"), threshold = NULL)
  
  expect_true(is.matrix(corr_mat))
  expect_equal(dim(corr_mat), c(3, 3))
  expect_equal(rownames(corr_mat), c("x", "y", "z"))
  expect_equal(colnames(corr_mat), c("x", "y", "z"))
  expect_true(all(abs(diag(corr_mat) - 1) < 1e-10))
})


test_that("run_correlation handles missing variables gracefully", {
  df <- data.frame(x = 1:10, y = 1:10)
  
  expect_error(
    run_correlation(df, vars = c("x", "y", "missing")),
    "The following variables are not in the data frame: missing"
  )
})


test_that("plot_correlation creates a ggplot object", {
  skip_if_not_installed("ggplot2")
  
  set.seed(123)
  df <- data.frame(
    x = rnorm(50),
    y = rnorm(50),
    z = rnorm(50)
  )
  
  corr_mat <- run_correlation(df, vars = c("x", "y", "z"), threshold = NULL)
  p <- plot_correlation(corr_mat, show_values = TRUE)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_correlation works with upper triangle option", {
  skip_if_not_installed("ggplot2")
  
  corr_mat <- matrix(c(1, 0.5, 0.3,
                       0.5, 1, 0.7,
                       0.3, 0.7, 1), nrow = 3)
  rownames(corr_mat) <- colnames(corr_mat) <- c("A", "B", "C")
  
  p <- plot_correlation(
    corr_mat,
    upper_triangle = TRUE,
    show_values = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})
