test_that("basic publication plotting helpers return ggplot objects", {
  set.seed(42)
  d <- data.frame(
    group = rep(c("A", "B", "C"), each = 20),
    subgroup = rep(c("Low", "High"), 30),
    x = rnorm(60),
    y = rnorm(60),
    value = rnorm(60),
    weight = runif(60, 0.5, 2)
  )
  heat <- aggregate(value ~ group + subgroup, d, mean)

  expect_s3_class(plot_heatmap(heat, x = "group", y = "subgroup", fill = "value"), "ggplot")
  expect_s3_class(plot_heatmap(heat, x = "group", y = "subgroup", fill = "value", show_values = TRUE), "ggplot")
  expect_s3_class(plot_stacked_bar(d, x = "group", fill = "subgroup"), "ggplot")
  expect_s3_class(plot_stacked_bar(d, x = "group", fill = "subgroup", weight = "weight", position = "stack"), "ggplot")
  expect_s3_class(plot_violin(d, x = "group", y = "value"), "ggplot")
  expect_s3_class(plot_scatter(d, x = "x", y = "y", color = "group"), "ggplot")
})

test_that("basic publication plotting helpers validate input columns", {
  d <- data.frame(x = 1:3, y = 1:3)
  expect_error(plot_heatmap(d, x = "x", y = "y", fill = "z"), "missing column")
  expect_error(plot_scatter(d, x = "x", y = "z"), "missing column")
})
