test_that("training-validation Cox workflow standardizes and compares results", {
  set.seed(101)
  n <- 120
  dat <- data.frame(
    time = rexp(n, 0.08),
    status = rbinom(n, 1, 0.35),
    age = rnorm(n, 60, 7),
    sex = factor(sample(c("F", "M"), n, TRUE)),
    olink_instance_0.gdf15 = rnorm(n),
    olink_instance_0.ager = rnorm(n),
    olink_instance_0.il6 = rnorm(n),
    olink_instance_0.tnf = rnorm(n)
  )
  proteins <- grep("^olink", names(dat), value = TRUE)

  res <- ukb_train_validation_cox(
    train_data = dat[1:80, ],
    validation_data = dat[81:120, ],
    main_vars = proteins,
    covariates = c("age", "sex"),
    endpoint = c("time", "status"),
    ties = "efron",
    add_protein_annotation = TRUE
  )

  expect_true(all(c("train_results", "validation_results", "comparison") %in% names(res)))
  expect_equal(nrow(res$scaling_parameters), length(proteins))
  expect_true(all(c("p_bh", "p_bonferroni", "significant_bonferroni") %in% names(res$train_results)))
  expect_true(all(c("train_logHR", "valid_logHR", "same_direction") %in% names(res$comparison)))
  expect_s3_class(plot_cox_loghr_correlation(res$comparison), "ggplot")
})

test_that("top HR helpers return source data and plot object", {
  d <- data.frame(
    variable = paste0("x", 1:6),
    HR = c(1.6, 1.3, 1.1, 0.9, 0.8, 0.7),
    lower95 = c(1.2, 1.1, 0.9, 0.7, 0.6, 0.5),
    upper95 = c(2.0, 1.6, 1.3, 1.1, 0.95, 0.9),
    pvalue = c(0.001, 0.01, 0.2, 0.3, 0.02, 0.001),
    p_bonferroni = c(0.006, 0.06, 1, 1, 0.12, 0.006),
    gene_symbol = paste0("G", 1:6)
  )

  top <- ukb_top_hr_results(d, n_each_direction = 1, p_col = "p_bonferroni")
  expect_equal(nrow(top), 2)
  expect_true(all(c("label", "direction") %in% names(top)))
  expect_s3_class(plot_top_hr_bars(top, facet_col = NULL), "ggplot")
})

test_that("sensitivity Cox comparison helpers support repeated analyses", {
  main <- data.frame(
    variable = paste0("x", 1:5),
    HR = c(1.5, 1.2, 0.8, 0.7, 1.1),
    pvalue = c(0.001, 0.01, 0.02, 0.20, 0.50)
  )
  sens <- rbind(
    data.frame(sensitivity = "Lag 2 years", variable = paste0("x", 1:5), HR = c(1.4, 1.1, 0.85, 0.75, 1.0), pvalue = c(0.002, 0.04, 0.03, 0.3, 0.8)),
    data.frame(sensitivity = "Lag 4 years", variable = paste0("x", 1:5), HR = c(1.3, 1.1, 0.9, 0.8, 1.05), pvalue = c(0.006, 0.05, 0.05, 0.4, 0.7))
  )

  cmp <- ukb_compare_sensitivity_cox(main, sens)
  expect_true(all(c("comparison", "correlation_summary") %in% names(cmp)))
  expect_equal(length(unique(cmp$comparison$sensitivity)), 2)
  expect_true(all(c("main_logHR", "sensitivity_logHR", "same_direction") %in% names(cmp$comparison)))
  expect_s3_class(plot_cox_sensitivity_correlation(cmp$comparison), "ggplot")
})

test_that("scaling helper accepts native and legacy parameter formats", {
  dat <- data.frame(x = 1:4, y = c(2, 4, 6, 8))
  native <- ukb_standardize_by_train(dat, variables = c("x", "y"))$parameters
  scaled_native <- ukb_scale_with_parameters(dat, native)
  expect_equal(round(mean(scaled_native$x), 10), 0)

  legacy <- data.frame(
    protein = rep(c("x", "y"), each = 2),
    statistic = rep(c("mean", "sd"), 2),
    value = c(mean(dat$x), stats::sd(dat$x), mean(dat$y), stats::sd(dat$y))
  )
  scaled_legacy <- ukb_scale_with_parameters(dat, legacy)
  expect_equal(unname(scaled_legacy$x), unname(scaled_native$x))
})
