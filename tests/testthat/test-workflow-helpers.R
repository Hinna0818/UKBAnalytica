test_that("RAP environment check and extraction manifest are lightweight", {
  env <- ukb_check_rap_env(verbose = FALSE)
  expect_s3_class(env, "ukb_rap_env")
  expect_true(all(c("check", "status", "message") %in% names(env$checks)))

  manifest <- ukb_create_extraction_manifest(
    field_id = c(31, 21022),
    variable_set = "clinical_core",
    variables = c("sex", "age"),
    purpose = "unit test"
  )
  expect_s3_class(manifest, "ukb_extraction_manifest")
  expect_true(nrow(manifest$fields) > 0)
  expect_true(all(c("source", "field_id", "ukb_col") %in% names(manifest$fields)))

  out <- tempfile(fileext = ".csv")
  expect_invisible(ukb_write_extraction_manifest(manifest, out))
  expect_true(file.exists(out))
  expect_true(file.exists(paste0(tools::file_path_sans_ext(out), "_summary.csv")))
})

test_that("participant flow records sequential exclusions and plots", {
  dat <- data.frame(
    eid = 1:6,
    age = c(50, 55, NA, 61, 70, 64),
    bmi = c(25, 26, 27, NA, 30, 24),
    status = c(0, 1, 0, 1, 0, 1),
    time = c(5, 4, 3, 2, 1, 6)
  )
  flow <- ukb_participant_flow(
    dat,
    steps = list(
      "Complete covariates" = c("age", "bmi"),
      "Valid follow-up" = ~ !is.na(time) & time >= 0
    ),
    id_col = "eid",
    outcome_col = "status"
  )
  expect_s3_class(flow, "ukb_participant_flow")
  expect_equal(flow$n_after[flow$step == "Complete covariates"], 4L)
  expect_equal(sum(attr(flow, "kept_index")), 4L)
  expect_s3_class(plot_participant_flow(flow), "ggplot")
})

test_that("sensitivity suite fits primary and sensitivity Cox models", {
  set.seed(123)
  n <- 120
  dat <- data.frame(
    time = pmin(rexp(n, 0.08), 8),
    exposure = rnorm(n),
    age = rnorm(n, 60, 6),
    sex = rbinom(n, 1, 0.5),
    bmi = rnorm(n, 27, 4),
    tdi = rnorm(n)
  )
  dat$status <- rbinom(n, 1, plogis(-1 + 0.6 * dat$exposure + 0.02 * (dat$age - 60)))
  dat$bmi[1:4] <- NA

  res <- ukb_sensitivity_suite(
    dat,
    exposure = "exposure",
    covariates = c("age", "sex", "bmi"),
    endpoint = c("time", "status"),
    early_event_years = c(1, 2),
    complete_case_covariates = c("age", "sex", "bmi"),
    additional_covariate_sets = list(plus_tdi = "tdi"),
    verbose = FALSE
  )

  expect_s3_class(res, "ukb_sensitivity_suite")
  expect_true(all(c("primary", "complete_case", "exclude_events_le_1y", "adjust_plus_tdi") %in% res$summary$analysis))
  expect_true("exposure" %in% res$summary$term)
  expect_true(all(res$summary$HR > 0, na.rm = TRUE))
})
