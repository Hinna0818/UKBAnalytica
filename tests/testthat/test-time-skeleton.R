if (!exists("ukb_time_skeleton", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "date_utils.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "time_skeleton.R"), local = FALSE)
}

test_that("ukb_time_skeleton builds administrative, death, and loss-to-follow-up times", {
  dat <- data.frame(
    eid = 1:5,
    p53_i0 = as.Date(c("2010-01-01", "2010-01-01", "2010-01-01", "2010-01-01", "2010-01-01")),
    p34 = c(1960, 1950, 1940, 1930, 1970),
    p52 = c(6, 7, 8, 9, 10),
    p21022 = c(49, 59, 69, 79, 39),
    p40000_i0 = as.Date(c(NA, "2015-01-01", "2025-01-01", NA, "2009-01-01")),
    p40000_i1 = as.Date(c(NA, NA, NA, NA, NA)),
    p191 = as.Date(c(NA, NA, "2014-01-01", NA, NA))
  )

  out <- ukb_time_skeleton(dat, admin_censor_date = as.Date("2020-12-31"))

  expect_s3_class(out, "ukb_time_skeleton")
  expect_equal(nrow(out), 5)
  expect_true(all(c(
    "baseline_date", "birth_date_approx", "age_at_baseline",
    "death_date", "lost_to_followup_date", "followup_end_date",
    "followup_end_reason", "followup_time_years", "valid_followup"
  ) %in% names(out)))

  expect_equal(out[eid == 1, followup_end_reason], "administrative_censoring")
  expect_equal(out[eid == 2, followup_end_reason], "death")
  expect_equal(out[eid == 3, followup_end_reason], "lost_to_followup")
  expect_equal(out[eid == 4, followup_end_date], as.Date("2020-12-31"))
  expect_true(out[eid == 5, death_before_baseline])
  expect_false(out[eid == 5, valid_followup])
})

test_that("ukb_time_skeleton works when optional birth and loss columns are absent", {
  dat <- data.frame(
    eid = 1:3,
    p53_i0 = as.Date(c("2010-01-01", "2011-01-01", "2012-01-01")),
    p21022 = c(50, 60, 70),
    p40000_i0 = as.Date(c(NA, "2014-01-01", NA))
  )

  out <- ukb_time_skeleton(dat, admin_censor_date = as.Date("2020-12-31"))

  expect_equal(nrow(out), 3)
  expect_true(all(is.na(out$birth_date_approx)))
  expect_equal(out$age_at_baseline, c(50, 60, 70))
  expect_equal(out[eid == 2, followup_end_reason], "death")
  expect_true(all(out$valid_followup))
})

test_that("ukb_time_skeleton accepts participant-prefixed RAP column names", {
  dat <- data.frame(
    participant.eid = 1:2,
    participant.p53_i0 = as.Date(c("2010-01-01", "2011-01-01")),
    participant.p21022 = c(50, 60),
    participant.p40000_i0 = as.Date(c(NA, "2013-01-01"))
  )

  out <- ukb_time_skeleton(dat, admin_censor_date = as.Date("2020-12-31"))

  expect_equal(out$eid, 1:2)
  expect_equal(out[eid == 2, followup_end_reason], "death")
  expect_equal(out[eid == 2, age_at_baseline_source], "participant.p21022")
})

test_that("ukb_time_skeleton validates required inputs", {
  dat <- data.frame(eid = 1:2, p21022 = c(50, 60))

  expect_error(
    ukb_time_skeleton(dat),
    "Baseline date column not found"
  )
  expect_error(
    ukb_time_skeleton(data.frame(p53_i0 = Sys.Date())),
    "ID column not found"
  )
})
