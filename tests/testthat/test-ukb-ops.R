if (!exists("ukb_clean_missing", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "ukb_ops.R"), local = FALSE)
}

test_that("ukb_clean_missing converts UKB non-response values", {
  dt <- data.frame(
    eid = 1:4,
    smoking = c("Never", "Do not know", " Prefer not to answer ", ""),
    alcohol = c(1, -1, -3, 2),
    stringsAsFactors = FALSE
  )

  cleaned <- suppressMessages(ukb_clean_missing(dt, verbose = TRUE))
  summary_dt <- attr(cleaned, "ukb_clean_missing_summary")

  expect_true(data.table::is.data.table(cleaned))
  expect_equal(cleaned$smoking, c("Never", NA, NA, NA))
  expect_equal(cleaned$alcohol, c(1, NA, NA, 2))
  expect_equal(sum(summary_dt$label_replaced), 2L)
  expect_equal(sum(summary_dt$numeric_code_to_na), 2L)
})

test_that("ukb_clean_missing can retain informative missing as Unknown", {
  dt <- data.frame(
    status = factor(c("Case", "Do not know", "Control")),
    stringsAsFactors = TRUE
  )

  cleaned <- suppressMessages(ukb_clean_missing(dt, action = "unknown"))

  expect_equal(cleaned$status, c("Case", "Unknown", "Control"))
})

test_that("ukb_snapshot records and retrieves cohort history", {
  suppressMessages(ukb_snapshot(id = "test-ops", reset = TRUE))

  dt1 <- data.frame(eid = 1:4, x = c(1, 2, NA, 4))
  dt2 <- dt1[!is.na(dt1$x), ]

  suppressMessages(ukb_snapshot(dt1, label = "raw", id = "test-ops"))
  hist <- suppressMessages(ukb_snapshot(dt2, label = "complete x", id = "test-ops"))

  expect_equal(nrow(hist), 2L)
  expect_equal(hist$label, c("raw", "complete x"))
  expect_equal(hist$nrow, c(4L, 3L))
  expect_equal(hist$row_delta, c(NA_integer_, -1L))
  expect_equal(hist$n_missing_cols, c(1L, 0L))

  retrieved <- ukb_snapshot(id = "test-ops", verbose = FALSE)
  expect_equal(retrieved, hist)

  suppressMessages(ukb_snapshot(id = "test-ops", reset = TRUE))
})
