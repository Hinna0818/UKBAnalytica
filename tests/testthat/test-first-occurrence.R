if (!exists("extract_cases_by_source", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "date_utils.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "first_occurrence.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "disease_definitions.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ICD_diagnose.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "self_report.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "death.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "opcs4.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "case_extraction.R"), local = FALSE)
}

test_that("create_disease_definition stores First Occurrence field metadata", {
  def <- create_disease_definition(
    name = "Myocardial Infarction",
    icd10_pattern = "^(I21|I22)",
    first_occurrence_fields = c(131298, 131300)
  )

  expect_equal(def$first_occurrence_fields, c(131298L, 131300L))
  expect_null(def$first_occurrence_source_fields)
})

test_that("FirstOccurrence source extracts earliest valid UKB date", {
  dt <- data.table::data.table(
    eid = 1:5,
    p53_i0 = as.Date(rep("2010-01-01", 5)),
    p131298_i0 = as.Date(c("2009-01-01", "2015-01-01", "1900-01-01", "2020-01-01", NA)),
    p131299_i0 = c("HES", "Death", "HES", "Death", NA),
    p131300 = as.Date(c(NA, NA, NA, "2018-01-01", NA)),
    p131301 = c(NA, NA, NA, "PrimaryCare", NA)
  )

  disease_definitions <- list(
    MI = create_disease_definition(
      name = "Myocardial Infarction",
      first_occurrence_fields = c(131298, 131300)
    )
  )

  res <- suppressMessages(extract_cases_by_source(
    dt = dt,
    disease_definitions = disease_definitions,
    sources = "FirstOccurrence",
    censor_date = as.Date("2022-01-01")
  ))

  expect_equal(res$eid, c(1L, 2L, 4L))
  expect_equal(res$earliest_date, as.Date(c("2009-01-01", "2015-01-01", "2018-01-01")))
  expect_equal(
    res$diagnosis_source,
    c("FirstOccurrence_HES", "FirstOccurrence_Death", "FirstOccurrence_PrimaryCare")
  )
  expect_equal(res$prevalent_case, c(TRUE, FALSE, FALSE))
  expect_equal(res$status, c(0L, 1L, 1L))
})

test_that("FirstOccurrence works without source columns", {
  dt <- data.table::data.table(
    eid = 1:2,
    p53_i0 = as.Date(c("2010-01-01", "2010-01-01")),
    p131382_i0 = as.Date(c("2012-06-01", NA))
  )

  disease_definitions <- list(
    AA = create_disease_definition(
      name = "Aortic Aneurysm",
      first_occurrence_fields = 131382
    )
  )

  res <- suppressMessages(extract_cases_by_source(
    dt = dt,
    disease_definitions = disease_definitions,
    sources = "FirstOccurrence"
  ))

  expect_equal(res$eid, 1L)
  expect_equal(res$diagnosis_source, "FirstOccurrence")
  expect_equal(res$status, 1L)
})
