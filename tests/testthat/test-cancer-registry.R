if (!exists("parse_cancer_registry", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "date_utils.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "cancer_registry.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "first_occurrence.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "disease_definitions.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ICD_diagnose.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "self_report.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "death.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "opcs4.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "case_extraction.R"), local = FALSE)
}

test_that("parse_cancer_registry extracts aligned cancer records", {
  dt <- data.table::data.table(
    eid = 1:3,
    p40006_i0 = c("C34.1", "C50.9 Breast", ""),
    p40005_i0 = as.Date(c("2015-01-01", "2016-01-01", NA)),
    p40011_i0 = c("8140", "8500", NA),
    p40012_i0 = c("3", "3", NA),
    p40006_i1 = c("D02.2", NA, "C34"),
    p40005_i1 = as.Date(c("2014-01-01", NA, "2018-05-01")),
    p40012_i1 = c("2", NA, "3")
  )

  res <- parse_cancer_registry(dt)

  expect_equal(nrow(res), 4L)
  expect_equal(res$cancer_icd10_code, c("C341", "C509", "D022", "C34"))
  expect_equal(res$cancer_behaviour, c(3L, 3L, 2L, 3L))
})

test_that("CancerRegistry source contributes cancer cases with behaviour filter", {
  dt <- data.table::data.table(
    eid = 1:4,
    p53_i0 = as.Date(rep("2010-01-01", 4)),
    p40006_i0 = c("C34.1", "C34.9", "C34.1", "C50.9"),
    p40005_i0 = as.Date(c("2009-01-01", "2015-01-01", "2018-01-01", "2017-01-01")),
    p40012_i0 = c(3L, 3L, 2L, 3L)
  )

  defs <- list(
    Lung_Cancer = create_disease_definition(
      name = "Lung Cancer",
      cancer_icd10_pattern = "^C34",
      cancer_behaviour = 3L
    )
  )

  res <- suppressMessages(extract_cases_by_source(
    dt = dt,
    disease_definitions = defs,
    sources = "CancerRegistry",
    censor_date = as.Date("2020-01-01")
  ))

  expect_equal(res$eid, c(1L, 2L))
  expect_equal(res$diagnosis_source, c("CancerRegistry", "CancerRegistry"))
  expect_equal(res$prevalent_case, c(TRUE, FALSE))
  expect_equal(res$status, c(0L, 1L))
})

test_that("CancerRegistry and ICD10 sources use earliest date across union", {
  dt <- data.table::data.table(
    eid = 1L,
    p53_i0 = as.Date("2010-01-01"),
    p41270 = "['C349']",
    p41280_a0 = as.Date("2018-01-01"),
    p40006_i0 = "C34.9",
    p40005_i0 = as.Date("2016-01-01"),
    p40012_i0 = 3L
  )

  defs <- list(Lung_Cancer = get_predefined_diseases()$Lung_Cancer)

  res <- suppressMessages(extract_cases_by_source(
    dt = dt,
    disease_definitions = defs,
    sources = c("ICD10", "CancerRegistry"),
    censor_date = as.Date("2020-01-01")
  ))

  expect_equal(nrow(res), 1L)
  expect_equal(res$earliest_date, as.Date("2016-01-01"))
  expect_equal(res$diagnosis_source, "CancerRegistry")
})
