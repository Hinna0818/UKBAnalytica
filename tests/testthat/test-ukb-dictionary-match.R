if (!exists("ukb_query_dictionary", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "field_metadata.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "rap_manifest.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ukb_dictionary_match.R"), local = FALSE)
}

make_test_official_dict <- function(path) {
  x <- data.frame(
    ID = c("eid", "p31", "p21022", "p4080_i0_a0", "p4079_i0_a0", "p41270"),
    title = c(
      "Participant ID",
      "Sex",
      "Age at recruitment",
      "Systolic blood pressure, automated reading | Instance 0 | Array 0",
      "Diastolic blood pressure, automated reading | Instance 0 | Array 0",
      "Diagnoses - ICD10"
    ),
    Instance = c("", "", "", "0", "0", ""),
    Array = c("", "", "", "0", "0", ""),
    name = c(
      "Participant ID",
      "Sex",
      "Age at recruitment",
      "Systolic blood pressure, automated reading",
      "Diastolic blood pressure, automated reading",
      "Diagnoses - ICD10"
    ),
    units = c("", "", "years", "mmHg", "mmHg", ""),
    folder_path = c(
      "Participant Information",
      "Population characteristics > Baseline characteristics",
      "Population characteristics > Baseline characteristics",
      "Assessment centre > Physical measures > Blood pressure",
      "Assessment centre > Physical measures > Blood pressure",
      "Health-related outcomes > Hospital inpatient > Summary Diagnoses"
    ),
    entity = "participant",
    linkout = "",
    type = c("string", "integer", "integer", "integer", "integer", "string"),
    coding_name = c("", "data_coding_9", "", "", "", "data_coding_19"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(x, path, row.names = FALSE)
}

make_test_zh_dict <- function(path) {
  x <- data.frame(
    类别1 = c("参与者信息", "人口特征", "评估中心"),
    类别2 = c("", "基线特征", "身体测量"),
    类别3 = c("", "", "血压"),
    类别4 = c("", "", ""),
    类别5 = c("", "", ""),
    类别6 = c("参与者ID", "性别", "收缩压自动读数|实例0|数组0"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

test_that("ukb_query_dictionary requires RAP by default", {
  official <- tempfile(fileext = ".csv")
  zh <- tempfile(fileext = ".csv")
  make_test_official_dict(official)
  make_test_zh_dict(zh)

  expect_error(
    ukb_query_dictionary("sex", official_dict = official, zh_dict = zh),
    "must be run inside a UK Biobank RAP environment"
  )
})

test_that("ukb_query_dictionary searches English, Chinese, field IDs, and columns", {
  official <- tempfile(fileext = ".csv")
  zh <- tempfile(fileext = ".csv")
  make_test_official_dict(official)
  make_test_zh_dict(zh)

  en <- ukb_query_dictionary("sex", official_dict = official, zh_dict = zh, require_rap = FALSE)
  expect_s3_class(en, "ukb_dictionary_query")
  expect_true(any(en$official_matches$FieldID == "31"))

  cn <- ukb_query_dictionary("收缩压", official_dict = official, zh_dict = zh, require_rap = FALSE)
  expect_true(nrow(cn$chinese_matches) > 0)
  expect_true(any(cn$official_matches$FieldID == "4080"))
  expect_true(all(c("match_method", "match_score", "review_status") %in% names(cn$official_matches)))

  id <- ukb_query_dictionary("21022", official_dict = official, zh_dict = zh, require_rap = FALSE)
  expect_true(any(id$official_matches$UKBColumn == "p21022"))

  col <- ukb_query_dictionary("participant.p31", official_dict = official, zh_dict = zh, require_rap = FALSE)
  expect_true(any(col$official_matches$FieldID == "31"))
})

test_that("ukb_validate_columns reports missing and matched columns", {
  dat <- data.frame(eid = 1:3, p31 = c(0, 1, 0), participant.p21022 = 50:52)
  out <- ukb_validate_columns(dat, c("eid", "participant.p31", "p21022", "p999"))

  expect_s3_class(out, "ukb_column_validation")
  expect_equal(out$present, c(TRUE, TRUE, TRUE, FALSE))
  expect_false(isTRUE(attr(out, "all_present")))
  expect_error(
    ukb_validate_columns(dat, c("p31", "p999"), error = TRUE),
    "Missing column"
  )
})
