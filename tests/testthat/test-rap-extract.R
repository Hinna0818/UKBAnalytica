if (!exists("rap_plan_extract", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "variable_preprocess.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "rap_extract.R"), local = FALSE)
}

make_rap_fields <- function() {
  data.frame(
    field_name = c(
      "participant.eid",
      "participant.p31",
      "participant.p53_i0",
      "participant.p21022",
      "participant.p21001_i0",
      "participant.p4080_i0_a0",
      "participant.p4080_i0_a1",
      "participant.p93_i0_a0",
      "participant.p20002_i0_a0"
    ),
    title = c(
      "Participant ID",
      "Sex",
      "Date of attending assessment centre | Instance 0",
      "Age at recruitment",
      "Body mass index | Instance 0",
      "Systolic blood pressure, automated reading | Instance 0 | Array 0",
      "Systolic blood pressure, automated reading | Instance 0 | Array 1",
      "Systolic blood pressure, manual reading | Instance 0 | Array 0",
      "Non-cancer illness code, self-reported | Instance 0 | Array 0"
    ),
    stringsAsFactors = FALSE
  )
}

test_that("rap field parser reads tab-separated dx output", {
  raw <- paste(
    "participant.p31\tSex",
    "participant.p53_i0\tDate of attending assessment centre | Instance 0",
    sep = "\n"
  )

  res <- .rap_parse_fields(raw)

  expect_equal(names(res), c("field_name", "title"))
  expect_equal(res$field_name, c("participant.p31", "participant.p53_i0"))
  expect_equal(res$title[[1]], "Sex")
})

test_that("rap_plan_extract matches field IDs without prefix collisions", {
  fields <- make_rap_fields()

  plan <- rap_plan_extract(
    field_id = c(31, 53, 4080),
    dataset = "app123.dataset",
    fields_df = fields
  )

  expect_s3_class(plan, "rap_extract_plan")
  expect_true("participant.eid" %in% plan$fields)
  expect_true(all(c("participant.p31", "participant.p53_i0") %in% plan$fields))
  expect_true(all(c("participant.p4080_i0_a0", "participant.p4080_i0_a1") %in% plan$fields))
  expect_false(any(grepl("^participant\\.p93", plan$fields)))
})

test_that("rap_plan_extract supports predefined UKBAnalytica variable names", {
  fields <- make_rap_fields()

  plan <- rap_plan_extract(
    variables = c("sex", "age", "bmi", "sbp_auto_1"),
    dataset = "app123.dataset",
    fields_df = fields
  )

  expect_true(all(c(
    "participant.p31",
    "participant.p21022",
    "participant.p21001_i0",
    "participant.p4080_i0_a0"
  ) %in% plan$fields))
  expect_false("participant.p4080_i0_a1" %in% plan$fields)
})

test_that("rap_submit_extract dry-run uses table-exporter field format", {
  fields <- make_rap_fields()

  plan <- rap_submit_extract(
    field_id = c(31, 4080),
    dataset = "app123.dataset",
    fields_df = fields,
    file = "baseline_core.csv",
    dry_run = TRUE
  )

  expect_s3_class(plan, "rap_extract_plan")
  expect_equal(plan$output, "baseline_core")
  expect_equal(plan$instance_type, "mem1_ssd1_v2_x4")
  expect_true(all(c("eid", "p31", "p4080_i0_a0", "p4080_i0_a1") %in% plan$fields))
  expect_false(any(grepl("^participant\\.", plan$fields)))
  expect_true(any(grepl("table-exporter", plan$command, fixed = TRUE)))
})
