if (!exists("ukb_metadata_setup", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "variable_preprocess.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "rap_extract.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "field_metadata.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "ukb_metadata.R"), local = FALSE)
}

if (!exists("make_rap_fields", mode = "function")) {
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
}

make_metadata_dict <- function(path) {
  lines <- c(
    "FieldID\tField\tCategory\tValueType\tUnits\tCoding\tInstances\tArray\tNotes",
    "31\tSex\tDemographics\tCategorical\t\t1\t1\t1\tBiological sex",
    "21022\tAge at recruitment\tDemographics\tContinuous\tyears\t\t1\t1\tBaseline age",
    "4080\tSystolic blood pressure, automated reading\tPhysical measures\tInteger\tmmHg\t\t1\t2\tAutomated blood pressure",
    "20116\tSmoking status\tLifestyle\tCategorical\t\t90\t1\t1\tSmoking status"
  )
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

make_codings_file <- function(path) {
  lines <- c(
    "Coding\tValue\tMeaning",
    "1\t0\tFemale",
    "1\t1\tMale",
    "90\t0\tNever",
    "90\t1\tPrevious",
    "90\t2\tCurrent"
  )
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

test_that("ukb_metadata_setup builds unified metadata", {
  dict_path <- tempfile(fileext = ".tsv")
  coding_path <- tempfile(fileext = ".tsv")
  make_metadata_dict(dict_path)
  make_codings_file(coding_path)

  meta <- ukb_metadata_setup(
    source = "files",
    data_dict = dict_path,
    codings = coding_path,
    fields_df = make_rap_fields()
  )

  expect_s3_class(meta, "ukb_metadata")
  expect_true(all(c("fields", "rap_fields", "codings") %in% names(meta)))
  expect_equal(nrow(meta$rap_fields[meta$rap_fields$field_id == 4080, ]), 2)
  expect_equal(meta$rap_fields$instance[meta$rap_fields$rap_column == "participant.p4080_i0_a0"], 0L)
  expect_equal(meta$rap_fields$array[meta$rap_fields$rap_column == "participant.p4080_i0_a0"], 0L)
  expect_equal(nrow(meta$codings[meta$codings$coding_id == "1", ]), 2)
})

test_that("ukb_search_fields and ukb_field_info support simple lookup", {
  dict_path <- tempfile(fileext = ".tsv")
  coding_path <- tempfile(fileext = ".tsv")
  make_metadata_dict(dict_path)
  make_codings_file(coding_path)
  meta <- ukb_metadata_setup(
    source = "files",
    data_dict = dict_path,
    codings = coding_path,
    fields_df = make_rap_fields()
  )

  res <- ukb_search_fields("blood pressure", metadata = meta)
  expect_s3_class(res, "ukb_search_result")
  expect_equal(res$field_id[[1]], 4080L)

  info <- ukb_field_info(4080, metadata = meta)
  expect_s3_class(info, "ukb_field_info")
  expect_equal(info$field$field_id[[1]], 4080L)
  expect_equal(nrow(info$rap_fields), 2)

  sex_info <- ukb_field_info("participant.p31", metadata = meta)
  expect_equal(sex_info$field$field_id[[1]], 31L)
  expect_equal(nrow(sex_info$codings), 2)
})

test_that("ukb_decode decodes values and column names", {
  dict_path <- tempfile(fileext = ".tsv")
  coding_path <- tempfile(fileext = ".tsv")
  make_metadata_dict(dict_path)
  make_codings_file(coding_path)
  meta <- ukb_metadata_setup(
    source = "files",
    data_dict = dict_path,
    codings = coding_path,
    fields_df = make_rap_fields()
  )

  dt <- data.frame(
    eid = 1:3,
    participant.p31 = c(0, 1, -1),
    participant.p4080_i0_a0 = c(120, 130, 125),
    check.names = FALSE
  )

  decoded_values <- ukb_decode_values(dt, metadata = meta, keep_raw = TRUE)
  expect_true("participant.p31_label" %in% names(decoded_values))
  expect_equal(decoded_values$participant.p31_label[1:2], c("Female", "Male"))

  decoded_names <- ukb_decode_column_names(decoded_values, metadata = meta)
  expect_true("sex" %in% names(decoded_names))
  expect_true("sex_label" %in% names(decoded_names))
  expect_true("systolic_blood_pressure_automated_reading_i0_a0" %in% names(decoded_names))
})
