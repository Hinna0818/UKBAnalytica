if (!exists("get_field_metadata", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "variable_preprocess.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "rap_extract.R"), local = FALSE)
  source(testthat::test_path("..", "..", "R", "field_metadata.R"), local = FALSE)
}

mock_live_field_html <- function() {
  paste(
    '<html><body>',
    '<h1>Data-Field 4080</h1>',
    '<table>',
    '<tr><td>Description:</td><td>Systolic blood pressure, automated reading</td></tr>',
    '<tr><td>Category:</td><td>Assessment centre > Physical measures > Blood pressure</td></tr>',
    '<tr><td>Participants</td><td>476,030</td><td>Value Type</td><td>Integer, mmHg</td><td>Sexed</td><td>Both sexes</td><td>Debut</td><td>Jan 2012</td></tr>',
    '<tr><td>Item count</td><td>1,155,218</td><td>Item Type</td><td>Data</td><td>Instances</td><td>Defined (4)</td><td>Version</td><td>Aug 2025</td></tr>',
    '<tr><td>Stability</td><td>Complete</td><td>Strata</td><td>Primary</td><td>Array</td><td>Yes (2)</td><td>Cost Tier</td><td>o1 s1</td></tr>',
    '</table>',
    '<p>Defined-instances run from 0 to 3. Array indices run from 0 to 1. Units of measurement are mmHg.</p>',
    '</body></html>',
    sep = ''
  )
}

make_showcase_dict <- function(path) {
  lines <- c(
    "FieldID\tField\tCategory\tValueType\tUnits\tCoding\tInstances\tArray\tNotes",
    "31\tSex\tDemographics\tCategorical\t\t1\t1\t1\tBiological sex",
    "21022\tAge at recruitment\tDemographics\tContinuous\tyears\t\t1\t1\tBaseline age",
    "4080\tSystolic blood pressure, automated reading\tPhysical measures\tInteger\tmmHg\t\t1\t2\tAutomated blood pressure"
  )
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

test_that("get_field_metadata structures RAP field metadata", {
  meta <- get_field_metadata(fields_df = make_rap_fields())

  expect_true(all(c(
    "field_id", "title", "field_name", "rap_field_names",
    "n_rap_columns", "source_backend"
  ) %in% names(meta)))

  row_4080 <- meta[meta$field_id == 4080, , drop = FALSE]
  expect_equal(nrow(row_4080), 1)
  expect_equal(row_4080$title[[1]], "Systolic blood pressure, automated reading")
  expect_equal(row_4080$n_rap_columns[[1]], 2L)
  expect_match(row_4080$rap_field_names[[1]], "participant.p4080_i0_a0", fixed = TRUE)
})

test_that("get_field_metadata reads showcase metadata and merges RAP fields", {
  dict_path <- tempfile(fileext = ".tsv")
  make_showcase_dict(dict_path)

  meta <- get_field_metadata(
    field_id = c(31, 4080),
    ukb_data_dict = dict_path,
    fields_df = make_rap_fields()
  )

  expect_equal(meta$field_id, c(31L, 4080L))
  expect_equal(meta$units[[2]], "mmHg")
  expect_equal(meta$category[[2]], "Physical measures")
  expect_equal(meta$n_rap_columns[[2]], 2L)
  expect_equal(meta$source_backend[[2]], "showcase+rap")
})

test_that("get_field_metadata supports keyword filtering", {
  dict_path <- tempfile(fileext = ".tsv")
  make_showcase_dict(dict_path)

  meta <- get_field_metadata(
    query = "blood pressure",
    ukb_data_dict = dict_path
  )

  expect_equal(nrow(meta), 1)
  expect_equal(meta$field_id[[1]], 4080L)
  expect_equal(meta$title[[1]], "Systolic blood pressure, automated reading")
})

test_that("get_field_info returns one row for a single field id", {
  dict_path <- tempfile(fileext = ".tsv")
  make_showcase_dict(dict_path)

  info <- get_field_info(
    4080,
    ukb_data_dict = dict_path,
    fields_df = make_rap_fields()
  )

  expect_equal(nrow(info), 1)
  expect_equal(info$field_id[[1]], 4080L)
  expect_equal(info$units[[1]], "mmHg")
  expect_equal(info$n_rap_columns[[1]], 2L)
})

test_that("live field page parser structures UKB showcase field metadata", {
  testthat::skip_if_not_installed("xml2")

  info <- .field_metadata_parse_live_html(
    mock_live_field_html(),
    field_id = 4080,
    url = "https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=4080"
  )

  expect_equal(info$field_id[[1]], 4080L)
  expect_equal(info$title[[1]], "Systolic blood pressure, automated reading")
  expect_equal(info$category[[1]], "Assessment centre > Physical measures > Blood pressure")
  expect_equal(info$value_type[[1]], "Integer")
  expect_equal(info$units[[1]], "mmHg")
  expect_equal(info$participants[[1]], "476,030")
  expect_equal(info$item_count[[1]], "1,155,218")
  expect_equal(info$instances[[1]], "Defined (4)")
  expect_equal(info$array[[1]], "Yes (2)")
  expect_equal(info$version[[1]], "Aug 2025")
  expect_equal(info$source_backend[[1]], "showcase_live")
})