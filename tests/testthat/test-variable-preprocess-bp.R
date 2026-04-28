if (!exists("calculate_blood_pressure", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "variable_preprocess.R"), local = FALSE)
}

test_that("calculate_blood_pressure defaults to automated readings with manual fallback", {
  dt <- data.frame(
    p4080_i0_a0 = c(120, NA, NA),
    p4080_i0_a1 = c(122, 118, NA),
    p93_i0_a0 = c(140, 130, 150),
    p93_i0_a1 = c(142, 132, 152),
    p4079_i0_a0 = c(80, NA, NA),
    p4079_i0_a1 = c(82, 78, NA),
    p94_i0_a0 = c(90, 88, 92),
    p94_i0_a1 = c(92, 90, 94)
  )

  sbp_res <- calculate_blood_pressure(dt, type = "sbp")
  dbp_res <- calculate_blood_pressure(dt, type = "dbp")

  expect_equal(sbp_res$sbp, c(121, 124, 151))
  expect_equal(dbp_res$dbp, c(81, 83, 93))
})

test_that("calculate_blood_pressure supports manual preference", {
  dt <- data.frame(
    p4080_i0_a0 = c(120, 126),
    p4080_i0_a1 = c(122, NA),
    p93_i0_a0 = c(140, NA),
    p93_i0_a1 = c(142, 132)
  )

  res <- calculate_blood_pressure(dt, type = "sbp", prefer = "manual")

  expect_equal(res$sbp, c(141, 129))
})

test_that("blood pressure field mapping matches UKB field definitions", {
  bp_info <- get_variable_info("blood_pressure")

  field_lookup <- setNames(bp_info$field_id, bp_info$variable)
  desc_lookup <- setNames(bp_info$description, bp_info$variable)

  expect_equal(unname(field_lookup["sbp_auto_1"]), 4080)
  expect_equal(unname(field_lookup["sbp_manual_1"]), 93)
  expect_equal(unname(field_lookup["dbp_auto_1"]), 4079)
  expect_equal(unname(field_lookup["dbp_manual_1"]), 94)

  expect_match(unname(desc_lookup["sbp_auto_1"]), "automated")
  expect_match(unname(desc_lookup["sbp_manual_1"]), "manual")
  expect_match(unname(desc_lookup["dbp_auto_1"]), "automated")
  expect_match(unname(desc_lookup["dbp_manual_1"]), "manual")
})
