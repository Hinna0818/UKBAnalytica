#' Generate a small synthetic UK Biobank-style demo dataset
#'
#' @description
#' Generates a small fully synthetic toy dataset for documentation and
#' smoke-test workflows. The data are created at runtime from parametric and
#' categorical toy distributions and are not stored in the package as
#' participant-level records.
#'
#' @param n Optional number of rows to return. If `NULL`, 500 rows are returned.
#' @param seed Integer random seed used to generate the toy data. The default
#'   provides reproducible examples. Use `NULL` to avoid setting a seed.
#'
#' @return A data.frame of synthetic cohort variables with missing values
#'   retained.
#' @export
#'
#' @examples
#' demo <- ukb_demo(100)
#' demo2 <- ukb_demo(100, seed = 1)
#' dim(demo)
#' names(demo)
ukb_demo <- function(n = NULL, seed = 20260618L) {
  if (is.null(n)) {
    n <- 500L
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("`n` must be a positive integer or NULL.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single numeric value or NULL.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  cols <- .ukb_demo_colnames()
  baseline <- as.Date("2010-01-01") + sample.int(365 * 3, n, replace = TRUE)
  age <- round(stats::rnorm(n, mean = 58, sd = 8), 1)
  sex <- sample(c(0L, 1L), n, replace = TRUE)
  bmi <- round(stats::rnorm(n, mean = 27, sd = 4.5), 1)
  ethnicity <- sample(c(1L, 1001L, 2001L, 3001L, 4001L), n, replace = TRUE,
                      prob = c(0.82, 0.06, 0.05, 0.04, 0.03))
  education <- sample(c(1L, 2L, 3L, 4L, 5L, 6L, -3L), n, replace = TRUE,
                      prob = c(0.28, 0.18, 0.25, 0.08, 0.11, 0.08, 0.02))
  smoking <- sample(c(0L, 1L, 2L, -3L), n, replace = TRUE,
                    prob = c(0.55, 0.32, 0.11, 0.02))
  drinking <- sample(c(0L, 1L, 2L, -3L), n, replace = TRUE,
                    prob = c(0.10, 0.62, 0.26, 0.02))
  tdi <- round(stats::rnorm(n, mean = 0, sd = 3), 2)

  no2_base <- stats::rnorm(n, mean = 25, sd = 6)
  nox_base <- stats::rnorm(n, mean = 42, sd = 12)
  pm25_base <- stats::rnorm(n, mean = 10, sd = 2)
  pm10_base <- stats::rnorm(n, mean = 18, sd = 4)

  out <- as.data.frame(
    setNames(replicate(length(cols), rep(NA_real_, n), simplify = FALSE), cols),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out$eid <- 90000L + seq_len(n)

  if ("p31" %in% cols) out$p31 <- sex
  if ("p21022" %in% cols) out$p21022 <- age
  for (col in grep("^p53_i[0-9]+$", cols, value = TRUE)) {
    instance <- as.integer(sub("^p53_i", "", col))
    out[[col]] <- baseline + round(instance * 365.25 * 2)
  }
  for (col in grep("^p21001_i[0-9]+$", cols, value = TRUE)) out[[col]] <- bmi + stats::rnorm(n, 0, 0.8)
  for (col in grep("^p21000_i[0-9]+$", cols, value = TRUE)) out[[col]] <- ethnicity
  for (col in grep("^p20116_i[0-9]+$", cols, value = TRUE)) out[[col]] <- smoking
  for (col in grep("^p20117_i[0-9]+$", cols, value = TRUE)) out[[col]] <- drinking
  for (col in grep("^p6138_i[0-9]+$", cols, value = TRUE)) out[[col]] <- education
  for (col in grep("^p20161_i[0-9]+$", cols, value = TRUE)) {
    out[[col]] <- sample(c(0L, 1L), n, replace = TRUE, prob = c(0.72, 0.28))
  }
  if ("p189" %in% cols) out$p189 <- tdi

  .assign_if_present <- function(col, value) {
    if (col %in% cols) out[[col]] <<- value
  }
  .assign_if_present("p24016", round(pmax(no2_base + stats::rnorm(n, 0, 1.2), 1), 2))
  .assign_if_present("p24018", round(pmax(no2_base + stats::rnorm(n, 0, 1.2), 1), 2))
  .assign_if_present("p24017", round(pmax(no2_base + stats::rnorm(n, 0, 1.2), 1), 2))
  .assign_if_present("p24003", round(pmax(no2_base + stats::rnorm(n, 0, 1.8), 1), 2))
  .assign_if_present("p24004", round(pmax(nox_base, 1), 2))
  .assign_if_present("p24006", round(pmax(pm25_base, 1), 2))
  .assign_if_present("p24019", round(pmax(pm10_base + stats::rnorm(n, 0, 1.0), 1), 2))
  .assign_if_present("p24005", round(pmax(pm10_base + stats::rnorm(n, 0, 1.4), 1), 2))

  for (col in grep("^p22009_a[0-9]+$", cols, value = TRUE)) {
    pc <- as.integer(sub("^p22009_a", "", col))
    out[[col]] <- stats::rnorm(n, mean = 0, sd = 1 / sqrt(pc + 1))
  }
  for (col in grep("^p3062_", cols, value = TRUE)) out[[col]] <- round(stats::rnorm(n, 2.7, 0.7), 2)
  for (col in grep("^p3063_", cols, value = TRUE)) out[[col]] <- round(stats::rnorm(n, 3.6, 0.8), 2)
  for (col in grep("^p3064_", cols, value = TRUE)) out[[col]] <- round(stats::rnorm(n, 78, 14), 1)

  icd10_codes <- rep("", n)
  icd10_dates <- as.Date(rep(NA, n), origin = "1970-01-01")
  if (n >= 10L) {
    prevalent_n <- max(1L, floor(0.04 * n))
    incident_n <- max(1L, floor(0.08 * n))
    prevalent_idx <- seq_len(prevalent_n)
    incident_idx <- seq(from = prevalent_n + 1L, length.out = incident_n)
    incident_idx <- incident_idx[incident_idx <= n]
    icd10_codes[prevalent_idx] <- "['I48']"
    icd10_dates[prevalent_idx] <- baseline[prevalent_idx] - sample(30:1200, length(prevalent_idx), replace = TRUE)
    icd10_codes[incident_idx] <- "['I48']"
    icd10_dates[incident_idx] <- baseline[incident_idx] + sample(30:3000, length(incident_idx), replace = TRUE)
    copd_idx <- sample(setdiff(seq_len(n), c(prevalent_idx, incident_idx)), max(1L, floor(0.03 * n)))
    icd10_codes[copd_idx] <- "['J44']"
    icd10_dates[copd_idx] <- baseline[copd_idx] + sample(30:2500, length(copd_idx), replace = TRUE)
  }
  .assign_if_present("p41270", icd10_codes)
  icd10_date_cols <- grep("^p41280_a[0-9]+$", cols, value = TRUE)
  if (length(icd10_date_cols) > 0L) out[[icd10_date_cols[[1L]]]] <- icd10_dates

  .assign_if_present("p41271", sample(c("", "4273", "496"), n, replace = TRUE, prob = c(0.92, 0.05, 0.03)))
  icd9_date_cols <- grep("^p41281_a[0-9]+$", cols, value = TRUE)
  if (length(icd9_date_cols) > 0L) {
    has_icd9 <- nzchar(out$p41271)
    out[[icd9_date_cols[[1L]]]] <- as.Date(rep(NA, n), origin = "1970-01-01")
    out[[icd9_date_cols[[1L]]]][has_icd9] <- baseline[has_icd9] + sample(-800:2200, sum(has_icd9), replace = TRUE)
  }

  sr_code_cols <- grep("^p20002_i[0-9]+_a[0-9]+$", cols, value = TRUE)
  sr_date_cols <- grep("^p20008_i[0-9]+_a[0-9]+$", cols, value = TRUE)
  if (length(sr_code_cols) > 0L) {
    sr_codes <- c(1065L, 1066L, 1074L, 1111L, 1112L, 1113L, 1220L, 1222L, 1472L)
    out[[sr_code_cols[[1L]]]] <- sample(c(sr_codes, NA_integer_), n, replace = TRUE,
                                        prob = c(rep(0.04, length(sr_codes)), 0.64))
    if (length(sr_date_cols) > 0L) {
      has_sr <- !is.na(out[[sr_code_cols[[1L]]]])
      out[[sr_date_cols[[1L]]]] <- as.Date(rep(NA, n), origin = "1970-01-01")
      out[[sr_date_cols[[1L]]]][has_sr] <- baseline[has_sr] - sample(30:2500, sum(has_sr), replace = TRUE)
    }
  }

  for (col in grep("^p40001_i[0-9]+$", cols, value = TRUE)) out[[col]] <- sample(c(0L, 1L, 2L, NA_integer_), n, TRUE, c(0.94, 0.03, 0.01, 0.02))
  cancer_code_cols <- grep("^p40002_i[0-9]+_a[0-9]+$", cols, value = TRUE)
  if (length(cancer_code_cols) > 0L) {
    out[[cancer_code_cols[[1L]]]] <- sample(c("", "C34", "C50", "C61", "C18"), n, TRUE, c(0.94, 0.015, 0.015, 0.015, 0.015))
  }
  for (col in grep("^p40000_i[0-9]+$", cols, value = TRUE)) {
    has_cancer <- if (length(cancer_code_cols) > 0L) nzchar(out[[cancer_code_cols[[1L]]]]) else rep(FALSE, n)
    out[[col]] <- as.Date(rep(NA, n), origin = "1970-01-01")
    out[[col]][has_cancer] <- baseline[has_cancer] + sample(-1200:2500, sum(has_cancer), replace = TRUE)
  }

  missing_cols <- intersect(
    c(grep("^p21001_", cols, value = TRUE), grep("^p20116_", cols, value = TRUE),
      grep("^p20117_", cols, value = TRUE), grep("^p21000_", cols, value = TRUE),
      grep("^p6138_", cols, value = TRUE), grep("^p240", cols, value = TRUE),
      grep("^p306", cols, value = TRUE), "p189"),
    cols
  )
  for (col in missing_cols) {
    miss_n <- floor(0.03 * n)
    if (miss_n > 0L) {
      out[sample.int(n, miss_n), col] <- NA
    }
  }

  message("This dataset is generated at runtime for documentation and testing. It does not contain UK Biobank participant data.")
  out[, cols, drop = FALSE]
}

#' Get column names of the synthetic UK Biobank-style demo dataset
#'
#' @description
#' Returns the column names generated by `ukb_demo()`. This is useful for
#' documentation examples that need RAP-style toy column names.
#'
#' @return A character vector of original demo-data column names.
#' @export
#'
#' @examples
#' get_ukb_demo_colnames()
get_ukb_demo_colnames <- function() {
  .ukb_demo_colnames()
}

.ukb_demo_colnames <- function() {
  if (exists("ukb_demo_colnames", inherits = TRUE)) {
    return(get("ukb_demo_colnames", inherits = TRUE))
  }
  c("eid", "p31", "p53_i0", "p21022", "p21001_i0", "p41270", "p41280_a0")
}
