# Internal utility: coerce mixed date inputs to Date without throwing hard errors.
.safe_as_date <- function(x, col_name = NULL, warn = TRUE) {
  if (inherits(x, "Date")) {
    return(x)
  }

  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  n <- length(x)
  out <- as.Date(rep(NA_character_, n))
  if (n == 0) {
    return(out)
  }

  orig_x <- x

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.numeric(x)) {
    x_num <- suppressWarnings(as.numeric(x))
    is_missing <- is.na(x_num)

    # UKB-like compact dates such as 20170131.
    x_int <- suppressWarnings(as.integer(abs(trunc(x_num))))
    compact_idx <- !is_missing & x_int >= 19000101L & x_int <= 21001231L
    if (any(compact_idx)) {
      compact_chr <- sprintf("%08d", x_int[compact_idx])
      out[compact_idx] <- as.Date(compact_chr, format = "%Y%m%d")
    }

    # Compact year-month values such as 201701; impute day=15.
    ym_compact_idx <- !is_missing & !compact_idx & x_int >= 190001L & x_int <= 210012L
    if (any(ym_compact_idx)) {
      ym_chr <- sprintf("%06d", x_int[ym_compact_idx])
      ym_mid <- paste0(substr(ym_chr, 1, 4), "-", substr(ym_chr, 5, 6), "-15")
      out[ym_compact_idx] <- as.Date(ym_mid, format = "%Y-%m-%d")
    }

    # R Date numeric origin fallback.
    # Keep this conservative to avoid interpreting year-like values (e.g. 2012)
    # as day offsets since 1970.
    origin_idx <- !is_missing & !compact_idx & !ym_compact_idx & x_num >= 5000 & x_num <= 100000
    if (any(origin_idx)) {
      out[origin_idx] <- as.Date(x_num[origin_idx], origin = "1970-01-01")
    }
  } else {
    x_chr <- trimws(as.character(x))
    x_chr[x_chr %in% c("", "NA", "N/A", "NULL", "null", "NaN")] <- NA_character_
    x_norm <- gsub("[./]", "-", x_chr)

    ymd_idx <- !is.na(x_norm) & grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", x_norm)
    if (any(ymd_idx)) {
      out[ymd_idx] <- as.Date(x_norm[ymd_idx], format = "%Y-%m-%d")
    }

    dmy_idx <- is.na(out) & !is.na(x_norm) & grepl("^\\d{1,2}-\\d{1,2}-\\d{4}$", x_norm)
    if (any(dmy_idx)) {
      out[dmy_idx] <- as.Date(x_norm[dmy_idx], format = "%d-%m-%Y")
    }

    compact_char_idx <- is.na(out) & !is.na(x_norm) & grepl("^\\d{8}$", x_norm)
    if (any(compact_char_idx)) {
      out[compact_char_idx] <- as.Date(x_norm[compact_char_idx], format = "%Y%m%d")
    }

    # Partial year-month values are common in some exports; impute day=15.
    ym_idx <- is.na(out) & !is.na(x_norm) & grepl("^\\d{4}-\\d{1,2}$", x_norm)
    if (any(ym_idx)) {
      ym_mid <- paste0(x_norm[ym_idx], "-15")
      out[ym_idx] <- as.Date(ym_mid, format = "%Y-%m-%d")
    }
  }

  if (isTRUE(warn)) {
    orig_chr <- trimws(as.character(orig_x))
    orig_missing <- is.na(orig_x) | orig_chr %in% c("", "NA", "N/A", "NULL", "null", "NaN")
    introduced_na <- !orig_missing & is.na(out)
    if (any(introduced_na)) {
      label <- if (!is.null(col_name) && nzchar(col_name)) col_name else "date input"
      bad_values <- unique(as.character(orig_x[introduced_na]))
      bad_preview <- paste(utils::head(bad_values, 3L), collapse = ", ")
      warning(
        sprintf(
          "[.safe_as_date] %d non-standard value(s) in '%s' were set to NA. Example(s): %s",
          sum(introduced_na),
          label,
          bad_preview
        ),
        call. = FALSE
      )
    }
  }

  out
}
