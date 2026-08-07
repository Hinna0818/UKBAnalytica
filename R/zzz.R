.onAttach <- function(libname, pkgname) {
  version <- tryCatch(
    as.character(packageVersion(pkgname)),
    error = function(e) NA_character_
  )
  version_label <- if (!is.na(version)) paste0(" ", version) else ""

  packageStartupMessage(
    paste0(
      "UKBAnalytica", version_label,
      ": RAP-oriented UK Biobank workflows. ",
      "Run ukb_check_rap_env() before RAP data extraction."
    )
  )
}
