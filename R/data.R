#' Chinese UK Biobank field-path dictionary
#'
#' @description
#' A field-path dictionary used by \code{\link{ukb_query_dictionary}} to support
#' Chinese-language lookup of UK Biobank variables. The table stores a six-level
#' translated category hierarchy and variable label. It does not contain
#' participant-level records and does not include official RAP data values.
#' Official UKB field IDs and RAP column names should still be resolved against
#' a project-specific RAP data dictionary generated inside RAP.
#'
#' @format A data frame with 34,953 rows and 6 translated hierarchy columns.
#' The original UTF-8 column names are preserved in the data object and can be
#' inspected with \code{names(ukb_dictionary_zh)} after loading the dataset.
#'
#' @source Curated UKBAnalytica Chinese field-path dictionary for metadata
#' lookup. This dataset contains metadata labels only.
"ukb_dictionary_zh"
