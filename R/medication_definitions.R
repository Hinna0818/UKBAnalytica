#' Create a medication definition object
#'
#' @description
#' Helper for defining medication code sets from UK Biobank self-reported
#' treatment/medication fields. The first implementation focuses on field
#' 20003 arrays and intentionally stores only medication codes and classes,
#' not copied source codelist descriptions.
#'
#' @param name Medication definition name.
#' @param codes Character or numeric medication codes.
#' @param source Source label. Defaults to `"Self-report 20003"`.
#' @param field_id UK Biobank field ID. Defaults to 20003.
#' @param medication_class Optional medication class label.
#' @param match_type Matching mode. Defaults to `"exact"`.
#'
#' @return A list describing the medication definition.
#' @export
#'
#' @examples
#' bp <- create_medication_definition("Any BP medication", c(1, 2, 3))
create_medication_definition <- function(name,
                                         codes,
                                         source = "Self-report 20003",
                                         field_id = 20003L,
                                         medication_class = NULL,
                                         match_type = "exact") {
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }
  codes <- as.character(codes)
  codes <- unique(codes[!is.na(codes) & nzchar(codes)])
  if (length(codes) == 0) {
    stop("`codes` must contain at least one non-missing code.", call. = FALSE)
  }
  if (!is.character(source) || length(source) != 1 || is.na(source) || !nzchar(source)) {
    stop("`source` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.numeric(field_id) || length(field_id) != 1 || is.na(field_id)) {
    stop("`field_id` must be a single numeric field ID.", call. = FALSE)
  }
  if (!is.character(match_type) || length(match_type) != 1 || is.na(match_type) || !nzchar(match_type)) {
    stop("`match_type` must be a single non-empty character string.", call. = FALSE)
  }

  list(
    name = name,
    codes = codes,
    source = source,
    field_id = as.integer(field_id),
    medication_class = medication_class,
    match_type = match_type
  )
}

#' Get predefined UK Biobank medication definitions
#'
#' @description
#' Returns curated field-20003 medication code sets for common self-reported
#' treatment groups. These definitions are designed for baseline covariate
#' derivation and sensitivity analyses, and are separate from disease endpoint
#' definitions returned by [get_predefined_diseases()].
#'
#' @return A named list of medication definition objects.
#' @export
#'
#' @examples
#' meds <- get_predefined_medications()
#' names(meds)
get_predefined_medications <- function() {
  bp <- list(
    ACE_Inhibitor = .ukb_split_codes(
      "1140851690,1140851692,1140860696,1140860706,1140860714,1140860728,1140860736,1140860738,1140860750,1140860752,1140860758,1140860764,1140860776,1140860784,1140860790,1140860802,1140860806,1140860878,1140860882,1140860892,1140860904,1140860912,1140860918,1140864618,1140864910,1140864952,1140881706,1140881712,1140881714,1140881716,1140888552,1140888556,1140888560,1140923712,1140923718,1141150328,1141150560,1141151382,1141153316,1141153328,1141164148,1141164154,1141165470,1141165476,1141167758,1141167822,1141169498,1141170544,1141170870,1141180592,1141180598,1141181186,1141188408,1141190934,1141199940,1141200698,1141200726"
    ),
    ARB = .ukb_split_codes(
      "1140916356,1140916362,1141145660,1141145668,1141151016,1141151018,1141152998,1141153006,1141156836,1141156846,1141166006,1141171336,1141171344,1141172492,1141172682,1141172686,1141179974,1141187788,1141187790,1141193282,1141193346,1141201038,1141201040"
    ),
    Beta_Blocker = .ukb_split_codes(
      "1140851480,1140851484,1140851492,1140851508,1140851522,1140851556,1140851576,1140860172,1140860180,1140860192,1140860194,1140860212,1140860220,1140860222,1140860230,1140860232,1140860244,1140860250,1140860266,1140860274,1140860278,1140860292,1140860294,1140860304,1140860308,1140860312,1140860314,1140860316,1140860318,1140860320,1140860322,1140860324,1140860328,1140860330,1140860332,1140860334,1140860336,1140860338,1140860340,1140860342,1140860348,1140860352,1140860356,1140860358,1140860362,1140860380,1140860382,1140860386,1140860390,1140860394,1140860396,1140860398,1140860400,1140860402,1140860404,1140860406,1140860410,1140860418,1140860422,1140860426,1140860434,1140860492,1140860498,1140861456,1140863724,1140864176,1140864410,1140864950,1140866704,1140866712,1140866724,1140866726,1140866738,1140866756,1140866758,1140866764,1140866766,1140866778,1140866782,1140866784,1140866798,1140866800,1140866802,1140866804,1140875808,1140879758,1140879760,1140879762,1140879818,1140879822,1140879824,1140879830,1140879834,1140879842,1140879854,1140879866,1140881722,1140883554,1140909368,1140910614,1140916342,1140916628,1140916730,1140916868,1140917076,1140922930,1140923336,1140923404,1141146124,1141146126,1141146128,1141146184,1141152076,1141156754,1141156808,1141164276,1141164280,1141168498,1141169516,1141171152,1141172742,1141180778,1141182904,1141182968,1141184324,1141184722,1141187048,1141187780,1141194804,1141194808,1141194810"
    ),
    Calcium_Channel_Blocker = .ukb_split_codes(
      "1140851730,1140851790,1140851794,1140851798,1140851800,1140861088,1140861090,1140861106,1140861110,1140861114,1140861120,1140861128,1140861130,1140861136,1140861138,1140861166,1140861176,1140861190,1140861194,1140861202,1140861276,1140861282,1140862728,1140866460,1140866466,1140866484,1140866546,1140866554,1140868036,1140872472,1140872568,1140879802,1140879806,1140879810,1140881692,1140881702,1140888510,1140888646,1140911088,1140911698,1140916930,1140917428,1140917452,1140923572,1140923618,1140926188,1140926778,1140926780,1140926954,1140926966,1140927934,1140927940,1140928212,1140928226,1140928234,1141145870,1141150500,1141150538,1141150926,1141151474,1141152600,1141153026,1141153032,1141153394,1141153400,1141153454,1141156656,1141157136,1141157140,1141162546,1141166752,1141167832,1141169096,1141169710,1141169730,1141171804,1141173766,1141174684,1141175224,1141180238,1141184390,1141185444,1141187056,1141187094,1141187774,1141187962,1141188152,1141188576,1141188836,1141188920,1141188936,1141190160,1141190548,1141199858,1141200400,1141200782,1141201814"
    ),
    Thiazide_Diuretic = .ukb_split_codes(
      "1140851332,1140851336,1140851338,1140851362,1140851364,1140851368,1140851430,1140851432,1140851436,1140851660,1140860562,1140864202,1140866072,1140866074,1140866078,1140866090,1140866092,1140866094,1140866096,1140866102,1140866104,1140866122,1140866128,1140866132,1140866136,1140866138,1140866140,1140866144,1140866146,1140866156,1140866158,1140866162,1140866164,1140866168,1140866324,1140866328,1140866330,1140866340,1140866352,1140866354,1140866360,1140866396,1140866400,1140866402,1140866404,1140866410,1140866416,1140866420,1140866422,1140866440,1140866446,1140866450,1140888918,1140888922,1140909706,1140910442,1140916870,1140917068,1140922324,1140923272,1140923276,1140923282,1141146378,1141180772,1141188636,1141194794,1141194800"
    )
  )

  dm <- list(
    Acarbose = .ukb_split_codes("1140868902,1140868908"),
    Glucagon = .ukb_split_codes("1140874826,1140923026,1141173144"),
    Insulin = .ukb_split_codes("1140883066"),
    Meglitinide = .ukb_split_codes("1141168660,1141168668,1141173786,1141173882"),
    Metformin = .ukb_split_codes("1140874686,1140884600,1140921964,1141189094,1141189090"),
    Sulfonamide = .ukb_split_codes("1140857500,1140857502"),
    Sulfonylurea = .ukb_split_codes(
      "1140857494,1140857496,1140857506,1140857584,1140857586,1140857590,1140874646,1140874650,1140874652,1140874658,1140874660,1140874664,1140874666,1140874674,1140874678,1140874680,1140874690,1140874706,1140874712,1140874716,1140874718,1140874724,1140874726,1140874728,1140874732,1140874736,1140874740,1140874744,1140874746,1141152590,1141156984,1141157284,1141169504,1141171508"
    ),
    TZD = .ukb_split_codes("1141153254,1141153262,1141171646,1141171652,1141177600,1141177606")
  )

  definitions <- list(
    Blood_Pressure_Medication = create_medication_definition(
      "Blood Pressure Medication",
      unique(unlist(bp, use.names = FALSE)),
      medication_class = "antihypertensive"
    ),
    Diabetes_Medication = create_medication_definition(
      "Diabetes Medication",
      unique(unlist(dm, use.names = FALSE)),
      medication_class = "glucose_lowering"
    )
  )

  for (nm in names(bp)) {
    definitions[[nm]] <- create_medication_definition(
      gsub("_", " ", nm),
      bp[[nm]],
      medication_class = tolower(nm)
    )
  }
  for (nm in names(dm)) {
    definitions[[nm]] <- create_medication_definition(
      gsub("_", " ", nm),
      dm[[nm]],
      medication_class = tolower(nm)
    )
  }
  definitions
}

#' Load UK Biobank field 20003 medication coding
#'
#' @description
#' Loads the UK Biobank coding 4 table used by field 20003
#' (treatment/medication code). This table is included as a lightweight
#' reference so users can inspect the meaning of medication codes used by
#' `get_predefined_medications()`.
#'
#' @param path Optional path to a local coding 4 TSV file. If NULL, the package
#'   copy in `inst/extdata/ukb_coding4_20003.tsv` is used.
#'
#' @return A data.frame with columns `coding` and `meaning`.
#' @export
#'
#' @examples
#' coding4 <- load_ukb_medication_coding()
#' head(coding4)
load_ukb_medication_coding <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "ukb_coding4_20003.tsv", package = "UKBAnalytica")
  }
  if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Medication coding file was not found: ", path, call. = FALSE)
  }

  out <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out) <- tolower(names(out))
  if (!all(c("coding", "meaning") %in% names(out))) {
    stop("Medication coding file must contain `coding` and `meaning` columns.", call. = FALSE)
  }
  out$coding <- as.character(out$coding)
  out
}

#' Extract self-reported medication indicators from field 20003
#'
#' @description
#' Matches UK Biobank treatment/medication code arrays (`p20003_i*_a*`) against
#' predefined or user-supplied medication definitions and appends binary
#' participant-level medication indicators.
#'
#' @param data A data.frame or data.table containing field 20003 array columns.
#' @param medications Optional medication definition names to extract. If NULL,
#'   all predefined definitions are used.
#' @param medication_definitions Named list of medication definitions. Defaults
#'   to `get_predefined_medications()`.
#' @param id_col Participant identifier column.
#' @param instance Optional UKB assessment instance. If NULL, all available
#'   instances are searched.
#' @param prefix Prefix for output variable names.
#' @param missing_as_zero Logical. If TRUE, participants with no valid 20003
#'   entries are coded as 0; otherwise they are coded as NA.
#' @param return_long Logical. If TRUE, return one row per participant and
#'   medication definition instead of appending wide columns.
#'
#' @return A data.table.
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   eid = 1:3,
#'   p20003_i0_a0 = c("1140883066", "1140874686", NA),
#'   p20003_i0_a1 = c(NA, "1140851690", NA)
#' )
#' extract_self_report_medications(dat, medications = c("Insulin", "Metformin"))
extract_self_report_medications <- function(data,
                                            medications = NULL,
                                            medication_definitions = get_predefined_medications(),
                                            id_col = "eid",
                                            instance = 0,
                                            prefix = "med20003",
                                            missing_as_zero = TRUE,
                                            return_long = FALSE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!id_col %in% names(data)) {
    stop(sprintf("`id_col` was not found in data: %s", id_col), call. = FALSE)
  }
  if (!is.list(medication_definitions) || length(medication_definitions) == 0) {
    stop("`medication_definitions` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(medication_definitions)) || any(!nzchar(names(medication_definitions)))) {
    stop("`medication_definitions` must be a named list.", call. = FALSE)
  }
  if (!is.null(medications)) {
    missing_defs <- setdiff(medications, names(medication_definitions))
    if (length(missing_defs) > 0) {
      stop(
        "Unknown medication definition(s): ",
        paste(missing_defs, collapse = ", "),
        ". Use names(get_predefined_medications()) to inspect available definitions.",
        call. = FALSE
      )
    }
    medication_definitions <- medication_definitions[medications]
  }

  dt <- data.table::as.data.table(data.table::copy(data))
  code_cols <- grep("^p20003_i[0-9]+_a[0-9]+$", names(dt), value = TRUE)
  if (!is.null(instance)) {
    if (!is.numeric(instance) || length(instance) != 1 || is.na(instance) || instance < 0) {
      stop("`instance` must be NULL or a single non-negative number.", call. = FALSE)
    }
    code_cols <- grep(paste0("^p20003_i", as.integer(instance), "_a[0-9]+$"), code_cols, value = TRUE)
  }
  if (length(code_cols) == 0) {
    warning("No field 20003 medication code columns found.", call. = FALSE)
    return(dt[])
  }

  code_mat <- as.data.frame(dt[, code_cols, with = FALSE], stringsAsFactors = FALSE)
  code_mat[] <- lapply(code_mat, as.character)
  valid_mat <- as.data.frame(
    lapply(code_mat, function(x) !is.na(x) & x != "" & !x %in% c("-1", "-3")),
    stringsAsFactors = FALSE
  )
  has_any_valid <- rowSums(valid_mat) > 0

  long_rows <- list()
  for (nm in names(medication_definitions)) {
    def <- medication_definitions[[nm]]
    if (is.null(def$codes)) {
      stop(sprintf("Medication definition '%s' does not contain `codes`.", nm), call. = FALSE)
    }
    matched <- Reduce(`|`, lapply(code_mat, function(x) x %in% def$codes))
    value <- as.integer(matched)
    if (!isTRUE(missing_as_zero)) {
      value[!has_any_valid] <- NA_integer_
    }
    out_col <- paste(prefix, .ukb_med_snake(nm), sep = "_")
    dt[, (out_col) := value]

    long_rows[[nm]] <- data.table::data.table(
      eid = dt[[id_col]],
      medication = nm,
      medication_name = def$name,
      medication_class = if (is.null(def$medication_class)) NA_character_ else def$medication_class,
      source = def$source,
      field_id = def$field_id,
      has_medication = value
    )
  }

  if (isTRUE(return_long)) {
    return(data.table::rbindlist(long_rows, use.names = TRUE, fill = TRUE))
  }
  dt[]
}

.ukb_split_codes <- function(x) {
  x <- unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE), use.names = FALSE)
  unique(trimws(x[nzchar(x)]))
}

.ukb_med_snake <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_+|_+$", "", x)
}
