# `ukbsci-rap-extract` — examples

Three copy-pastable end-to-end examples. Run them from inside an
**authenticated RAP R session**.

Every example prepends the same header:

```r
###############################################################################
# UKBAnalytica skill: ukbsci-rap-extract
# Citation: He N. UKBAnalytica R package. https://github.com/Hinna0818/UKBAnalytica_v2
###############################################################################

library(UKBAnalytica)
library(data.table)

stopifnot(nzchar(Sys.which("dx")))                  # dx must be installed
ok <- system2("dx", "whoami", stdout = TRUE)
stopifnot(length(ok) > 0)                           # dx must be logged in
```

---

## Example A — Small synchronous extract (demographics + BP)

Use when: ≤ ~500 output columns, ≤ ~1M rows, exploratory.

```r
# Step 1. Discover dataset
dataset <- rap_find_dataset()

# Step 2. Search for the fields you want
meta <- ukb_metadata_setup(source = "auto", dataset = dataset, cache = TRUE)
hits <- ukb_search_fields(query = "blood pressure", metadata = meta,
                          max_results = 30)
print(hits[, c("field_id", "title", "value_type", "instances", "array")])

# Step 3. Plan the extraction (numeric field IDs)
plan <- rap_plan_extract(
  field_id = c(31,     # Sex
               21003,  # Age at recruitment
               21001,  # BMI
               4079,   # Diastolic BP automated
               4080),  # Systolic BP automated
  dataset  = dataset,
  manifest = "/mnt/project/extracts/manifest_demo_bp.csv"
)

stopifnot(length(plan$unmatched) == 0)
cat("Plan will produce", plan$n_columns, "columns.\n")

# Step 4. Run it synchronously, output stays in RAP storage
dt <- rap_extract_pheno(
  field_id            = plan$matched$field_id,
  dataset             = dataset,
  output              = "/mnt/project/extracts/demo_bp.csv",
  read                = TRUE,
  strip_entity_prefix = TRUE
)

# Step 5. Decode columns + coded values
dt <- ukb_decode(dt, metadata = meta,
                 decode_names = TRUE, decode_values = TRUE,
                 keep_raw = TRUE)

str(dt)
```

---

## Example B — Large asynchronous extract (full hospital history)

Use when: full cohort + many array fields (`p41270` / `p41280_a*` for ICD-10,
`p20002_i*_a*` for self-report, etc.).

```r
dataset <- rap_find_dataset()

# Plan first — this is what tells you "go async"
plan <- rap_plan_extract(
  field_id = c(31, 21003, 53,        # baseline scaffolding
               41270, 41280,         # ICD-10 codes + dates
               41271, 41281,         # ICD-9
               41272, 41282,         # OPCS-4 procedures
               40000, 40001, 40002,  # death dates + primary + secondary causes
               20002, 20008),        # self-report illness + year
  dataset = dataset
)
cat("Planned columns:", plan$n_columns, "\n")
if (plan$n_columns > 500) message("Routing to async via rap_submit_extract().")

# Submit the table-exporter job
job <- rap_submit_extract(
  field_id      = plan$matched$field_id,
  dataset       = dataset,
  file          = "pheno_full_baseline_records",
  instance_type = NULL,           # auto-pick
  priority      = "low",
  manifest      = "/mnt/project/extracts/manifest_full.csv"
)
job$job_id

# In a separate terminal (still on RAP), monitor:
#   dx watch <job-id>
#   dx describe <job-id>

# When the job is done, load the result back into R **inside RAP**:
out_path <- file.path("/mnt/project", paste0(job$output, ".csv"))
dt <- fread(out_path)
# Do not print data rows or share row-level previews with the local agent.

# Decode columns + values
meta <- ukb_metadata_setup(source = "auto", dataset = dataset, cache = TRUE)
dt <- ukb_decode_column_names(dt, metadata = meta, style = "snake")
```

> **Do not** `dx download` `pheno_full_baseline_records.csv` to a laptop.
> If the user needs the data locally, push back — they should run the
> downstream analysis on RAP and only export aggregate results.

---

## Example C — Offline metadata + variable-set planning

Use when: you want to scope a study **before** you have RAP access (e.g.
drafting an analysis plan, or sanity-checking field IDs from a paper).

```r
# Offline mode — needs the two public TSVs (download from UKB Showcase manually)
meta <- ukb_metadata_setup(
  source    = "files",
  data_dict = "Data_Dictionary_Showcase.tsv",
  codings   = "Codings.tsv"
)
print(meta)

# Look up a single field
info <- ukb_field_info("21003", metadata = meta)
print(info$summary)
info$codings

# Plan against a fields_df instead of a live dx call
fields_df <- data.frame(
  field_name = c("participant.eid", "participant.p31_i0",
                 "participant.p21003_i0", "participant.p4080_i0_a0"),
  title      = c("eid", "Sex | Instance 0",
                 "Age at recruitment | Instance 0",
                 "Systolic blood pressure, automated reading | Instance 0 | Array 0")
)
plan_offline <- rap_plan_extract(
  field_id  = c(31, 21003, 4080),
  fields_df = fields_df,
  dataset   = "stub.dataset"
)
plan_offline$matched
plan_offline$n_columns
```

When the user later moves to a real RAP session, the same `rap_plan_extract()`
call (without `fields_df`) reuses the same field IDs and validates them
against the live dataset.

---

## Snippet — choosing sync vs async programmatically

The skill should generate this decision helper when the user is unsure:

```r
choose_extract_mode <- function(plan,
                                threshold_cols = 500,
                                heavy_prefixes = c("^p41", "^p20002",
                                                   "^p23", "^p30", "^p87", "^p126")) {
  big   <- plan$n_columns > threshold_cols
  heavy <- any(unlist(lapply(heavy_prefixes, function(p) any(grepl(p, plan$fields)))))
  if (big || heavy) "async" else "sync"
}

mode <- choose_extract_mode(plan)
cat("Suggested mode:", mode, "\n")
```

---

## Snippet — reproducibility footer (aggregate/shareable)

Append to every script. Captures provenance with **no participant rows**.

```r
log_path <- "/mnt/project/results/extract_session_info.txt"
sink(log_path)
cat("UKBAnalytica version: ", as.character(utils::packageVersion("UKBAnalytica")), "\n")
cat("R version          : ", R.version.string, "\n")
cat("dx version         : ", system2("dx", "--version", stdout = TRUE), "\n")
cat("Dataset            : ", dataset, "\n")
cat("Field IDs          : ", paste(sort(unique(plan$matched$field_id)), collapse = ", "), "\n")
cat("Columns extracted  : ", plan$n_columns, "\n")
cat("Extraction time    : ", format(Sys.time(), tz = "UTC"), " UTC\n")
sink()
```
