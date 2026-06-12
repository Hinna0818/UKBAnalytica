# `ukbsci-rap-extract` — function reference

All signatures reflect `UKBAnalytica` 1.0.0
(`R/rap_extract.R`, `R/rap_manifest.R`, `R/ukb_metadata.R`,
`R/field_metadata.R`).

Symbols:
- **RAP-required** — calls `dx`; the R session must be inside an authenticated
  RAP node.
- **Offline** — pure-R, no network or `dx` requirement.

---

## 0. `ukb_check_rap_env()` — Offline + optional dx checks

```r
ukb_check_rap_env(output_dir = NULL,
                  require_rap = FALSE,
                  require_dx = FALSE,
                  check_auth = FALSE,
                  check_write = FALSE,
                  verbose = TRUE)
```

Detects RAP environment signals (`DX_PROJECT_CONTEXT_ID`, workspace/job IDs,
`/mnt/project`) and dx availability. With `require_rap = TRUE` or
`require_dx = TRUE`, stops early if a RAP-only script is being run locally.

**Returns:** list with `is_rap`, `dx_available`, RAP signal details, auth/write
check status, and a checks table.

---

## 1. `rap_find_dataset()` — RAP-required

```r
rap_find_dataset(refresh = FALSE, timeout = 30)
```

Locate the `.dataset` file in the current RAP project root.

| Arg | Type | Default | Meaning |
|-----|------|---------|---------|
| `refresh` | logical | `FALSE` | Bypass the session cache and re-run `dx ls`. |
| `timeout` | int | `30` | Seconds to wait for `dx ls`. |

**Returns:** `character(1)` — the dataset file name (e.g.
`"app99999_20231001154123.dataset"`). Cached in
`.rap_extract_cache$dataset`.

**Errors when:** `dx` not on PATH, multiple `.dataset` files, or none found.

---

## 2. `rap_list_fields()` — RAP-required

```r
rap_list_fields(
  dataset  = NULL,
  pattern  = NULL,
  entity   = "participant",
  refresh  = FALSE,
  timeout  = 120
)
```

List approved RAP dataset fields. Optionally regex-filter at list time.

| Arg | Type | Default | Meaning |
|-----|------|---------|---------|
| `dataset` | char | `NULL` | Dataset file name. If `NULL`, uses `rap_find_dataset()`. |
| `pattern` | char | `NULL` | Regex applied to `field_name` and `title`. |
| `entity` | char | `"participant"` | RAP entity to query. |
| `refresh` | logical | `FALSE` | Bypass per-dataset cache. |
| `timeout` | int | `120` | Seconds for the dx call. |

**Returns:** `data.frame` with columns:

- `field_name` (char) — e.g. `"participant.p31_i0"`
- `title` (char) — e.g. `"Sex | Instance 0"`

Cached per `dataset:entity` pair.

---

## 3. `rap_plan_extract()` — Offline (no dx call)

```r
rap_plan_extract(
  field_id        = NULL,
  field_names     = NULL,
  variables       = NULL,
  dataset         = NULL,
  fields_df       = NULL,
  entity          = "participant",
  include_eid     = TRUE,
  table_exporter  = FALSE,
  manifest        = NULL
)
```

Plan an extraction without running it. Validates and expands the field set.

| Arg | Type | Meaning |
|-----|------|---------|
| `field_id` | int / char vector | UKB numeric field IDs. All instances + arrays included. |
| `field_names` | char vector | Exact RAP column names (`"participant.p31"`, `"p31"`). |
| `variables` | char vector | Predefined variable names from `get_variable_info()`. |
| `dataset` | char | RAP dataset; resolved via `rap_find_dataset()` if `NULL`. |
| `fields_df` | data.frame | Output of `rap_list_fields()` for a faster plan. |
| `entity` | char | RAP entity, default `"participant"`. |
| `include_eid` | logical | Always include the participant ID column. |
| `table_exporter` | logical | Format `fields` for `table-exporter` instead of `dx extract_dataset`. |
| `manifest` | char | Optional CSV path to record matched fields (RAP-resident). |

**Returns:** S3 `rap_extract_plan` (list):

| Element | Type | Description |
|---------|------|-------------|
| `fields` | char vector | Expanded RAP column names to extract |
| `matched` | data.frame | Columns: `request_type`, `request`, `field_id`, `title`, `n_cols`, `field_names` |
| `unmatched` | char vector | Requests that did not resolve |
| `n_columns` | int | Total output columns (sync vs async decision input) |
| `dataset`, `entity`, `field_id`, `variables`, `table_exporter` | metadata |

**Best practice:** check `length(plan$unmatched) == 0` before executing.

---

## 3b. Manifest helpers — Offline

```r
ukb_create_extraction_manifest(field_id = NULL,
                               variable_set = NULL,
                               variables = NULL,
                               dataset = NULL,
                               entity = "participant",
                               output = NULL,
                               include_eid = TRUE,
                               purpose = NULL,
                               notes = NULL)

ukb_write_extraction_manifest(manifest, path,
                              format = c("csv", "rds"))
```

`ukb_create_extraction_manifest()` creates an S3 `ukb_extraction_manifest`
list. The field-level table is stored in `manifest$fields`; it contains no
participant rows. `ukb_write_extraction_manifest()` writes the manifest to CSV
or RDS.

---

## 4. `rap_extract_pheno()` — RAP-required (synchronous)

```r
rap_extract_pheno(
  field_id            = NULL,
  field_names         = NULL,
  variables           = NULL,
  dataset             = NULL,
  output              = NULL,
  read                = TRUE,
  strip_entity_prefix = FALSE,
  dry_run             = FALSE,
  timeout             = 300,
  ...
)
```

Run `dx extract_dataset` synchronously. Use for small/medium extracts only.

| Arg | Type | Default | Meaning |
|-----|------|---------|---------|
| `output` | char | `NULL` | CSV path **inside** RAP storage. If `NULL`, a temp path is used. |
| `read` | logical | `TRUE` | Load CSV into R as `data.table`. If `FALSE`, return the path. |
| `strip_entity_prefix` | logical | `FALSE` | Drop `"participant."` from column names. |
| `dry_run` | logical | `FALSE` | Return the plan only — no dx call. |
| `timeout` | int | `300` | Seconds for the dx call. |

**Returns:**
- `read = TRUE` → `data.table` with `attr(., "rap_extract_plan")` set.
- `read = FALSE` → `character(1)` path.
- `dry_run = TRUE` → S3 `rap_extract_plan`.

**Caveat:** blocks the R session. Do not use for whole-cohort hospital records.

---

## 5. `rap_submit_extract()` — RAP-required (asynchronous)

```r
rap_submit_extract(
  field_id      = NULL,
  field_names   = NULL,
  variables     = NULL,
  dataset       = NULL,
  file          = NULL,
  instance_type = NULL,
  priority      = c("low", "high"),
  dry_run       = FALSE,
  manifest      = NULL,
  ...
)
```

Submit a `table-exporter` job. Preferred for large extracts.

| Arg | Type | Default | Meaning |
|-----|------|---------|---------|
| `file` | char | auto `"ukba_pheno_YYYYMMDD_HHMMSS"` | Output stem on RAP. |
| `instance_type` | char | `NULL` (auto) | DNAnexus instance type. |
| `priority` | char | `"low"` | `"low"` or `"high"`. |
| `dry_run` | logical | `FALSE` | Build the command, do not upload/submit. |
| `manifest` | char | `NULL` | Optional CSV path (RAP-resident). |

**Returns:** S3 `rap_extract_job` (list):

| Element | Type | Description |
|---------|------|-------------|
| `job_id` | char | DNAnexus job ID |
| `dataset` | char | Source dataset |
| `output` | char | File stem on RAP |
| `fields_file_id` | char | Uploaded fields file ID |
| `instance_type` | char | Worker instance type used |
| `priority` | char | `"low"` / `"high"` |
| `n_columns` | int | Output column count |
| `matched` | data.frame | Resolved field requests |
| `unmatched` | char | Unresolved requests |
| `fields` | char vector | Extracted column names |

Monitor outside R via `dx watch <job-id>` or `dx describe <job-id>`.

---

## 6. `ukb_metadata_setup()` — Offline / hybrid

```r
ukb_metadata_setup(
  source     = c("auto", "files", "rap"),
  data_dict  = NULL,
  codings    = NULL,
  fields_df  = NULL,
  dataset    = NULL,
  entity     = "participant",
  cache      = FALSE,
  cache_dir  = NULL,
  refresh    = FALSE,
  quiet      = FALSE
)
```

Build the searchable metadata object that powers `ukb_search_fields()`,
`ukb_field_info()`, and `ukb_decode*()`.

| Arg | Type | Meaning |
|-----|------|---------|
| `source` | char | `"auto"` tries RAP + files; `"files"` is offline; `"rap"` requires RAP. |
| `data_dict` | char | Path to `Data_Dictionary_Showcase.tsv`. |
| `codings` | char | Path to `Codings.tsv`. |
| `fields_df` | data.frame | Output of `rap_list_fields()` (for offline injection). |
| `dataset` | char | RAP dataset file (RAP source only). |
| `cache` | logical | Persist to `tools::R_user_dir("UKBAnalytica","cache")/ukb_metadata.rds`. |

**Returns:** S3 `ukb_metadata` (list) with:

| Element | Description |
|---------|-------------|
| `fields` | Standardized field metadata (field_id, title, description, category, value_type, units, coding, instances, array, participants, …) |
| `rap_fields` | RAP-specific field info (entity, field_id, instance, array, rap_column, field_name, title, title_clean) |
| `codings` | Coding table (coding_id, value, meaning, parent_value, selectable, source) |
| `categories` | Placeholder data.frame |
| `source` | Construction metadata (mode, data_dict, codings, dataset, entity, built_at, cache_dir, rap_error) |

Has a `print()` method.

---

## 7. `ukb_search_fields()` — Offline

```r
ukb_search_fields(
  query       = NULL,
  field_id    = NULL,
  metadata    = NULL,
  max_results = 50,
  search_in   = c("title", "description", "category",
                  "field_name", "rap_field_names", "coding_id"),
  ...
)
```

Keyword / exact-ID search across the metadata object.

**Returns:** `data.frame` (S3 `ukb_search_result`) with the metadata columns
plus `score` (numeric) and `matched_on` (char). `attr(., "query")` records the
input.

---

## 8. `ukb_field_info()` — Offline (with optional `live`)

```r
ukb_field_info(
  x,
  by       = c("auto", "field_id", "title", "rap_column", "variable"),
  metadata = NULL,
  live     = FALSE,
  ...
)
```

Inspect one UKB field's full metadata.

**Returns:** S3 `ukb_field_info` (list):

| Element | Description |
|---------|-------------|
| `field` | 1-row data.frame of field metadata |
| `rap_fields` | All matching RAP columns |
| `codings` | Coding labels (if applicable) |
| `matched_by` | How the field was resolved |
| `query` | Original input |
| `summary` | Pretty-printed character vector |

Has a `print()` method.

---

## 9. `ukb_extract_fields()` — RAP-required for sync/job

```r
ukb_extract_fields(
  x        = NULL,
  field_id = NULL,
  metadata = NULL,
  mode     = c("plan", "sync", "job"),
  top_n    = NULL,
  dataset  = NULL,
  entity   = "participant",
  ...
)
```

Convenience wrapper that consumes a `ukb_search_result` or `ukb_field_info`
and dispatches to `rap_plan_extract` / `rap_extract_pheno` /
`rap_submit_extract`.

| Mode | Delegates to | Returns |
|------|--------------|---------|
| `"plan"` | `rap_plan_extract()` | S3 `rap_extract_plan` |
| `"sync"` | `rap_extract_pheno()` | `data.table` or `character` path |
| `"job"` | `rap_submit_extract()` | S3 `rap_extract_job` |

---

## 10. `get_field_metadata()` — Offline / hybrid

```r
get_field_metadata(
  field_id      = NULL,
  query         = NULL,
  ukb_data_dict = NULL,
  dataset       = NULL,
  fields_df     = NULL,
  entity        = "participant"
)
```

Lower-level than `ukb_metadata_setup()` — returns just the joined metadata
`data.frame`. Errors if *none* of `ukb_data_dict`, `fields_df`, `dataset` are
provided.

**Returns:** `data.frame` with `source_backend` ∈ `{"showcase", "rap", "showcase+rap"}`.

---

## 11. `get_field_info()` — Offline (with optional `live`)

```r
get_field_info(
  field_id,
  ukb_data_dict = NULL,
  dataset       = NULL,
  fields_df     = NULL,
  entity        = "participant",
  live          = FALSE,
  timeout       = 30
)
```

One-shot single-field lookup. If `live = TRUE`, fetches from the public UKB
Showcase web page (requires internet + `xml2`).

---

## 12. `ukb_decode*()` — Offline

```r
ukb_decode(
  data,
  metadata      = NULL,
  decode_names  = TRUE,
  decode_values = TRUE,
  keep_raw      = TRUE,
  suffix        = "_label",
  ...
)

ukb_decode_column_names(
  data,
  metadata       = NULL,
  style          = c("snake", "title", "field_id_title"),
  keep_instance  = TRUE,
  keep_array     = TRUE,
  max_nchar      = 80,
  ...
)

ukb_decode_values(
  data,
  metadata      = NULL,
  keep_raw      = TRUE,
  suffix        = "_label",
  missing_to_na = TRUE,
  ...
)
```

`ukb_decode()` is the combined helper. The two specialized variants give
finer control.

**Returns:** same class as input (`data.frame` or `data.table`).

---

## 13. `ukb_snapshot()` — Offline

Quick column-by-column snapshot (class, NA count, range/distinct, sample
values). Useful right after extraction to detect mis-coded fields or
unexpectedly empty columns.

---

## Internal helpers (not exported, but invoked by all of the above)

- `.rap_require_dx()` — stops with a clear message if `dx` is missing.
- `.rap_check_logical()` — argument validation.
- `.rap_dx_env()` — sets `PYTHONIOENCODING=utf-8` for child dx calls.
- `.safe_as_date()` — coerces mixed UKB date inputs (`YYYYMMDD`, `YYYYMM`,
  R numeric origin, character variants) to `Date`. Returns `NA` and logs
  warnings for non-standard values.

You should not call these directly, but they explain why some apparent
"errors" surface as messages — they are upstream sanitation warnings.

---

## Source backends (`source_backend`)

The metadata `data.frame` reports where each row came from:

| Value | Means |
|-------|-------|
| `"showcase"` | Public UKB Data Dictionary only |
| `"rap"` | RAP-listed field only (no showcase entry) |
| `"showcase+rap"` | Joined / present in both |

In an analysis script, prefer rows with `"showcase+rap"` for completeness.
