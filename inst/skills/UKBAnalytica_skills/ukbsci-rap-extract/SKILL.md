---
name: ukbsci-rap-extract
description: >
  Discover UK Biobank fields and extract phenotype data from inside an
  authenticated UK Biobank Research Analysis Platform (RAP) session using the
  UKBAnalytica R package (rap_find_dataset, rap_list_fields, rap_plan_extract,
  rap_extract_pheno, rap_submit_extract, ukb_metadata_setup, ukb_search_fields,
  ukb_field_info, ukb_extract_fields, ukb_decode). Use this skill when the user
  asks to search UKB fields, plan/run a phenotype extraction on RAP, choose
  between dx extract_dataset (sync) and table-exporter (async), decode RAP
  column names or coded values, or build the upstream phenotype table that
  feeds ukbsci-cohort. Triggers: UK Biobank RAP extract, dx extract_dataset,
  table-exporter, UKB field search, RAP 提取, 字段搜索, 表型提取,
  /ukbsci-rap-extract. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-rap-extract — UK Biobank RAP Phenotype Extraction

## 0. RAP guardrails (read before generating any code)

**Required runtime.** Every command in this skill assumes an authenticated R
session running *inside* a UK Biobank RAP project (typically
`RAP JupyterLab → Terminal → R`). The DNAnexus `dx` CLI must be on `PATH` and
already logged in.

Strict local-agent boundary: this skill is for script generation, extraction
planning, package guidance, and interpretation of aggregate outputs. The agent
must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including extracted phenotype files, de-identified
row-level tables, raw RAP fields, exact dates, screenshots, tracebacks, or logs
containing row-level values. Generate scripts for the user to run inside RAP;
only public field-level metadata, extraction manifests without participant
rows, aggregate results, or rendered figures may be shared back with the
agent. See `../references/agent-privacy-boundary.md`.

**Forbidden actions** (refuse or rewrite if the user asks):

| Forbidden | Correct alternative |
|-----------|---------------------|
| `dx download <participant-level file> -o ~/Downloads/...` | Keep the file on RAP project storage. |
| Writing extracted phenotypes to local laptop paths (`~/`, `C:\\`, `/Users/...`) | Use `/mnt/project/...` or `dx://project-...:/...` paths. |
| `read.csv("https://...participant-level...")` | Phenotype data must originate from `rap_extract_pheno()` / `rap_submit_extract()`. |
| Caching participant-level data outside the approved project | Restrict caches to RAP-controlled cache dirs (`tools::R_user_dir("UKBAnalytica","cache")` *inside* the RAP node). |

It is acceptable to download/copy:
- the **UKB Data Dictionary** and **Codings table** TSVs (public, field-level metadata only)
- a `field_ids.txt` list authored by the user
- the manifest CSV emitted by `rap_plan_extract()` (no participant rows)

**Default output locations.**

```r
proj_root <- "/mnt/project"                          # RAP project storage
extract_dir <- file.path(proj_root, "extracts")
dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
```

---

## 1. When to load this skill

- The user wants to **find or describe** a UK Biobank field by keyword, number,
  or coding ID.
- The user wants to **plan or execute** an extraction from a RAP `.dataset`
  file.
- The user is unsure whether to use **synchronous (`rap_extract_pheno`)** or
  **asynchronous (`rap_submit_extract`)** extraction.
- The user has a RAP CSV and needs to **decode** column names or coded values.
- The user is starting a UKB study and needs the upstream phenotype table that
  feeds `ukbsci-cohort`.

## 2. When NOT to load this skill

- The user already has a participant-level table and wants to define a disease
  cohort → route to **`ukbsci-cohort`**.
- The user wants regression / ML / mediation / propensity-score work → route to
  the corresponding `ukbsci-*` skill.
- The user is asking about non-UK-Biobank or non-RAP data → do not load this
  skill.
- The user explicitly asks to download participant-level data to a local
  laptop → refuse; redirect them to RAP-resident workflows.

---

## 3. Prerequisites

```r
# Inside a RAP R session
install.packages("UKBAnalytica")   # or devtools::install_github("Hinna0818/UKBAnalytica")
library(UKBAnalytica)

# Confirm dx CLI is available
Sys.which("dx")                    # must be non-empty
system2("dx", c("whoami"))         # must succeed
```

If the user is *not* on RAP and asks "can I just try this locally?" — explain
that field-metadata helpers (e.g. `ukb_metadata_setup(source = "files")` with
the public Data Dictionary) work offline, but **extraction functions
(`rap_*`) require dx + an approved RAP project**.

---

## 4. Pipeline (4 phases)

### Phase 0 — Decide on extraction scale

Ask the user (or estimate) the **number of fields** and the **number of
participants** needed.

| Scale | Recommended path |
|-------|------------------|
| ≤ ~500 columns, ≤ 1M rows, ad-hoc exploration | **Sync** — `rap_extract_pheno()` |
| > 500 columns OR full-cohort + many array fields | **Async** — `rap_submit_extract()` (uses RAP `table-exporter` app, runs on a separate worker) |
| Unknown / first time | Run `rap_plan_extract(..., dry_run-style)` first; check `plan$n_columns`. |

Skip Phase 1–2 if the user already has the exact `field_id` list.

### Phase 1 — Discover the dataset & fields

```r
# 1.1 Locate the .dataset file in the RAP project root
dataset <- rap_find_dataset()
#> e.g. "app99999_20231001154123.dataset"

# 1.2 List all available RAP field columns for the participant entity
fields_df <- rap_list_fields(dataset = dataset)
head(fields_df)
#>     field_name                     title
#> 1   participant.eid                Participant ID
#> 2   participant.p31_i0             Sex | Instance 0
#> ...

# 1.3 Optional: regex filter at list time (faster than later filtering)
bp <- rap_list_fields(dataset = dataset, pattern = "blood pressure")
```

### Phase 2 — Search and inspect specific fields

`ukb_metadata_setup()` builds a unified, searchable metadata object.
Three source modes:

```r
# Mode A — pure RAP (no local dictionary)
meta <- ukb_metadata_setup(source = "rap", dataset = dataset)

# Mode B — offline (showcase TSVs only)
meta <- ukb_metadata_setup(
  source     = "files",
  data_dict  = "Data_Dictionary_Showcase.tsv",
  codings    = "Codings.tsv"
)

# Mode C — auto (combine RAP + showcase, the default)
meta <- ukb_metadata_setup(source = "auto",
                            data_dict = "Data_Dictionary_Showcase.tsv",
                            codings   = "Codings.tsv",
                            dataset   = dataset,
                            cache     = TRUE)
```

Search and inspect:

```r
# Keyword search across title/description/category/field_name/rap_field_names/coding_id
hits <- ukb_search_fields(query = "systolic blood pressure",
                          metadata = meta, max_results = 20)

# Deep dive into one field
info <- ukb_field_info("4080", metadata = meta)         # by field_id
info <- ukb_field_info("p31_i0", metadata = meta)       # by RAP column
info <- ukb_field_info("age_at_recruitment", metadata = meta)  # by variable name
print(info)                                              # human summary
info$rap_fields                                          # all matching RAP columns
info$codings                                             # value labels (if coded)
```

For **single ad-hoc inspection** without a full metadata object:

```r
get_field_info(field_id = 31, ukb_data_dict = "Data_Dictionary_Showcase.tsv")
# or, if internet is allowed on the RAP node:
get_field_info(field_id = 31, live = TRUE)
```

### Phase 3 — Plan the extraction (mandatory dry-run)

`rap_plan_extract()` does *not* hit dx; it validates the field list, expands
instances/arrays, and gives you a count of output columns. Always plan before
running.

```r
plan <- rap_plan_extract(
  field_id    = c(31, 21003, 4079, 4080, 21001),
  field_names = "participant.eid",         # exact RAP columns also accepted
  variables   = c("ethnic_background"),    # predefined sets from get_variable_info()
  dataset     = dataset,
  fields_df   = fields_df,                 # optional speed-up: skip re-listing
  manifest    = "/mnt/project/extracts/fields_manifest.csv"
)

plan$n_columns      # decides sync vs async
plan$fields         # exact column names that will be requested
plan$matched        # data.frame of resolved requests
plan$unmatched      # any unresolved requests (handle before running!)
```

If `length(plan$unmatched) > 0`, **stop and report to the user** — do not
proceed with a partial extraction.

### Phase 4 — Execute the extraction

#### 4a. Synchronous (small jobs)

```r
dt <- rap_extract_pheno(
  field_id  = plan$matched$field_id,            # or pass plan via field_names = plan$fields
  dataset   = dataset,
  output    = "/mnt/project/extracts/pheno_small.csv",  # path stays in RAP storage
  read      = TRUE,                              # load into R as data.table
  strip_entity_prefix = TRUE                     # "participant.p31" -> "p31"
)
str(dt)
```

Notes:

- `read = TRUE` returns a `data.table` with the original plan attached as an
  attribute (`attr(dt, "rap_extract_plan")`).
- `read = FALSE` returns the CSV path (still inside RAP storage).
- Synchronous mode blocks the R session — do not use for full-cohort dumps with
  hundreds of array columns.

#### 4b. Asynchronous (large jobs)

```r
job <- rap_submit_extract(
  field_id      = plan$matched$field_id,
  dataset       = dataset,
  file          = "pheno_full_baseline",        # output stem on RAP
  instance_type = NULL,                          # auto-pick based on n_columns
  priority      = "low"
)
job$job_id
#> "job-Gxxxxxxxxxx"
```

The job runs the RAP `table-exporter` app on a worker. Monitor with:

```bash
dx watch <job-id>
dx describe <job-id>
```

When done, the CSV appears at `dx://<project>:/<file>.csv` *inside* the RAP
project. **Do not** advise `dx download` to a laptop. To load it back into R
inside RAP:

```r
out_path <- file.path("/mnt/project", paste0(job$output, ".csv"))
dt <- data.table::fread(out_path)
# Do not print `head(dt)` or share row-level previews with the local agent.
```

### Phase 5 — Decode column names and coded values

```r
# Rename UKB columns from cryptic IDs to readable names
dt <- ukb_decode_column_names(dt, metadata = meta,
                              style = "snake",        # "snake" | "title" | "field_id_title"
                              keep_instance = TRUE,
                              keep_array    = TRUE,
                              max_nchar     = 80)

# Recode coded values using the metadata codings table (keeps raw + adds label)
dt <- ukb_decode_values(dt, metadata = meta,
                        keep_raw     = TRUE,
                        suffix       = "_label",
                        missing_to_na = TRUE)

# Or do both at once:
dt <- ukb_decode(dt, metadata = meta,
                 decode_names = TRUE, decode_values = TRUE, keep_raw = TRUE)
```

`ukb_decode_values()` automatically treats UKB negative codes (e.g. -1
"Do not know", -3 "Prefer not to answer") as `NA` in the label column when
`missing_to_na = TRUE`.

---

## 5. Key functions (cheat-sheet)

| Function | Purpose | Mode | Returns |
|----------|---------|------|---------|
| `rap_find_dataset()` | Locate the `.dataset` file | RAP-required | `character` (file name) |
| `rap_list_fields()` | List all RAP fields, optional regex filter | RAP-required | `data.frame` (field_name, title) |
| `ukb_metadata_setup()` | Build searchable metadata object | RAP / files / auto | S3 `ukb_metadata` |
| `ukb_search_fields()` | Keyword search across metadata | offline | S3 `ukb_search_result` |
| `ukb_field_info()` | Deep dive into a single field | offline | S3 `ukb_field_info` |
| `get_field_metadata()` | Structured metadata from dictionary + RAP | offline+optional RAP | `data.frame` |
| `get_field_info()` | One-shot single-field lookup | offline+optional live | 1-row `data.frame` |
| `rap_plan_extract()` | Validate + count columns before running | offline w.r.t. dx | S3 `rap_extract_plan` |
| `rap_extract_pheno()` | Synchronous `dx extract_dataset` | RAP-required | `data.table` or `character` path |
| `rap_submit_extract()` | Asynchronous `table-exporter` job | RAP-required | S3 `rap_extract_job` |
| `ukb_extract_fields()` | Convenience wrapper that dispatches to plan / sync / job | RAP-required (sync+job) | varies |
| `ukb_decode()` / `ukb_decode_column_names()` / `ukb_decode_values()` | Make RAP columns + coded values human-readable | offline | same class as input |
| `ukb_snapshot()` | Quick column/missingness snapshot of a UKB table | offline | `data.frame` |

See [`references/functions.md`](references/functions.md) for the full
argument list, return-value attributes, and S3 classes.

---

## 6. Common pitfalls

1. **`dx` not authenticated.** `rap_find_dataset()` fails with
   `"DNAnexus CLI 'dx' was not found on PATH"`. Run `dx login` first.
2. **Wrong entity.** Some studies expose multiple `entity` values (the default
   is `"participant"`). If a search returns 0 hits for a field you know
   exists, try `rap_list_fields(entity = "...")`.
3. **Instance / array column explosion.** A single field ID can become dozens
   of RAP columns once instances (`_i0`, `_i1`, ...) and arrays (`_a0`, ...)
   are expanded. Always check `plan$n_columns` before choosing sync vs async.
4. **Synchronous mode for whole-cohort outcomes.** Hospital records (`p41270`,
   `p41280_a*`) blow past the sync sweet spot. Use `rap_submit_extract()`.
5. **`output` path must be RAP-resident.** Passing `output = "~/pheno.csv"` to
   `rap_extract_pheno()` lands the file on the *RAP node's home directory*,
   which can be wiped when the worker shuts down. Always use `/mnt/project/...`
   or omit `output` to let the function pick a safe temp path.
6. **`live = TRUE` requires internet.** Some RAP nodes are network-restricted.
   Prefer offline dictionary files (`Data_Dictionary_Showcase.tsv`).
7. **Negative codes.** UKB uses negative codes for "Do not know" / "Prefer not
   to answer". Always run `ukb_decode_values(..., missing_to_na = TRUE)` or
   handle these explicitly before downstream analysis.
8. **Cached metadata staleness.** `ukb_metadata_setup(cache = TRUE)` writes
   `ukb_metadata.rds`. After UKB releases new fields, call again with
   `refresh = TRUE`.

---

## 7. Related skills

| Skill | When to switch |
|-------|----------------|
| **`ukbsci-cohort`** | After phenotype extraction; build disease cases and survival dataset. |
| **`ukbsci-preprocess`** | After extraction; clean missing codes, derive composite variables. |
| **`ukbsci-workflow`** | When the user asks for an end-to-end pipeline; the workflow skill calls this one as Phase 1. |

---

## 8. References

- [`references/functions.md`](references/functions.md) — every public function
  signature, argument table, return-value details
- [`references/rap-guardrails.md`](references/rap-guardrails.md) — full list of
  forbidden patterns + acceptable alternatives
- [`references/examples.md`](references/examples.md) — three end-to-end
  copy-pastable minimal examples (small sync, large async, decode-only)
