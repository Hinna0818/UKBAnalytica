# RAP guardrails for `ukbsci-rap-extract`

These rules are **non-negotiable**. The agent must refuse — or rewrite — any
request that violates them. This file exists so the agent can quote authority
back to the user when pushing back.

---

## 1. Why these rules exist

UK Biobank participant-level data is governed by:

- the **UK Biobank Material Transfer Agreement (MTA)** between UKB and your
  institution / approved application,
- the **DNAnexus RAP terms of service** for the project containing the
  approved dataset, and
- your institution's Research Ethics / Data Governance approval.

All three converge on one principle: **participant-level records must remain
inside the approved RAP project and must not be exposed to a local AI agent**.
The agent should generate scripts and interpret aggregate outputs only. Real
rows may not be moved to a personal laptop, an unapproved bucket, an external
collaborator's machine, or an external LLM/API context.

---

## 2. Forbidden → permitted matrix

| User request | Agent response |
|--------------|----------------|
| "Download the extracted CSV to my laptop." | **Refuse.** Explain rule. Offer to: (a) keep file in `/mnt/project/...`, (b) compute aggregate summary inside RAP and export only that. |
| `dx download <participant-level file>` to local disk | **Refuse.** Same as above. |
| `read.csv("https://.../participant.csv")` (any URL pointing at raw rows) | **Refuse.** Phenotypes must originate from `rap_extract_pheno()` / `rap_submit_extract()`. |
| Writing extracted phenotypes to `~/`, `C:\\Users\\...`, `/Users/...`, `/tmp/` (laptop), etc. | **Refuse.** Force the output path to `/mnt/project/...`. |
| `scp`, `rsync`, `gsutil`, `aws s3 cp` of participant-level files out of RAP | **Refuse.** |
| Pushing participant-level CSV to GitHub / GitLab / Google Drive / Box | **Refuse.** |
| "Send me a sample of the data so I can prototype locally." | **Refuse for real data.** Offer to generate a *simulated* toy dataset with the same column schema and zero real rows. |
| "Cache the extract to my personal cache dir on my laptop." | **Refuse.** Caches must stay on the RAP node. |
| "Open this RDS/CSV so you can debug it." | **Refuse for real RAP data.** Ask for sanitized structure, column names, aggregate counts, package versions, and error text with row values removed. |
| "Here is a screenshot/log with data rows." | **Refuse to inspect row-level content.** Ask for a redacted log or aggregate summary. |

### Permitted (no participant rows involved)

- Downloading the public **UKB Data Dictionary** (`Data_Dictionary_Showcase.tsv`)
  and **Codings** TSV files.
- Downloading the `field_ids.txt` list the user authored.
- Downloading the **manifest CSV** from `rap_plan_extract(..., manifest = ...)`
  (it lists fields, not rows).
- Downloading **aggregate** outputs: regression coefficients, Table 1 summary,
  Kaplan-Meier curves built from group-level counts, ROC/PR data for ML
  models, SHAP global summary plots (per-feature aggregates).
- Pushing **code** (R scripts, Quarto/Rmd source) between local and RAP.

When in doubt, ask the user: *"Is this an aggregate summary or does it contain
participant-level rows or row-level values?"* Aggregate → OK. Rows → keep on
RAP and out of the agent context.

---

## 3. Path discipline

All extraction outputs default to RAP project storage:

```r
# Canonical RAP paths
proj_root <- "/mnt/project"
extract_dir <- file.path(proj_root, "extracts")
result_dir  <- file.path(proj_root, "results")
```

Anti-patterns the agent must not emit:

```r
# WRONG — laptop home directories
output = "~/pheno.csv"
output = "/Users/<name>/Downloads/pheno.csv"
output = "C:/Users/<name>/Documents/pheno.csv"

# WRONG — RAP worker scratch (wiped when worker dies)
output = "/tmp/pheno.csv"
output = "./pheno.csv"

# WRONG — explicitly downloading off the project
system("dx download project-XXXX:/extracts/pheno.csv -o ~/Downloads/")
```

Correct patterns:

```r
output = file.path("/mnt/project/extracts", "pheno_baseline.csv")
# or, equivalently, the DNAnexus URL form:
output = "dx://project-XXXX:/extracts/pheno_baseline.csv"
```

---

## 4. Environment checks the agent must run

Before generating any `rap_*` call, the agent should ensure these checks have
run at least once in the session:

```r
# 4.1 dx is on PATH
stopifnot(nzchar(Sys.which("dx")))

# 4.2 dx is logged in
ok <- system2("dx", "whoami", stdout = TRUE, stderr = TRUE)
if (!length(ok)) stop("dx not authenticated; run `dx login` first.")

# 4.3 We are in the approved project
proj <- system2("dx", c("env", "--bash"), stdout = TRUE)
# Expect a line like: export DX_PROJECT_CONTEXT_ID=project-XXXXXXXXXXX
```

If the user is *not* on RAP and asks "can you just run this locally?":

> No — `rap_*` functions need the DNAnexus `dx` CLI authenticated against an
> approved UK Biobank RAP project. You can still build / inspect *field
> metadata* locally with `ukb_metadata_setup(source = "files", ...)` using the
> public Data Dictionary + Codings TSVs.

---

## 5. Async vs sync — the boundary

`rap_extract_pheno()` blocks the R session. If the user is heading for a
multi-thousand-column extract (hospital codes, all touchscreen instances,
imaging fields, proteomics, metabolomics), force them onto
`rap_submit_extract()` instead. Concrete heuristic:

```r
plan <- rap_plan_extract(field_id = ids, dataset = dataset)
if (plan$n_columns > 500 || any(grepl("^p41|^p20002|^p23|^p30", plan$fields))) {
  message("Use rap_submit_extract() instead of rap_extract_pheno() ",
          "— ", plan$n_columns, " columns.")
}
```

---

## 6. Coded-value sanitation

UK Biobank uses negative codes for special non-answers:

| Code | Meaning |
|------|---------|
| `-1` | Do not know |
| `-2` | Not applicable |
| `-3` | Prefer not to answer |
| `-7` | None of the above (in self-report items) |
| `-10` | Less than one |

Treat these as `NA` before downstream analysis. The standard incantation:

```r
dt <- ukb_decode_values(dt, metadata = meta, missing_to_na = TRUE)
```

For numeric continuous fields where these codes appear without a coding table,
sanitize explicitly:

```r
sentinel <- c(-1, -2, -3, -7, -10)
for (col in numeric_cols) {
  dt[[col]][dt[[col]] %in% sentinel] <- NA
}
```

---

## 7. Reproducibility footer

Every script emitted by the agent should end with a footer that captures the
extraction provenance — but **only metadata**, never participant rows:

```r
# Reproducibility footer
sink("/mnt/project/results/extract_session_info.txt")
cat("UKBAnalytica version: ", as.character(utils::packageVersion("UKBAnalytica")), "\n")
cat("R version          : ", R.version.string, "\n")
cat("dx version         : ", system2("dx", "--version", stdout = TRUE), "\n")
cat("Dataset            : ", dataset, "\n")
cat("Field IDs requested: ", paste(plan$matched$field_id, collapse = ", "), "\n")
cat("Columns extracted  : ", plan$n_columns, "\n")
cat("Extraction time    : ", format(Sys.time(), tz = "UTC"), " UTC\n")
sink()
```

This footer may be shared with the local agent because it contains no
participant rows.
