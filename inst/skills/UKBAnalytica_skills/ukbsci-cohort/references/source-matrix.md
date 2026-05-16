# Source × UKB-field matrix

Cross-reference for which UKB fields each diagnostic source needs. Use this to:

1. Check the extract from `ukbsci-rap-extract` actually contains the columns
   the user's chosen sources require.
2. Decide which sources to add when expanding a cohort definition.
3. Spot quick fixes when a source returns 0 cases unexpectedly.

---

## Required RAP columns per source

| Source key | Required field IDs / RAP columns | Notes |
|------------|---------------------------------|-------|
| `"ICD10"` | `p41270`, `p41280_a*` | Diagnosis array + date array, hospital inpatient |
| `"ICD9"` | `p41271`, `p41281_a*` | Historical hospital records (pre-2001 in most UKB centres) |
| `"OPCS4"` | `p41272` or `p41272_a*`, `p41282_a*` | Operative procedures; expanded or list-string layout supported |
| `"Death"` | `p40000_i*`, `p40001_i*`, `p40002_i*_a*` | Death date + primary + secondary causes |
| `"CancerRegistry"` | `p40005_i*`, `p40006_i*`, `p40011_i*`, `p40012_i*` | Date + ICD-10 + histology + behaviour |
| `"FirstOccurrence"` | `p<131xxx>` + `p<131xxx + 1>` (configurable) | UKB Category 1712 |
| `"Algorithm"` | `p<algo_date_field>` + (`p<algo_source_field>`) | UKB Category 42 outcomes |
| `"Self-report"` | `p20002_i*_a*`, `p20008_i*_a*` | Touchscreen illness code + fractional year |

**Always required:** `eid`, `p53_i0` (baseline assessment date).

---

## Sentinel value handling

The parsers automatically drop UKB sentinel dates that mean "unknown" or
"before records started":

| Sentinel | Source | Handling |
|----------|--------|----------|
| `1900-01-01` | Algorithm (Cat 42) | Dropped (unknown date) |
| `1900-01-01`, `1901-01-01`, `1902-02-02`, `1903-03-03`, `2037-07-07` | First Occurrence | Dropped (UKB "before / after available linkage") |
| Self-report year `-1`, `-3` | Self-report | Dropped (do-not-know / prefer not to answer) |
| Self-report year `< 1900` | Self-report | Dropped (likely encoding error) |

The user should **not** override this — these dates are not meaningful for
event times.

---

## Field-ID quick lookup (Category 1712 First Occurrence)

A subset most-commonly used in `first_occurrence_fields`:

| Disease region | Date field IDs |
|----------------|----------------|
| Ischaemic heart disease (I20–I25) | 131296, 131298, 131300, 131302, 131304, 131306 |
| Cerebrovascular (I60–I69) | 131360, 131362, 131364, 131366, 131368 |
| Atherosclerosis & PAD (I70–I79) | 131380, 131382, 131384, 131386 |
| Diabetes (E10–E14) | 130706, 130708, 130710, 130712, 130714 |
| COPD (J40–J47) | 131490, 131492, 131494, 131496 |
| Dementia (F00–F03, G30) | 130836, 130838, 130840, 131036 |
| Parkinson's (G20–G22) | 131022, 131024, 131026 |
| CKD (N17–N19) | 132036, 132038, 132040 |

The `_source_field` is always `date_field + 1` unless explicitly overridden.

---

## Recommended source combinations

| Use case | `prevalent_sources` | `outcome_sources` | Why |
|----------|---------------------|-------------------|-----|
| Default Cox study | `c("ICD10","ICD9","Self-report","Death")` | `c("ICD10","ICD9","Death")` | Broad exclusion of prevalent; tight (precise-date) outcomes |
| Strict hospital-only sensitivity | `c("ICD10","ICD9")` | `c("ICD10","ICD9")` | Recovers smaller-but-cleaner cohort |
| Algorithm-validated outcome | `c("ICD10","ICD9","Self-report","Death")` | `c("Algorithm")` | When UKB algorithm exists for the disease (MI, Stroke, COPD, Dementia, Parkinson, Asthma) |
| Cancer-only outcome | `c("ICD10","ICD9","Self-report","Death","CancerRegistry")` | `c("CancerRegistry","Death")` | Registry is gold standard for cancer dates |
| First-occurrence-centric | `c("ICD10","ICD9","Self-report","Death","FirstOccurrence")` | `c("FirstOccurrence","Death")` | Catches GP records that hospital ICD doesn't |

When in doubt, run `compare_data_sources(dt, defs)` first and let the user
see the overlap before locking in the combination.

---

## Pre-flight column check (snippet for the agent)

The agent should run this before generating any `extract_cases_by_source()`
call. It is a *check*, not an extraction.

```r
need_for_source <- list(
  ICD10           = c("p41270", "p41280_a0"),
  ICD9            = c("p41271", "p41281_a0"),
  OPCS4           = c("p41272", "p41282_a0"),
  Death           = c("p40000_i0", "p40001_i0"),
  CancerRegistry  = c("p40005_i0", "p40006_i0", "p40011_i0", "p40012_i0"),
  `Self-report`   = c("p20002_i0_a0", "p20008_i0_a0")
)

audit_sources <- function(dt, sources) {
  cols <- colnames(dt)
  problems <- list()
  for (s in sources) {
    needs <- need_for_source[[s]]
    if (is.null(needs)) next
    pat <- sub("_a0$", "_a", needs)
    missing <- !sapply(pat, function(p) any(startsWith(cols, p)))
    if (any(missing)) problems[[s]] <- needs[missing]
  }
  problems
}

audit_sources(dt, c("ICD10", "ICD9", "Self-report"))
```

If `audit_sources()` returns a non-empty list, **stop** and route the user
back to `ukbsci-rap-extract` to pull the missing columns. Don't silently run
an extraction that will return zero cases.
