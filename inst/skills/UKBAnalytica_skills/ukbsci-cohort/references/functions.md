# `ukbsci-cohort` — function reference

All signatures reflect `UKBAnalytica` 1.0.0.

---

## 1. `create_disease_definition()`

```r
create_disease_definition(
  name                           = NULL,
  icd10_pattern                  = NULL,
  icd9_pattern                   = NULL,
  sr_codes                       = NULL,
  death_icd10                    = NULL,
  opcs4_pattern                  = NULL,
  first_occurrence_fields        = NULL,
  first_occurrence_source_fields = NULL,
  cancer_icd10_pattern           = NULL,
  cancer_histology               = NULL,
  cancer_behaviour               = NULL,
  algo_date_field                = NULL,
  algo_source_field              = NULL,
  # Deprecated aliases (kept for backwards compatibility):
  icd10                          = NULL,
  icd9                           = NULL,
  self_report                    = NULL
)
```

| Arg | Type | Meaning |
|-----|------|---------|
| `name` | char | Disease label (default `"Custom disease"`). |
| `icd10_pattern` | char (regex) | ICD-10 codes; auto-anchored with `^` if not already. |
| `icd9_pattern` | char (regex) | ICD-9 codes. |
| `sr_codes` | int vector | UKB self-report illness codes (Field 20002). |
| `death_icd10` | char (regex) | Defaults to `icd10_pattern` if `NULL`. |
| `opcs4_pattern` | char (regex) | OPCS-4 procedure codes. |
| `first_occurrence_fields` | int vector | UKB Category 1712 date field IDs (e.g. `131298` for I21). |
| `first_occurrence_source_fields` | int vector | Source field IDs; defaults to `first_occurrence_fields + 1`. |
| `cancer_icd10_pattern` | char (regex) | Cancer registry ICD-10 (Field 40006). |
| `cancer_histology` | int vector | Tumour histology codes (Field 40011). |
| `cancer_behaviour` | int vector | Tumour behaviour codes (Field 40012). Use `3L` for malignant only. |
| `algo_date_field` | int | UKB algorithmic outcome date field (Category 42, e.g. `42016` for COPD). |
| `algo_source_field` | int | Algorithmic outcome source field (e.g. `42017` for COPD). |

**Returns:** plain list with all keys named as above.

**Validation:** if `first_occurrence_source_fields` is supplied, its length
must equal `first_occurrence_fields`.

---

## 1b. `ukb_time_skeleton()`

```r
ukb_time_skeleton(data,
                  id_col = "eid",
                  baseline_col = "p53_i0",
                  birth_year_col = "p34",
                  birth_month_col = "p52",
                  age_col = "p21022",
                  death_date_cols = "^(participant\\.)?p40000_i[0-9]+$",
                  lost_to_followup_col = "p191",
                  admin_censor_date = as.Date("2023-10-31"),
                  keep_source_dates = TRUE)
```

Builds a reusable participant-level time backbone. It standardizes baseline
date, approximate birth date, age at baseline, death date, loss-to-follow-up
date, administrative censoring date, participant-specific `followup_end_date`,
and follow-up time.

**Returns:** `data.table`, one row per participant. It does not define disease
outcomes; pass it to `build_survival_dataset(time_skeleton = ...)`.

---

## 2. `get_predefined_diseases()`

```r
get_predefined_diseases(source = c("curated", "pomegranate", "both"),
                        merge_type = c("intersection", "union"),
                        disease = NULL,
                        supported_only = TRUE)
```

Returns a **named list** of disease definitions. Defaults to 86 curated
UKBAnalytica definitions. `source = "pomegranate"` returns 313
Pomegranate-derived definitions that current parsers can use.
`source = "both", merge_type = "union"` returns a 331-definition expanded
library with matched definitions merged under curated keys.

| Group | Members |
|-------|---------|
| Aortic | `AA`, `TAA`, `AAA` |
| Cardio | `CVD`, `MI`, `STEMI`, `NSTEMI`, `HF`, `Angina`, `Stroke`, `Ischaemic_Stroke`, `Intracerebral_Haemorrhage`, `Subarachnoid_Haemorrhage`, `Stroke_TIA`, `Hypertension`, `Vascular_Disease`, `Arrhythmia`, `Atrial_Fibrillation`, `Ventricular_Arrhythmia`, `AV_Block`, `Intraventricular_Block`, `SVT`, `PAD`, `VTE` |
| Endocrine / Metabolic | `Diabetes`, `T1DM`, `T2DM`, `Hyperlipidemia`, `Thyroid_Disorders` |
| Respiratory | `Asthma`, `COPD`, `Bronchiectasis` |
| Renal | `CKD`, `ESRD`, `Renal_Disease` |
| Cancer | `Lung_Cancer`, `Breast_Cancer`, `Prostate_Cancer`, `Colorectal_Cancer`, `Colon_Cancer`, `Rectal_Cancer`, `Melanoma`, `Non_Melanoma_Skin_Cancer`, `Ovarian_Cancer`, `Uterus_Cancer`, `Oesophageal_Cancer`, `Stomach_Cancer` |
| Neuro / Mental | `Parkinsons`, `Parkinsonism`, `Progressive_Supranuclear_Palsy`, `Multiple_System_Atrophy`, `Dementia`, `Alzheimers_Disease`, `Vascular_Dementia`, `Frontotemporal_Dementia`, `Motor_Neurone_Disease`, `Multiple_Sclerosis`, `Epilepsy`, `Migraine`, `Depression`, `Anxiety`, `Schizophrenia_Bipolar`, `Severe_Mental_Illness`, `Alcohol_Use_Disorder`, `Substance_Use_Disorder` |
| GI | `Dyspepsia`, `Irritable_Bowel_Syndrome`, `Inflammatory_Bowel_Disease`, `Diverticular_Disease`, `Treated_Constipation`, `Chronic_Liver_Disease` |
| Musculo / Bone | `Osteoarthritis`, `Rheumatoid_Arthritis`, `Systemic_Lupus_Erythematosus`, `Fracture` |
| ENT / Eye | `Menieres_Disease`, `Glaucoma`, `Cataract`, `AMD` |
| Other | `Pernicious_Anaemia`, `Psoriasis_Eczema`, `Prostate_Disorders`, `Erectile_Dysfunction`, `MACE` |

Use `names(get_predefined_diseases())` to enumerate at runtime.

---

## 3. Catalog inspection helpers

```r
get_disease_catalog(source = c("all", "curated", "pomegranate"),
                    disease = NULL, code_system = NULL,
                    supported_only = FALSE)
get_pomegranate_diseases(disease = NULL, supported_only = TRUE)
get_pomegranate_source_manifest()
load_pomegranate_portal_coding(path)
```

Use `get_disease_catalog()` to inspect code-level evidence before building a
cohort. It returns `definition_id`, `disease_name`, `source`, `code_system`,
`field_id`, `code`, `match_rule`, `validation_status`, and provenance fields.
The portal audit CSV is optional and must be supplied explicitly through
`load_pomegranate_portal_coding(path = ...)`; it is not required for endpoint
construction.

---

## 4. `combine_disease_definitions(..., name = "Combined")`

```r
combine_disease_definitions(..., name = "Combined")
```

Union multiple disease definitions into one composite endpoint.
Patterns concatenated with `|`, code vectors unioned.

---

## 5. Parsers (one row per diagnosis event)

| Function | Required `dt` columns | Output columns |
|----------|-----------------------|----------------|
| `parse_icd10_diagnoses(dt)` | `eid`, `p41270`, `p41280_a*` | `eid`, `icd10_code`, `diag_date`, `source = "ICD10"` |
| `parse_icd9_diagnoses(dt)` | `eid`, `p41271`, `p41281_a*` | `eid`, `icd9_code`, `diag_date`, `source = "ICD9"` |
| `parse_opcs4_procedures(dt)` | `eid`, `p41272` or `p41272_a*`, `p41282_a*` | `eid`, `opcs4_code`, `diag_date`, `source = "OPCS4"` |
| `parse_cancer_registry(dt)` | `eid`, `p40005_i*`, `p40006_i*`, `p40011_i*`, `p40012_i*` | `eid`, `cancer_icd10_code`, `diag_date`, `cancer_histology`, `cancer_behaviour`, `source = "CancerRegistry"` |
| `parse_death_records(dt)` | `eid`, `p40000_i*`, `p40001_i*`, `p40002_i*_a*` | `eid`, `death_code`, `death_date`, `source = "Death"`, `cause_type` (`"primary"` / `"secondary"`) |
| `get_death_dates(dt)` | `eid`, `p40000_i*` | `eid`, `death_date` (earliest per participant) |
| `parse_self_reported_illnesses(dt, baseline_col = "p53_i0")` | `eid`, `p20002_i*_a*`, `p20008_i*_a*` | `eid`, `sr_code`, `diag_date` (month-precision; UKB sentinels removed), `source = "Self-report"`, `instance`, `array_idx` |

All parsers return empty `data.table` (with the right column schema) when the
required columns are missing — they emit a `message()` describing the
omission. The agent must check `nrow()` before relying on a parser result.

---

## 6. `extract_cases_by_source()` — flexible primitive

```r
extract_cases_by_source(
  dt,
  disease_definitions,
  sources       = c("ICD10", "ICD9", "Self-report", "Death"),
  censor_date   = as.Date("2023-10-31"),
  baseline_col  = "p53_i0"
)
```

| Arg | Type | Meaning |
|-----|------|---------|
| `dt` | data.table | UKB phenotype table (auto-coerced from data.frame). |
| `disease_definitions` | named list | Each element a `create_disease_definition()` result. |
| `sources` | char vector | Any subset of `c("ICD10", "ICD9", "Self-report", "Death", "OPCS4", "CancerRegistry", "FirstOccurrence", "Algorithm")`. |
| `censor_date` | Date | Administrative censoring date. |
| `baseline_col` | char | Baseline date column name. |

**Returns:** `data.table` with columns:

| Column | Meaning |
|--------|---------|
| `eid` | Participant ID |
| `disease` | Disease key from `disease_definitions` |
| `status` | `1` = event, `0` = censored |
| `surv_time` | Years from baseline to event-or-censor |
| `prevalent_case` | logical |
| `earliest_date` | Date — earliest valid diagnosis from `sources` |
| `diagnosis_source` | char — winning source. For Algorithm, refined to `"Algorithm_<source-code>"`. |

**Behavior:**

- Silently skips disease definitions missing the field needed for a requested
  source (e.g. asks for `"Self-report"` but `sr_codes` is `NULL`).
- For `"FirstOccurrence"`: excludes UKB special date codes (`1900-01-01`,
  `1901-01-01`, `2037-07-07`, etc.).
- For `"Algorithm"`: requires `algo_date_field`; further filters out the
  `1900-01-01` sentinel ("unknown date").

---

## 6b. `extract_disease_diagnosis()` — selected disease diagnosis table

```r
extract_disease_diagnosis(dt,
                          disease,
                          disease_definitions = NULL,
                          sources = c("ICD10", "ICD9", "Self-report", "Death"),
                          censor_date = as.Date("2023-10-31"),
                          baseline_col = "p53_i0",
                          include_all = TRUE)
```

Extracts diagnosis status for one or more selected diseases. If
`disease_definitions = NULL`, it uses `get_predefined_diseases()`.

**Returns:** `data.table` with `eid`, `disease`, `diagnosed`,
`prevalent_case`, `incident_case`, `earliest_date`, `diagnosis_source`,
`status`, and `surv_time`.

---

## 7. `extract_disease_history()`

```r
extract_disease_history(
  dt,
  diseases,
  disease_definitions = NULL,
  sources              = "ICD10",
  baseline_col         = "p53_i0"
)
```

Adds one boolean column per disease — `<Disease>_history` — to the wide
`dt`. Sources can be any subset of the eight valid values.

**Validation:** stops if any name in `diseases` is missing from
`disease_definitions`.

---

## 8. `extract_disease_history_sensitivity()`

```r
extract_disease_history_sensitivity(
  dt,
  diseases,
  disease_definitions = NULL,
  baseline_col         = "p53_i0"
)
```

Convenience: adds three columns per disease

| Column | Sources |
|--------|---------|
| `<Disease>_history_ICD10` | `c("ICD10")` |
| `<Disease>_history_hospital` | `c("ICD10", "ICD9")` |
| `<Disease>_history_all` | `c("ICD10", "ICD9", "Self-report")` |

Used to power downstream `ukbsci-subgroup-sensitivity` analyses.

---

## 9. `compare_data_sources()`

```r
compare_data_sources(
  dt,
  disease_definitions,
  baseline_col = "p53_i0"
)
```

Returns a per-disease summary `data.table` with counts and overlaps across
ICD-10, ICD-9, hospital (ICD-10 ∪ ICD-9), self-report, and all-source
combinations. Useful before deciding `prevalent_sources` / `outcome_sources`
for `build_survival_dataset()`.

---

## 10. `build_survival_dataset()` — main entry point

```r
build_survival_dataset(
  dt,
  disease_definitions,
  prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources   = c("ICD10", "ICD9", "Death"),
  censor_date       = as.Date("2023-10-31"),
  baseline_col      = "p53_i0",
  time_skeleton     = NULL,
  primary_disease   = NULL,
  output            = c("wide", "long"),
  include_all       = TRUE,
  show_flow         = TRUE,
  dt_threads        = NULL
)
```

| Arg | Type | Meaning |
|-----|------|---------|
| `prevalent_sources` | char vector | Sources used to mark prevalent cases. Broader by default. |
| `outcome_sources` | char vector | Sources used for incident outcome. Excludes Self-report by default (imprecise dates). |
| `time_skeleton` | data.frame / NULL | Optional output from `ukb_time_skeleton()`. |
| `primary_disease` | char vector | One or more disease keys. One key drives `outcome_*`; multiple keys create disease-specific status/time columns. |
| `output` | char | `"wide"` (default) or `"long"`. |
| `include_all` | logical | If `output = "long"`, include non-cases with `status = 0`. |
| `show_flow` | logical | Print/attach participant attrition table. |
| `dt_threads` | int | Temporarily set `data.table::setDTthreads()`; restores on exit. |

**Returns (wide):** the original `data.table` plus, **per disease**:

- `<Disease>_history` (int 0/1)
- `<Disease>_incident` (int 0/1)

For one `primary_disease`:

- `outcome_status` (int / NA — NA = prevalent / not at risk)
- `outcome_surv_time` (numeric years / NA)

For multiple `primary_disease` keys:

- `<Disease>_status` (int / NA)
- `<Disease>_surv_time` (numeric years / NA)

If `show_flow = TRUE`, the attrition table is stored at
`attr(result, "participant_flow")`.

**Returns (long):** one row per participant–disease with columns from
`extract_cases_by_source()` plus `eid`.

---

## 11. `select_incident_by_years()`

```r
select_incident_by_years(
  df,
  n_years        = 5,
  time_col       = "outcome_surv_time",
  status_col     = "outcome_status",
  baseline_col   = "p53_i0",
  event_date_col = "earliest_date",
  group_col      = "incident_timing",
  output         = c("combined", "split"),
  copy           = TRUE,
  verbose        = TRUE
)
```

Stratifies incident cases by years since enrollment.

| `output` | Returns |
|----------|---------|
| `"combined"` (default) | One `data.frame`/`data.table` with `group_col` taking values `"within_<n>_years"` / `"after_<n>_years"`. |
| `"split"` | Named list with elements `within_<n>_years` and `after_<n>_years`. |

If `time_col` is absent in `df`, it is computed from `baseline_col` +
`event_date_col`. If `status_col` is absent, it is inferred from
`event_date_col` being non-NA.

**Validation:** errors on `n_years <= 0`, on `status_col` containing values
other than `{0, 1}`, on invalid date formats.

---

## Notes for the agent

- Always inspect `attr(cohort, "participant_flow")` after
  `build_survival_dataset()` and present it to the user — the attrition table
  may be shared with the local agent because it contains only aggregate counts.
- For sensitivity comparisons, prefer `extract_disease_history_sensitivity()`
  + a downstream `runmulti_*()` call to a custom loop.
- For composite endpoints, use `combine_disease_definitions()` rather than
  manually concatenating regex patterns — it handles the `name` and source
  field bookkeeping for you.
