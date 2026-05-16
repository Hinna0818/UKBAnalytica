# Predefined variable catalogue

Snapshot of `get_variable_info()` — verify in-session for the live state.

| Category | Examples |
|----------|----------|
| `demographics` | `age`, `sex`, `ethnic_background`, `assessment_date`, `assessment_centre` |
| `anthropometrics` | `bmi`, `waist`, `hip`, `whr`, `weight`, `height`, `body_fat_percentage` |
| `lifestyle` | `smoking_status`, `alcohol_intake`, `physical_activity`, `sleep_duration` |
| `socioeconomic` | `townsend_index`, `qualifications`, `household_income`, `employment` |
| `blood_pressure` | `sbp_auto`, `dbp_auto`, `sbp_manual`, `dbp_manual`, `pulse_rate` |
| `medications` | `bp_medication`, `cholesterol_medication`, `diabetes_medication`, `hormone_therapy` |
| `biomarkers` | `hba1c`, `ldl_cholesterol`, `hdl_cholesterol`, `triglycerides`, `crp`, `alt`, `ast`, `egfr` |
| `pollution` | `no2_2010`, `pm10_2010`, `pm25_2010`, `nox_2010` |
| `diet` | `fruit_intake`, `vegetable_intake`, `fish_intake`, `meat_intake`, `cereal_intake`, `milk_type` |

Use `get_variable_set("clinical_core")` and `get_variable_set("air_pollution")`
for curated bundles used across most epidemiological studies.

Each set returns columns `set, category, variable, field_id, ukb_col, label,
role, notes`. `role` is one of `"exposure"`, `"covariate"`, `"outcome"`,
`"identifier"`, or `"derived"`.

## Negative-code semantics (UKB convention)

| Code | Meaning |
|------|---------|
| `-1` | Do not know |
| `-2` | Not applicable |
| `-3` | Prefer not to answer |
| `-7` | None of the above |
| `-10` | Less than one |
| `-11` | Less than once a year |
| `-13` | Never went to school |
| `-15` | Less than once a week |

`ukb_clean_missing()`'s default `invalid_codes` covers the "non-answer"
subset; the field-specific ones (`-10`, `-11`, `-13`, `-15`) are included
because they distort means but are domain-specific — confirm with the user
before stripping.
