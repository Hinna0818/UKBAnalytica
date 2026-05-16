# `ukbsci-preprocess` — function reference

All signatures verbatim from `R/variable_preprocess.R` + `R/variable_sets.R`.

---

## `ukb_clean_missing()`

```r
ukb_clean_missing(df, cols = NULL,
                  invalid_codes = c(-1, -3, -7, -9, -11, -13, -15))
```

| Arg | Meaning |
|-----|---------|
| `cols` | `NULL` → all numeric columns; or char vector |
| `invalid_codes` | UKB sentinel set (defaults match "Prefer not to answer", "Do not know", etc.) |

**Returns:** data.table with NAs substituted in.

---

## `calculate_blood_pressure()`

```r
calculate_blood_pressure(df, type = c("sbp","dbp"),
                         prefer = c("auto","manual"))
```

Computes mean across two readings; falls back to manual when automated
readings missing (or vice versa with `prefer = "manual"`).

**Returns:** data.table with a `sbp` or `dbp` column added.

---

## `calculate_air_pollution()`

```r
calculate_air_pollution(df, pollutants = c("NO2","PM10","PM2.5","NOx"))
```

Averages each pollutant across time-stamped UKB exposure fields.

---

## `calculate_diet_score()`

```r
calculate_diet_score(df,
                     components  = c("fruit","vegetable","fish","meat",
                                     "cereal","milk"),
                     na_handling = c("strict","partial"))
```

Healthy-diet score on 0–7 scale (or normalised to 7 when
`na_handling = "partial"` averages over available components).

---

## `preprocess_baseline()`

```r
preprocess_baseline(df, variables,
                    custom_mapping = NULL,
                    missing_action = c("keep","drop"),
                    invalid_codes  = c(-1, -3))
```

Unified pipeline using the predefined variable map. Errors when a name in
`variables` is not in the catalogue (unless `custom_mapping` provides it).

---

## `prepare_analysis_dataset()`

```r
prepare_analysis_dataset(df,
                         preprocess_vars = NULL,
                         impute          = FALSE,
                         impute_method   = "pmm",
                         impute_m        = 5,
                         seed            = NULL)
```

Convenience: preprocessing + optional one-shot MI via `run_imputation()`.
Use `ukbsci-imputation` for full MI pooling.

---

## Catalogue helpers

```r
get_variable_info(category = "all")
get_variable_set(set,
                 output = c("data.frame","field_id","ukb_col"))
get_variable_sets(set = NULL, category = NULL)
```

`get_variable_info()` valid categories: `"all", "demographics",
"anthropometrics", "lifestyle", "socioeconomic", "blood_pressure",
"medications", "biomarkers", "pollution", "diet"`.
