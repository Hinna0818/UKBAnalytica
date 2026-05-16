# `ukbsci-preprocess` — examples

```r
library(UKBAnalytica); library(data.table)
```

## A. Standard pipeline

```r
dt <- ukb_clean_missing(dt, invalid_codes = c(-1, -3, -7, -9, -11, -13, -15))
dt <- calculate_blood_pressure(dt, type = "sbp")
dt <- calculate_blood_pressure(dt, type = "dbp")
dt <- calculate_air_pollution(dt, pollutants = c("NO2","PM10","PM2.5","NOx"))
dt <- calculate_diet_score(dt, na_handling = "strict")
```

## B. Unified `preprocess_baseline`

```r
dt <- preprocess_baseline(
  dt,
  variables = c("age","sex","bmi","sbp","dbp","ldl_cholesterol","hba1c",
                "crp","smoking_status","alcohol_intake","townsend_index"),
  missing_action = "keep",
  invalid_codes  = c(-1, -3)
)
```

## C. One-shot prep + imputation

```r
ready <- prepare_analysis_dataset(
  dt,
  preprocess_vars = c("age","sex","bmi","sbp","ldl_cholesterol","hba1c"),
  impute = TRUE, impute_method = "pmm", impute_m = 5, seed = 1234
)
# `ready` is a list of 5 imputed data.tables
```

## D. Look up a variable set before extraction

```r
clinical <- get_variable_set("clinical_core")
clinical_ids <- get_variable_set("clinical_core", output = "field_id")
# feed into ukbsci-rap-extract
```
