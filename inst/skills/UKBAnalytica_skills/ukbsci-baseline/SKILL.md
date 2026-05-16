---
name: ukbsci-baseline
description: >
  Create a stratified baseline characteristics table (Table 1) for a
  UK Biobank Research Analysis Platform (RAP) cohort built with
  UKBAnalytica. Wraps create_baseline_table() over the tableone package,
  producing a printable / CSV-exportable summary with appropriate
  chi-square / Fisher / t-test / Wilcoxon comparisons across a case-control
  or exposure strata column. Use this skill when the user asks for a Table
  1, baseline characteristics, demographic summary, case-vs-control
  comparison, or stratified descriptive statistics. Triggers: Table 1,
  baseline table, baseline characteristics, demographics summary, UKB
  Table 1, 基线表, 基线特征, /ukbsci-baseline. Hard rule: do not export UKB
  RAP raw participant-level data or direct identifiers; aggregate baseline
  summaries are safe to export.
---

# ukbsci-baseline — Table 1 generator for UKBAnalytica cohorts

## 0. RAP guardrails

Shared privacy boundary: do not export UK Biobank RAP individual-level raw
data, direct identifiers (`eid`), exact dates, raw RAP fields, or row-level
source tables that can be linked back to participants. De-identified analytical
figures and aggregate summaries (curves, coefficients, metrics, feature-level
or bin-level source tables) are generally exportable when no identifying or raw
participant-level fields accompany them.

Input: cohort `data.table` from `ukbsci-cohort` (or its preprocessed
version). Output: an aggregate Table 1 (counts, means, medians, p-values)
— **safe to export** off RAP.

---

## 1. When to load

- The user wants a "Table 1" / baseline characteristics summary.
- Stratifying summary statistics by case status, exposure, or treatment arm.
- Comparing means / proportions across two groups with appropriate tests.

## 2. When NOT to load

- Preprocessing missing values / building composite variables →
  `ukbsci-preprocess`.
- Statistical inference on associations beyond simple group comparisons →
  `ukbsci-regression`.
- Balance assessment after propensity-score matching →
  `ukbsci-propensity` (uses standardized mean differences, not p-values).

---

## 3. Prerequisites

```r
library(UKBAnalytica)
library(tableone)        # Suggests: dependency
library(data.table)
```

---

## 4. Pipeline

### Phase 1 — Pick the stratifier

A `case_col` must be **binary** (0/1 or two-level factor). Common choices:

- `<primary>_history` from `build_survival_dataset()` (prevalent vs free).
- `<exposure>_above_median` derived in preprocessing.
- The post-matching treatment indicator from `ukbsci-propensity`.

### Phase 2 — List variables

Decide which columns are factors (categorical) and which are continuous.
`create_baseline_table()` does **not** auto-classify — pass both lists
explicitly.

```r
factor_cols <- c("sex", "ethnic_background", "smoking_status",
                 "Hypertension_history", "Diabetes_history")
cont_cols   <- c("age", "bmi", "sbp", "dbp",
                 "ldl_cholesterol", "hba1c", "crp")
```

### Phase 3 — Generate the table

```r
tab1 <- create_baseline_table(
  data             = cohort,
  case_col         = "AA_history",
  factor_cols      = factor_cols,
  continuous_cols  = cont_cols,
  test             = TRUE          # enable group tests
)

# tab1 is a tableone::TableOne S3 object — `print()` accepts tableone options:
mat <- print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE,
             smd = FALSE, missing = TRUE)
# `mat` is a character matrix — write it directly to CSV
write.csv(mat, "/mnt/project/<area>/04-results/01-baseline_table.csv")
```

### Phase 4 — Statistical tests (when `test = TRUE`)

`tableone` chooses tests automatically:

| Variable type | Default test | Alternative (rare) |
|---------------|--------------|---------------------|
| Continuous, normal-ish | t-test | `oneway.test` (≥ 3 groups) |
| Continuous, skewed | Wilcoxon / Kruskal-Wallis | (manual override via `nonnormal`) |
| Categorical | Chi-square | Fisher's exact (sparse cells) |

To force a Wilcoxon test for specific variables, pass `nonnormal` to
`print()`:

```r
print(tab1, nonnormal = c("crp", "alt"), quote = FALSE, noSpaces = TRUE)
```

---

## 5. Common pitfalls

1. **`case_col` not binary.** A three-level factor will trigger a
   one-way-ANOVA / chi-square across all levels — fine if you want it, but
   confirm with the user.
2. **NA in `case_col`.** Rows with `NA` are dropped by `tableone`. For
   incident-only Table 1, filter on `outcome_status %in% c(0, 1)` *before*
   calling.
3. **Mixed variable lists.** A column listed in both `factor_cols` and
   `continuous_cols` is treated as the one that appears last in the
   `tableone::CreateTableOne` call order — surprises follow. Don't overlap.
4. **`test = TRUE` for very large cohorts.** Tiny clinically irrelevant
   differences become "p < 0.001". Report effect sizes (SMD via `smd = TRUE`
   in `print()`) alongside p-values.
5. **Factor reference / level ordering.** Set with `factor(..., levels = ...)`
   *before* the call — `tableone` orders rows by the factor's `levels()`.
6. **Skewed continuous variables.** Forgetting to flag them via `nonnormal`
   in `print()` yields misleading means and t-tests. Recommend Wilcoxon /
   medians + IQR for biomarkers like CRP, ALT, triglycerides.
7. **`SMD` column.** `print(tab1, smd = TRUE)` adds a standardized mean
   difference column — useful when the cohort is the matched output from
   `ukbsci-propensity` (where p-values are uninformative post-matching).

---

## 6. Key functions

```r
create_baseline_table(data, case_col,
                      factor_cols      = NULL,
                      continuous_cols  = NULL,
                      test             = FALSE)
```

Returns a `tableone::TableOne` S3 object. Print / export options:

| `print()` arg | Effect |
|---------------|--------|
| `quote = FALSE`, `noSpaces = TRUE` | tidy CSV-friendly output |
| `showAllLevels = TRUE` | display every level of categorical vars |
| `nonnormal = c(...)` | use Wilcoxon / median + IQR for listed vars |
| `smd = TRUE` | report standardized mean differences |
| `missing = TRUE` | add a missingness column |
| `exact = c(...)` | force Fisher's exact for listed categoricals |

---

## 7. Related skills

| Skill | When to switch |
|-------|----------------|
| `ukbsci-cohort` | Build cohort + `<disease>_history` first. |
| `ukbsci-preprocess` | Negative-code sanitation + composite variables before Table 1. |
| `ukbsci-propensity` | After PSM, summarize matched cohort with `smd = TRUE`. |
| `ukbsci-plot` | If Table 1 needs companion bar / box plots. |

---

## 8. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md) — three examples
  (vanilla, with SMD post-PSM, with `nonnormal` skewed vars)
