# UKBAnalytica News

## UKBAnalytica 0.2.0 (2026-01-25)

### Major changes
- Refactored the primary analysis interface to return a cohort-retaining **wide-format** dataset by default.
- Added a `primary_disease` argument to compute `outcome_status` and `outcome_surv_time` for a single primary endpoint.
- Added `prevalent_sources` and `outcome_sources` argument into `build_survival_dataset` function to manage self-report bias.

### New features
- Multi-source phenotyping support with configurable `sources` (ICD-10, ICD-9, self-report, death).
- Cohort-level follow-up time computation with administrative censoring and death censoring.
- Expanded predefined disease definitions to cover common conditions for rapid prototyping.

### Data acquisition (RAP)
- Added Python utilities under `inst/python/` to download:
  - Demographic fields (user-specified UKB field IDs; optional ID file input).
  - Metabolomics (all fields; plus non-ratio subset driven by `inst/extdata/metabolites_non_ratio.txt`).
  - Proteomics (batch download with optional merge).

### Documentation
- Updated README to prioritize data acquisition (RAP) before survival endpoint construction.
- Added a package overview figure in `man/figures/`.

## UKBAnalytica 0.1.0

- Initial release: parsing UKB RAP exports and generating survival-analysis-ready datasets.