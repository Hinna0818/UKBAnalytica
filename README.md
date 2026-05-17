# UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for UK Biobank RAP Data

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/Hinna0818/UKBAnalytica_v2?style=flat)](https://github.com/Hinna0818/UKBAnalytica_v2/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/Hinna0818/UKBAnalytica_v2)](https://github.com/Hinna0818/UKBAnalytica_v2/commits/main)
[![Visits](https://hits.sh/github.com/Hinna0818/UKBAnalytica_v2.svg)](https://hits.sh/github.com/Hinna0818/UKBAnalytica_v2/)
<!-- badges: end -->

<img src="man/figures/logo.png" align="right" height="139" alt="UKBAnalytica logo" />

**UKBAnalytica** is a high-performance R package for working with UK Biobank
Research Analysis Platform (RAP) data inside approved RAP projects. It focuses on standardized
phenotyping, survival-ready datasets, scalable preprocessing, and downstream analysis.
The package does not ship UK Biobank participant-level source records; examples
use field IDs, simulated toy data, or user-provided tables that remain within
RAP-controlled storage.

**For details, please visit**: [Full documentation for UKBAnalytica](https://hinna0818.github.io/UKBAnalytica_v2/)

![](docs/assets/workflow.png)


## Installation
You can install the development version of `UKBAnalytica` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Hinna0818/UKBAnalytica_v2")
```

## Quick start

```r
library(UKBAnalytica)
library(data.table)

## suppose a participant-level table is available within your approved RAP project
ukb_data <- fread("population.csv")

diseases <- get_predefined_diseases()[
  c("AA", "Hypertension", "Diabetes")
]

analysis_dt <- build_survival_dataset(
  dt = ukb_data,
  disease_definitions = diseases,
  prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death"),
  outcome_sources = c("ICD10", "ICD9", "Death"),
  primary_disease = "AA",
  show_flow = TRUE,
  dt_threads = 8
)
```

## Disease Definition Sources

UKBAnalytica builds disease phenotypes by taking the earliest valid evidence
from multiple UK Biobank sources. Predefined disease definitions are available
through `get_predefined_diseases()`, and custom definitions can be created with
`create_disease_definition()`.

Supported sources include:

- `ICD10`: hospital inpatient diagnosis codes.
- `ICD9`: historical hospital inpatient diagnosis codes.
- `Self-report`: touchscreen/verbal interview illness codes.
- `Death`: primary or contributory death-cause ICD-10 codes.
- `OPCS4`: hospital operative procedure codes.
- `CancerRegistry`: cancer registry diagnosis dates and ICD-10 morphology information.
- `FirstOccurrence`: UKB first-occurrence date fields derived from linked health records.
- `Algorithm`: UKB algorithmically-defined outcomes, such as myocardial infarction, stroke, dementia, asthma, COPD, and Parkinson's disease.

If a selected source is not defined for a disease, it is ignored automatically.
For example, cancer registry evidence is used only when the disease definition
has a `cancer_icd10_pattern`; procedure evidence is used only when
`opcs4_pattern` is defined.

## Minimal RAP Extraction Workflow

Run the following inside a UK Biobank RAP R session. Participant-level data
should remain inside approved RAP projects and RAP-controlled storage.

```r
library(UKBAnalytica)

dataset <- rap_find_dataset()
fields <- rap_list_fields()

meta <- ukb_metadata_setup(fields_df = fields)

ids <- get_variable_set("clinical_core", output = "field_id")

dt <- ukb_extract_fields(
  field_id = ids,
  metadata = meta,
  mode = "sync",
  strip_entity_prefix = FALSE
)

dt <- ukb_decode(dt, metadata = meta)
dt <- ukb_clean_missing(dt, action = "na")

ukb_snapshot(dt, "Clinical core extracted and cleaned")
```


## Basic Survival Dataset

```r
disease_defs <- get_predefined_diseases()[c("Hypertension", "Diabetes", "Stroke")]

analysis_dt <- build_survival_dataset(
  dt = analysis_input,
  disease_definitions = disease_defs,
  primary_disease = "Stroke",
  prevalent_sources = c("ICD10", "ICD9", "Self-report", "Death", "FirstOccurrence"),
  outcome_sources = c("ICD10", "ICD9", "Death", "FirstOccurrence"),
  baseline_col = "p53_i0",
  show_flow = TRUE
)
```

Here `analysis_input` should already contain baseline variables and the
diagnosis/date columns required by the selected sources. The output contains one
`*_history` column per disease, one `*_incident` column per disease, and for the
primary disease a standard survival endpoint: `outcome_status` and
`outcome_surv_time`.

## Skills for UKB data analyses

UKBAnalytica also provides a set of optional analysis skills under
`inst/skills/UKBAnalytica_skills/`. These skills are designed to help users
write and revise RAP-ready R scripts, plan UKB analysis workflows, and improve
the reporting of aggregate results. They cover common UKB tasks such as RAP
field extraction, baseline preprocessing, disease phenotyping, cohort assembly,
regression, survival analysis, machine learning, proteomics analysis, and
publication-oriented visualization.

The skills follow a strict privacy boundary. They must not be used to expose
UK Biobank participant-level data to an agent or any external AI system. Users
should describe only the dataframe schema, such as column names, variable roles,
and the intended analysis; the agent can then generate synthetic toy data for
interface testing and produce final scripts for execution on RAP. Real UKB
rows, `eid` values, exact dates, raw RAP fields, row-level predictions, SHAP
matrices, screenshots, and logs containing row-level values should never be
shared with the agent. Only aggregate outputs, such as flow counts, baseline
tables, regression summaries, model metrics, enrichment results, and
non-identifying rendered figures, should be returned for interpretation.

All real-data workflows generated by these skills are intended to be run inside
the official UK Biobank RAP environment and approved RAP project storage.


## Supplementary Materials
Here we provide some learning materials for UK Biobank in which you may be interested:
- [UK Biobank database browser](https://biobank.ndph.ox.ac.uk/ukb/index.cgi)
- [UK Biobank RAP platform](https://ukbiobank.dnanexus.com/landing)
- [UK Biobank GitHub resources](https://github.com/UK-Biobank)
- [UK Biobank learning guides supported by our team](https://hinna0818.github.io/Bioinfo-SMU/Epidemiology/UK_Biobank/) 

## Star History

<p align="center">
  <a href="https://www.star-history.com/?repos=Hinna0818%2FUKBAnalytica_v2&type=date&legend=top-left">
    <img 
      src="https://api.star-history.com/svg?repos=Hinna0818/UKBAnalytica_v2&type=date&legend=top-left"
      width="600"
      alt="Star History Chart"
    />
  </a>
</p>
