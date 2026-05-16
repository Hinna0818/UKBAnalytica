# ukbsci-cohort

Agent skill for building UK Biobank disease cohorts and Cox-ready survival
datasets, using the `UKBAnalytica` R package.

## What this skill does

1. **Picks or builds** a disease definition (`get_predefined_diseases`,
   `create_disease_definition`, `combine_disease_definitions`).
2. **Parses** UKB diagnostic source tables (ICD-10, ICD-9, OPCS-4, death
   registry, cancer registry, first-occurrence, self-report, algorithmic
   outcomes).
3. **Identifies** cases per disease and per source
   (`extract_cases_by_source`, `extract_disease_history`).
4. **Compares** sources for sensitivity analyses
   (`compare_data_sources`, `extract_disease_history_sensitivity`).
5. **Builds** the master Cox-ready survival dataset via
   `build_survival_dataset()`, with prevalent / incident split, follow-up time,
   and participant attrition.
6. **Stratifies** incident cases by follow-up window (`select_incident_by_years`).

## What this skill does NOT do

- Extract raw phenotype tables from RAP → that is `ukbsci-rap-extract`.
- Fit Cox / logistic / linear models → that is `ukbsci-regression` / `ukbsci-survival`.
- Suggest moving participant-level rows out of RAP project storage.

## Compatibility

Provider-agnostic. The YAML front-matter of `SKILL.md` exposes the standard
`name` + `description` keys, parseable by Claude Code skills, OpenAI
Assistants, LangChain tool-style adapters, and any agent runtime that handles
Markdown + YAML.

## File map

```
ukbsci-cohort/
├── SKILL.md
├── README.md
├── references/
│   ├── functions.md          ← full signatures + return values
│   ├── disease-catalog.md    ← curated phenotype list
│   ├── source-matrix.md      ← which source needs which UKB field IDs
│   └── examples.md           ← three end-to-end examples
└── evals/
    └── evals.json
```
