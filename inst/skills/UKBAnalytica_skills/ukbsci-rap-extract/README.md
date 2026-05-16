# ukbsci-rap-extract

Agent skill for UK Biobank phenotype extraction on the **Research Analysis
Platform (RAP)** using the `UKBAnalytica` R package.

## What this skill does

1. **Discovers** the `.dataset` file in your approved RAP project.
2. **Searches** fields by keyword, field ID, or coding ID.
3. **Plans** an extraction (validates field IDs, expands instances/arrays,
   counts output columns) before any dx call.
4. **Executes** the extraction either synchronously (`dx extract_dataset` via
   `rap_extract_pheno`) or asynchronously (`table-exporter` app via
   `rap_submit_extract`).
5. **Decodes** RAP column names and coded values for downstream analysis.

## What this skill does NOT do

- Define disease phenotypes from the extracted table → that is
  `ukbsci-cohort`.
- Suggest downloading participant-level data to a local laptop. Ever.
- Run regressions, ML, propensity scoring, etc. → those are other skills.

## Compatibility

Works inside any agent runtime that loads `SKILL.md` + `references/*.md` as
plain Markdown plus YAML front-matter (`name`, `description`). Tested mental
models for Claude Code skills and OpenAI Assistants.

## File map

```
ukbsci-rap-extract/
├── SKILL.md                  ← loaded by the router
├── README.md                 ← this file
├── references/
│   ├── functions.md          ← full signatures + return values
│   ├── rap-guardrails.md     ← forbidden-actions matrix
│   └── examples.md           ← three end-to-end examples
└── evals/
    └── evals.json            ← router-recall test cases
```
