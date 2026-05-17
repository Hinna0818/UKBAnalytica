# Agent Privacy Boundary for UKBAnalytica Skills

These skills are designed for **script generation, workflow planning, package
usage guidance, and interpretation of aggregate outputs**. They are not a
mechanism for a local AI agent to inspect or process real UK Biobank
participant-level RAP data.

## Strict Boundary

The agent may:

- Generate R scripts intended to be run inside the approved UK Biobank RAP
  environment.
- Review package source code, documentation, function signatures, and toy
  examples.
- Use simulated data or small synthetic tables for smoke tests.
- Interpret aggregate results already produced inside RAP, such as baseline
  tables, participant-flow counts, model coefficients, AUCs, enrichment
  results, and rendered figures.
- Help write manuscript methods/results based on aggregate tables and figures.

The agent must not:

- Open, read, summarize, or transform real RAP participant-level files
  (`.csv`, `.tsv`, `.rds`, `.fst`, `.parquet`, `.feather`, `.arrow`, etc.).
- Inspect R objects containing real participant rows, even if direct
  identifiers were removed.
- Request or display `head(data)`, row samples, per-participant predictions,
  row-level SHAP matrices, raw Olink/proteomics matrices, raw RAP fields, or
  exact participant dates.
- Receive copied console output, tracebacks, screenshots, or logs that include
  row-level values from real UK Biobank data.
- Move RAP data to local paths, cloud drives, GitHub, or any external LLM/API
  context.

## Practical Workflow

1. The local agent writes scripts and analysis plans.
2. The user runs those scripts inside RAP.
3. Only aggregate outputs or rendered figures are returned to the agent for
   interpretation or manuscript writing.
4. If debugging is needed, the user should share sanitized error messages and
   object structures (`names()`, `str()` with values removed, counts, and
   package versions), not row-level data.

## Schema-Only Prompt Workflow

When the user describes a real UKB RAP dataframe, the agent may use the
description of column names, intended variable roles, and approximate variable
types to generate code. The agent must not request or inspect real rows,
`head(data)`, row samples, `eid` values, exact dates, raw RAP fields, per-row
predictions, or row-level SHAP matrices.

Recommended workflow:

1. Infer the schema from the user's description.
2. Generate a small synthetic toy dataset with the same column names and
   approximate variable types.
3. Run smoke tests only on the synthetic dataset to verify syntax, function
   signatures, expected output columns, and plotting calls.
4. Clearly separate the synthetic test block from the final RAP script.
5. Instruct the user to run the final script inside RAP on the real dataset.
6. Ask the user to share back only aggregate outputs or rendered
   non-identifying figures.

Synthetic data must be random or toy data. It must not be sampled from UKB
participants and must not be derived from copied real rows.

## Acceptable Aggregate Outputs

Examples of acceptable material to share with the agent:

- Participant-flow counts.
- Missingness percentages by variable.
- Baseline Table 1 summaries.
- Regression coefficient tables.
- Model metric tables.
- ROC/PR/calibration/DCA curve summaries.
- Gene/protein-level association results without participant rows.
- Enrichment and PPI network summaries.
- Rendered manuscript figures, provided they contain no identifiers, exact
  dates, raw RAP fields, row labels, or data previews.

When uncertain, treat the object as participant-level and keep it inside RAP.
