---
name: ukb-research
description: Design and evaluate UK Biobank research topics from user-provided keywords, exposures, outcomes, omics layers, or clinical domains. Use this skill to perform literature-informed topic scouting, summarize existing work and its limitations, identify research gaps, propose gap-resolving epidemiological or omics study designs, verify PMID/DOI references, identify required software and UKBAnalytica functions, and return a structured academic Markdown research brief. Triggers: UK Biobank research idea, UKB topic selection, feasibility analysis, PubMed scoping review, literature gap analysis, study design, research proposal, UKB选题, 研究设计, 文献调研, gap分析, 可行性分析, /UKB-research.
---

# UKB Research — topic scouting and study-design brief

## 0. Scope and privacy boundary

Use this skill to help users turn a broad UK Biobank idea into a feasible,
literature-informed study plan. The output should be a polished Markdown
research brief written in academic English.

Privacy boundary:

- Do not request, read, inspect, summarize, or process real UK Biobank
  participant-level data, screenshots, row-level logs, exact dates, `eid`
  values, or per-row model outputs.
- It is acceptable to use public literature, UKB Showcase metadata, public
  codelists, package documentation, and user-provided variable names or field
  IDs.
- Generate scripts or workflow plans for execution inside the official UKB RAP
  environment. Interpret only aggregate tables, figures, and summary metrics.

## 1. When to use

Use this skill when the user asks to:

- identify a promising UK Biobank research topic;
- evaluate feasibility for a keyword, exposure, phenotype, omics layer, or
  disease area;
- summarize PubMed or related literature for study design;
- compare existing studies, identify their limitations, and define the
  research gap the new study should address;
- propose a complete UKB analysis workflow;
- map the workflow to R packages, Python libraries, external tools, and
  UKBAnalytica functions;
- generate a Markdown research proposal or scoping report.

Do not use this skill for direct code execution on real participant-level
data. Route implementation details to the relevant UKBAnalytica workflow skill
after the study design is agreed.

## 2. Required inputs

Ask for the minimum missing information only. A good prompt can contain:

- disease or phenotype keyword, e.g. "COPD", "atrial fibrillation";
- exposure or feature domain, e.g. "air pollution", "diet", "proteomics";
- preferred study type, e.g. prospective cohort, cross-sectional, prediction,
  mediation, MR, omics profiling;
- target output, e.g. methods section, case-study plan, grant idea, figure plan.

If the user provides only broad keywords, proceed with a scoping design and
state assumptions explicitly.

## 3. Literature search workflow

When the user asks for evidence or feasibility, perform a focused web/literature
search unless they explicitly prohibit browsing.

Prioritize sources in this order:

1. PubMed-indexed papers and PMC full text when available.
2. Journal pages or supplementary material for phenotype/codelist details.
3. UK Biobank Showcase, RAP documentation, and public UKB codelist resources.
4. GitHub repositories only when they are linked by a publication or are clearly
   maintained as public phenotype-code resources.

For each search:

- use 2-4 targeted queries combining the keyword, "UK Biobank", exposure,
  outcome, and study type;
- record the publication year, cohort, exposure, outcome, main method, main
  finding, and whether codes/variables are reusable;
- distinguish established findings from unresolved gaps or speculative
  opportunities;
- cite sources with PMID/DOI/URL when available.

For reusable search patterns, read
[`references/search-strategies.md`](references/search-strategies.md). Use it
when the user gives only broad keywords or when the topic spans omics,
prediction, mediation, GWAS, or MR.

## 4. Existing-work and gap analysis

Every report must move beyond a generic literature summary. For the user's
research idea, explicitly answer:

1. **What has already been done?** Identify closely related UKB and non-UKB
   studies, their design, exposure/outcome definitions, data modality, and
   main conclusions.
2. **What are the limitations?** Note gaps in temporality, phenotype
   ascertainment, exposure measurement, confounding control, omics coverage,
   external validation, interpretability, causal inference, and reproducibility.
3. **What is the actionable gap?** Define a gap that can plausibly be addressed
   using UK Biobank and available software.
4. **How can the new study address it?** Propose concrete design choices,
   methods, sensitivity analyses, and expected outputs.

Citation rule: any claim about prior work, limitation of prior work, or gap
must be supported by a verified PMID, DOI, or source URL. If a reference cannot
be verified, mark it as "unverified" and do not use it as primary support.

## 5. Reference verification

Before finalizing a Markdown report, run the bundled verifier when a file is
available:

```bash
python3 inst/skills/ukb-research/scripts/verify_references.py report.md \
  --output reference_audit.json
```

If network access is unavailable, run:

```bash
python3 inst/skills/ukb-research/scripts/verify_references.py report.md \
  --offline --output reference_audit.json
```

Use the audit to catch malformed or invented PMIDs/DOIs. Do not fabricate
references. Prefer PubMed IDs and DOIs over vague author-year citations. If
the user needs reference-manager files, export BibTeX or RIS:

```bash
python3 inst/skills/ukb-research/scripts/verify_references.py report.md \
  --output reference_audit.json \
  --bibtex references.bib \
  --ris references.ris
```

## 6. Feasibility assessment

Score feasibility qualitatively as **High**, **Moderate**, or **Low** across:

- phenotype ascertainment in UKB;
- exposure availability and measurement quality;
- sample size and expected event count;
- temporality and censoring suitability;
- covariate availability;
- omics or imaging availability, if relevant;
- reproducibility and public-code support;
- novelty relative to prior UKB literature.

Flag common blockers:

- rare outcomes with too few incident cases;
- exposure measured after outcome onset;
- missing or ambiguous GP-only phenotypes;
- unavailable biomarker/omics platform for the intended population;
- high risk of reverse causation or collider bias;
- analysis requiring external GWAS/MR pipelines not yet implemented in
  UKBAnalytica.

## 7. Study design rules

Propose one primary design and, where useful, one or two alternatives.

For each design, specify:

- population and exclusions;
- exposure definition;
- outcome definition and data sources;
- time zero, follow-up, event, and censoring rules;
- covariates and confounder rationale;
- primary model;
- secondary, sensitivity, and subgroup analyses;
- how the design addresses the literature gap;
- expected direction of results or plausible null scenarios;
- key figures and tables.

Prefer designs that can be implemented with UKBAnalytica on RAP:

- `rap_plan_extract()` / `ukb_create_extraction_manifest()` for field planning;
- `get_variable_info()` / `get_variable_sets()` for variable discovery;
- `preprocess_baseline()` for baseline covariate harmonization;
- `get_predefined_diseases()` / `get_disease_catalog()` for endpoint codes;
- `build_survival_dataset()` for prevalent/incident separation and follow-up;
- `create_baseline_table()` for baseline summaries;
- `runmulti_cox()`, `runmulti_logit()`, `runmulti_lm()`, `runmulti_glm()`,
  `runmulti_negbin()`, `runmulti_gam()`, `run_rcs()` for modeling;
- `ukb_ml_flow()` / `ukb_ml_compare()` for machine-learning workflows;
- proteomics, metabolomics, enrichment, PPI, forest, volcano, ROC, SHAP, KM,
  calibration, and DCA plotting helpers when relevant.

## 8. Software mapping

Always include an implementation table with:

- **UKBAnalytica** functions for RAP extraction planning, preprocessing,
  phenotyping, survival data construction, regression, ML, omics, and plotting;
- **R packages** such as `data.table`, `survival`, `tableone`, `pROC`,
  `ggplot2`, `rms`, `clusterProfiler`, `STRINGdb`, `xgboost`, or `shapviz`
  when relevant;
- **Python libraries** such as `pandas`, `scikit-learn`, `lifelines`,
  `statsmodels`, or `shap` only if the analysis is better suited to Python;
- **external tools** such as PLINK, SAIGE, BOLT-LMM, METAL, LDSC, TwoSampleMR
  resources, or GWAS Catalog when the proposal requires large-scale genotype,
  GWAS, or MR analysis.

State which parts can be handled directly by UKBAnalytica and which require
external software.

## 9. Required output

Return a complete Markdown document. If the user asks for a file, write a
`.md` file to the requested path or a sensible project path.

Use the report structure in
[`references/report-template.md`](references/report-template.md). Keep the
main text concise, with evidence and workflow tables for detail.

Quality bar:

- academic English;
- explicit assumptions;
- citations or source links for literature claims;
- explicit "existing work → limitations → gap → proposed solution" logic;
- reference audit completed when possible;
- no unsupported claims of novelty;
- no participant-level data exposure;
- actionable next steps.

Before modifying the verifier script, run its offline smoke test:

```bash
python3 inst/skills/ukb-research/scripts/smoke_test_reference_verifier.py
```

## 10. Handoff to other skills

After the research brief is approved:

- use `ukbsci-rap-extract` for extraction manifests;
- use `ukbsci-preprocess` and `ukbsci-cohort` for baseline variables and
  phenotype definition;
- use `ukbsci-regression`, `ukbsci-survival`, `ukbsci-ml`, `ukbsci-proteomics`,
  or `ukbsci-metabolomics` for implementation;
- use `ukbsci-plot` for manuscript-ready visualization.
