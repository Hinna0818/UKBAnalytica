# `ukbsci-README.md` — Phase 8 wrap-up template

After the last specialist skill finishes, the workflow skill **must** write
this file at `/mnt/project/<area>/06-note/ukbsci-README.md`. The wrap-up is
fact-based: it lists the actual outputs and the headline numbers — never a
restatement of the plan.

The skill should:

1. Run `list.files()` against the project tree to enumerate real files.
2. Run `wc -l` against `07-log/*.log` to verify non-empty runs without
   copying log contents into the agent context.
3. Read the main coefficient table to extract the headline effect estimate.
4. Read the baseline table to extract Ns.
5. Generate the file from the template below.

---

```markdown
###############################################################################
# UKBAnalytica Citation:
# He N. UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for
# UK Biobank RAP Data. R package version <pkg_version>.
# https://github.com/Hinna0818/UKBAnalytica
###############################################################################

# <Study Title> — Analysis Summary

Final report generated <UTC timestamp>.

## 1. Analysis summary

- **Project**: <title>
- **Hypothesis**: <one-sentence claim>
- **Cohort size**: <n_total> participants enrolled; <n_prevalent> excluded as
  prevalent; <n_analytic> at-risk in the analytic dataset
- **Follow-up**: median <m> years (IQR <q1>–<q3>); maximum <max_y> years
- **Events**: <n_events> incident <primary_disease> cases
- **Primary model**: <description, e.g. "Multivariable Cox proportional hazards
  adjusting for age, sex, BMI, smoking, diabetes history">
- **Sensitivity variants run**: <list, e.g. "ICD-10 only, Algorithm-validated,
  excluded < 1y events">

## 2. Key findings

1. <Headline 1 — e.g. "Each 10-mmHg increase in baseline SBP was associated
   with HR 1.18 (95% CI 1.12–1.24) for incident AA.">
2. <Headline 2 — subgroup or sensitivity confirmation>
3. <Headline 3 — any unexpected pattern or null result>
4. <Headline 4 — ML model performance if Phase 7 ran>
5. <Caveat or limitation — e.g. "Algorithmic-outcome variant attenuated the
   estimate by ~12%; main results should be interpreted alongside the
   sensitivity table.">

## 3. Output tree (actual files produced)

\`\`\`
/mnt/project/<area>/
├── 01-script/
│   ├── 01-rap_extract.R
│   ├── 02-cohort_build.R
│   ├── ...
│   └── 11-plots.R
├── 02-extract/
│   ├── manifest.csv
│   └── pheno.csv                       (RAP-resident — DO NOT export)
├── 03-cohort/
│   ├── cohort.csv                       (RAP-resident)
│   ├── cohort_clean.csv                 (RAP-resident)
│   ├── cohort_flow.csv                  (aggregate — shareable)
│   └── source_compare.csv
├── 04-results/
│   ├── 01-baseline_table.csv
│   ├── 02-cox_main.csv
│   ├── 03-cox_subgroup.csv
│   ├── 04-cox_sensitivity.csv
│   ├── (...)
│   └── 08-ml_shap_summary.csv
├── 05-figs/
│   ├── Fig01-flow.pdf / .png
│   ├── Fig02-forest_main.pdf / .png
│   ├── Fig03-km.pdf / .png
│   ├── Fig04-forest_subgroup.pdf / .png
│   ├── Fig05-volcano.pdf / .png
│   ├── Fig06-shap_summary.pdf / .png
│   └── data/
│       ├── Fig02-forest_main.csv
│       ├── Fig03-km.csv
│       └── ...
├── 06-note/
│   ├── ukbsci-pipeline.md               (Phase 0 plan)
│   └── ukbsci-README.md                 (this file)
└── 07-log/
    ├── 01-rap_extract.log               (lines: <n>)
    ├── 02-cohort_build.log              (lines: <n>)
    └── ...
\`\`\`

## 4. Verification

- All scripts ran successfully (every log file non-empty): <yes/no>
- Citation header present in every script: <yes/no>
- Aggregate-only files identified for export: <count> files (see §3 markers)
- Decisions D0–D8 recorded in `06-note/ukbsci-pipeline.md`: <yes/no>

## 5. Reproducibility footer

\`\`\`r
sessionInfo()
\`\`\`

Save this footer with the wrap-up:

\`\`\`bash
Rscript -e 'capture.output(sessionInfo(), file = "/mnt/project/<area>/06-note/sessionInfo.txt")'
\`\`\`

## 6. Next steps

- [ ] Manuscript draft (use `04-results/01-baseline_table.csv` for Table 1)
- [ ] Submit figures (use `05-figs/*.pdf` editable sources)
- [ ] Archive scripts + manifest to project Git inside RAP
- [ ] Confirm aggregate-only export bundle with PI before any download
```

---

## How to populate from R

The workflow skill can generate this file programmatically:

```r
template <- readLines(
  system.file("skills/UKBAnalytica_skills/ukbsci-workflow/references/wrap-up-template.md",
              package = "UKBAnalytica")
)

# Trim the inline ```markdown wrapper, substitute placeholders, write out
body <- template[grep("^```markdown$", template)[1] + 1:length(template)]
body <- body[1:(grep("^```$", body)[1] - 1)]

# Substitute <title>, <pkg_version>, etc.
body <- gsub("<title>", "Baseline SBP and incident AA", body, fixed = TRUE)
body <- gsub("<pkg_version>",
             as.character(utils::packageVersion("UKBAnalytica")),
             body, fixed = TRUE)
# ... (other substitutions from saved Phase outputs)

writeLines(body, "/mnt/project/<area>/06-note/ukbsci-README.md")
```
