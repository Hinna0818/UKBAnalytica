# UKBAnalytica Skill Pack (`UKBAnalytica_skills`)

Agent-runtime-agnostic skill bundle for the [`UKBAnalytica`](https://github.com/Hinna0818/UKBAnalytica_v2)
R package. The pack works with any agent stack that supports a "skill" / "tool
description" loaded from YAML front-matter + Markdown (Claude Code skills,
OpenAI Assistants instructions / function-tool descriptions, custom RAG-driven
agents, etc.).

The pack lives at `skills/UKBAnalytica_skills/` in the repository root.
Load it via `file.path(getwd(), "skills/UKBAnalytica_skills")` or see
[INSTALL.md](INSTALL.md) for per-runtime loader snippets.

---

## Charter (read first)

All skills in this pack enforce the following non-negotiable rules:

1. **RAP-only execution.** Every code example assumes an authenticated session
   inside the UK Biobank Research Analysis Platform (RAP), typically
   `RAP JupyterLab → Terminal → R`. Skills never instruct the user to download
   UK Biobank participant-level raw data to a local laptop.
2. **RAP privacy boundary.** Do not export UK Biobank participant-level raw
   data, direct identifiers (e.g. `eid`, exact dates), or row-level source
   tables that can be linked back to participants. Aggregate results (model
   metrics, regression summaries, figures), and de-identified analytical
   outputs are generally acceptable for export. For individual-level
   visualizations (e.g. SHAP force plots), ensure no `eid`, raw RAP fields,
   exact dates, or re-identifiable row-level data are included.
3. **Large extracts go async.** When the requested field count is large,
   skills route the user through `rap_submit_extract()` (DNAnexus
   `table-exporter`) rather than pulling everything into the R session.
4. **`UKBAnalytica` is the sole executor.** Skills wrap real exported
   functions only — every function name in a skill must be present in
   `NAMESPACE`.
5. **No journal-brand styling.** Plotting guidance uses neutral palette names
   (`ukbsci_clinical`, `ukbsci_diverging`, `ukbsci_sequential`) and avoids
   referencing specific top-tier journal brands.

---

## Provider compatibility

| Runtime | Loader pattern |
|---------|----------------|
| Claude Code skills | Drop `skills/UKBAnalytica_skills/` into `~/.claude/skills/` (or workspace `.claude/skills/`); each `skills/<name>/SKILL.md` is auto-discovered via its `name` + `description` front-matter. |
| OpenAI Assistants / Responses | Concatenate `SKILL.md` files into the assistant `instructions` (or attach as file-search documents); use the `description` field as the routing summary. |
| LangChain / LlamaIndex agents | Treat each `skills/<name>/` as a Tool whose docstring = `description`, whose body = `SKILL.md` + `references/*.md`. |
| Generic JSON tool list | Read each `SKILL.md` YAML header → `{name, description}`; load body lazily when the agent decides to invoke. |

The front-matter is intentionally minimal (`name`, `description`) so it parses
identically across providers. No provider-specific keys are used.

---

## Pack layout

```
skills/UKBAnalytica_skills/
├── README.md                       ← this file
├── INSTALL.md                      ← per-provider install snippets
├── MANIFEST.json                   ← machine-readable index (name, description, path)
└── skills/
    ├── ukbsci-rap-extract/         (P2) RAP discover / plan / extract
    ├── ukbsci-cohort/              (P2) disease definitions + survival cohort
    ├── ukbsci-workflow/            (P2) end-to-end orchestrator
    ├── ukbsci-regression/          (P3) batch lm / logit / Cox + extensions
    ├── ukbsci-survival/            (P3) KM + risk table + log-rank
    ├── ukbsci-baseline/            (P3) tableone Table 1
    ├── ukbsci-propensity/          (P4) PS / PSM / IPTW / balance
    ├── ukbsci-mediation/           (P4) regmedint 4-way decomposition
    ├── ukbsci-subgroup-sensitivity/(P4) subgroup × interaction + sensitivity
    ├── ukbsci-imputation/          (P4) mice + Rubin pooling
    ├── ukbsci-proteomics/          (P5) Olink / STRING / GO / KEGG / PPI
    ├── ukbsci-ml/                  (P5) classification + survival ML + SHAP
    ├── ukbsci-preprocess/          (P5) variable cleaning + composites
    └── ukbsci-plot/                (P6) forest / volcano / calibration / theme
```

Each skill directory contains:

```
SKILL.md           ← YAML front-matter (name, description) + body
README.md          ← human-readable overview
references/
  functions.md     ← every exported function signature + caveats
  rap-guardrails.md← what is forbidden in this module
  examples.md      ← copy-pastable minimal examples
evals/
  evals.json       ← trigger-recall test cases for the skill router
```

---

## Trigger-word routing table

Every `SKILL.md` description ends with the phrase **"UK Biobank RAP" or "UKBAnalytica"**
to keep the router from confusing these skills with generic R / statistics
skills. Triggers below are written into the `description` field verbatim.

| Skill | English triggers | Chinese triggers |
|-------|------------------|------------------|
| `ukbsci-rap-extract` | UK Biobank RAP extract, dx extract_dataset, table-exporter, UKB field search, phenotype extraction | RAP 提取, 字段下载, table-exporter, UKB 字段搜索, 表型提取 |
| `ukbsci-cohort` | UKB cohort, disease definition, prevalent vs incident, survival dataset, ICD10 phenotyping | UKB 队列, 病例定义, prevalent vs incident, 生存数据集, ICD10 表型 |
| `ukbsci-workflow` | end-to-end UKB analysis, UKB pipeline, RAP-to-publication, full study plan | 端到端 UKB 分析, RAP 到论文, 完整流程, 项目计划 |
| `ukbsci-regression` | UKB regression, Cox model, logistic / linear, batch regression, PH diagnostics, competing risks, p_trend | UKB 回归, Cox 模型, 批量回归 |
| `ukbsci-survival` | UKB KM curve, Kaplan-Meier, log-rank, risk table | UKB 生存曲线, KM 曲线 |
| `ukbsci-baseline` | Table 1, baseline characteristics, demographics summary | 基线表, 基线特征 |
| `ukbsci-propensity` | propensity score, PSM, IPTW, ATE / ATT, Love plot, covariate balance | 倾向评分, 倾向得分匹配 |
| `ukbsci-mediation` | mediation, indirect effect, natural direct effect, TNIE, PNDE, proportion mediated | 中介分析 |
| `ukbsci-subgroup-sensitivity` | subgroup analysis, interaction, effect modification, sensitivity, complete-case, lag | 亚组分析, 敏感性分析 |
| `ukbsci-imputation` | multiple imputation, MI, mice, Rubin's rules, FMI | 多重插补, 插补合并 |
| `ukbsci-proteomics` | UKB proteomics, Olink, UKB-PPP, STRING PPI, GO ORA, KEGG ORA, MCODE | 蛋白组分析, 通路富集 |
| `ukbsci-ml` | UKB ML, XGBoost, random forest, SHAP, C-index, calibration, decision curve, AUC | UKB 机器学习, SHAP, 生存 ML |
| `ukbsci-preprocess` | UKB preprocessing, variable cleaning, negative code, derive BP / air pollution / diet score | UKB 变量预处理 |
| `ukbsci-plot` | UKB plotting, forest plot, volcano, calibration, manuscript figure | ukbsci 画图, 论文级图 |

---

## Status

All 14 skills shipped (v1.0.0). Phase grouping reflects the writing order:

| Phase | Skills | State |
|-------|--------|-------|
| P2 | `ukbsci-rap-extract`, `ukbsci-cohort`, `ukbsci-workflow` | shipped |
| P3 | `ukbsci-regression`, `ukbsci-survival`, `ukbsci-baseline` | shipped |
| P4 | `ukbsci-propensity`, `ukbsci-mediation`, `ukbsci-subgroup-sensitivity`, `ukbsci-imputation` | shipped |
| P5 | `ukbsci-proteomics`, `ukbsci-ml`, `ukbsci-preprocess` | shipped |
| P6 | `ukbsci-plot` | shipped |

See [`supp/UKBAnalytica-skill-roadmap.md`](../supp/UKBAnalytica-skill-roadmap.md)
for the full implementation roadmap.

---

## Citation

When agents emit R scripts based on this pack, the recommended citation header
at the top of every generated script is:

```r
###############################################################################
# UKBAnalytica Citation:
# He N. UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for
# UK Biobank RAP Data. R package version 0.6.2.2.
# https://github.com/Hinna0818/UKBAnalytica_v2
###############################################################################
```
