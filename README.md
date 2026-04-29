# UKBAnalytica: Scalable Phenotyping and Statistical Pipeline for UK Biobank RAP Data

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/Hinna0818/UKBAnalytica_v2?style=flat)](https://github.com/Hinna0818/UKBAnalytica_v2/stargazers)
[![GitHub last commit](https://img.shields.io/github/last-commit/Hinna0818/UKBAnalytica_v2)](https://github.com/Hinna0818/UKBAnalytica_v2/commits/main)
[![Visits](https://hits.sh/github.com/Hinna0818/UKBAnalytica_v2.svg)](https://hits.sh/github.com/Hinna0818/UKBAnalytica_v2/)
<!-- badges: end -->

<img src="man/figures/logo.png" align="right" height="139" alt="UKBAnalytica logo" />

**UKBAnalytica** is a high-performance R package for processing UK Biobank 
Research Analysis Platform (RAP) data exports. It focuses on standardized
phenotyping, survival-ready datasets, scalable preprocessing, and downstream analysis.

**For details, please visit**: [Full documentation for UKBAnalytica](https://hinna0818.github.io/UKBAnalytica_v2/)

![](docs/image.png)


## Installation
You can install the development version of `UKBAnalytica` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Hinna0818/UKBAnalytica_v2")
```

Sometimes due to the network problem, it is not easy to use `devtools` to install, so you can install in this way:
```{r}
# install.packages("pak")
pak::pkg_install("Hinna0818/UKBAnalytica_v2")
```

Or just clone this repo and intall it locally:
```{bash}
git clone https://github.com/Hinna0818/UKBAnalytica_v2.git
cd UKBAnalytica
R CMD INSTALL .
```

## Quick start

```r
library(UKBAnalytica)
library(data.table)

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

head(analysis_dt[, .(
  eid,
  AA_history,
  Hypertension_history,
  Diabetes_history,
  outcome_status,
  outcome_surv_time
)])

# Optional: retrieve participant flow table printed by show_flow
flow_dt <- attr(analysis_dt, "participant_flow")
if (!is.null(flow_dt)) print(flow_dt)
```

## Phenotyping with ICD-10 + OPCS4

Predefined disease definitions now support operative procedure evidence through
`opcs4_pattern`. This is useful for surgical phenotypes and procedure-augmented
endpoints such as arrhythmia.

```r
rhythm_defs <- get_predefined_diseases()[
  c("Arrhythmia", "Atrial_Fibrillation", "Ventricular_Arrhythmia")
]

arrhythmia_dt <- build_survival_dataset(
  dt = ukb_data,
  disease_definitions = rhythm_defs,
  prevalent_sources = c("ICD10", "OPCS4"),
  outcome_sources = c("ICD10", "OPCS4"),
  primary_disease = "Arrhythmia"
)

head(arrhythmia_dt[, .(
  eid,
  Arrhythmia_history,
  Atrial_Fibrillation_history,
  Ventricular_Arrhythmia_history,
  outcome_status,
  outcome_surv_time
)])
```

If a disease definition does not provide `opcs4_pattern`, operative procedure
data are ignored even if `OPCS4` is included in `prevalent_sources` or
`outcome_sources`.

## Recent phenotyping updates (v0.6.2.1)

- Added `OPCS4` operative procedure support using `p41272` + `p41282_a*`.
- Added `opcs4_pattern` to `create_disease_definition()` for procedure-aware phenotyping.
- Added predefined arrhythmia endpoints: `Arrhythmia`, `Ventricular_Arrhythmia`, `AV_Block`, `Intraventricular_Block`, and `SVT`.
- Extended `Atrial_Fibrillation` to support combined ICD-10 + OPCS4 ascertainment.

## Recent survival updates (v0.6.2)

- Added optional `show_flow` in `build_survival_dataset()` to print participant attrition in terminal and attach a reusable flow table via `attr(result, "participant_flow")`.
- Added optional `dt_threads` in `build_survival_dataset()` to temporarily control `data.table` threads for large runs, with automatic restoration on exit.
- Added algorithm-source column compatibility for both `p{field}_i0` and `p{field}` naming styles.
- Improved date robustness in ICD/self-report/death parsing to prevent malformed date values from interrupting full-pipeline execution.

## What this package covers

### Core Functionality
- R-native RAP phenotype extraction via `dx extract_dataset` and `table-exporter`
- Legacy RAP data download helpers (Python scripts)
- Baseline preprocessing with standardized mappings
- Multi-source disease definitions (ICD-10, ICD-9, OPCS4, self-report, death, algorithm)
- Survival analysis datasets with prevalent/incident classification
- Baseline Table 1 summaries and multiple imputation

### Advanced Analysis Modules (v0.5.0+ supports)
- **Subgroup Analysis**: Stratified analysis with interaction p-values
- **Propensity Score Methods**: PSM matching and IPTW weighting
- **Mediation Analysis**: Causal mediation using regmedint backend
- **MI Pooling**: Multiple imputation result combining (Rubin's Rules)
- **Sensitivity Analysis Preprocessing**: Exclude early events or rows with missing covariates before regression
- **Machine Learning**: Unified ML interface with SHAP interpretation
- **Visualization**: Forest plots, K-M curves, balance plots, SHAP plots

## Machine Learning Example

```r
library(UKBAnalytica)

# Train a Random Forest classifier
ml_rf <- ukb_ml_model(
  diabetes ~ age + bmi + sbp + smoking,
  data = ukb_data,
  model = "rf",
  task = "classification",
  seed = 42
)

# Model evaluation
print(ml_rf)  # AUC, accuracy, etc.
ukb_ml_metrics(ml_rf, ci = TRUE)

# ROC curve
roc <- ukb_ml_roc(ml_rf)
plot(roc)

# SHAP interpretation (requires fastshap)
shap <- ukb_shap(ml_rf, sample_n = 1000)
plot_shap_summary(shap)
plot_shap_dependence(shap, feature = "age")

# Model comparison
ml_xgb <- ukb_ml_model(diabetes ~ ., data, model = "xgboost")
comparison <- ukb_ml_compare(ml_rf, ml_xgb)
plot(comparison)

# Survival ML
surv_rf <- ukb_ml_survival(
  Surv(time, event) ~ age + sex + bmi,
  data = ukb_data,
  model = "rsf"
)
print(surv_rf)  # C-index
```

## Advanced Analysis Example

```r
# Subgroup analysis
results <- run_subgroup_analysis(
  data = dt, exposure = "treatment", outcome = "event",
  subgroup_var = "age_group", model_type = "cox",
  endpoint = c("time", "status")
)
plot_forest(results)

# Multiple imputation pooling
pooled <- pool_mi_models(
  datasets = mi_datasets,
  formula = Surv(time, status) ~ treatment + age + sex,
  model_type = "cox"
)
summary(pooled)

# Mediation analysis
med <- run_mediation(
  data = dt, exposure = "treatment", mediator = "biomarker",
  outcome = "event", outcome_type = "cox"
)
plot_mediation(med, type = "effects")
```

## Sensitivity Analysis Example

```r
# Remove events occurring in the first 2 years of follow-up
dt_sens1 <- sensitivity_exclude_early_events(
  data = analysis_dt,
  endpoint = c("outcome_surv_time", "outcome_status"),
  n_years = 2
)

# Remove rows with any missing adjustment covariate
dt_sens2 <- sensitivity_exclude_missing_covariates(
  data = dt_sens1,
  covariates = c("age", "sex", "bmi", "smoking")
)

# Pass directly to the standard regression interface
cox_sens <- runmulti_cox(
  data = dt_sens2,
  main_var = c("bmi", "sbp"),
  covariates = c("age", "sex", "bmi", "smoking"),
  endpoint = c("outcome_surv_time", "outcome_status")
)
```

## Basic Workflow Demonstration for Chinese Users
请各位用户注意，**UKB的数据是不能下载到本地的**，这个R包开发的目的是提高研究人员在官方的RAP平台上的数据分析效率。当你在做一个UKB研究的时候，首先要获取你想要的数据（注意这里说的download指的是从云平台获取你研究相关的数据，不是下载到本地）。现在推荐在RAP的R会话中直接使用R函数提取表型数据：

```r
library(UKBAnalytica)

# 先检查项目里有哪些可用字段
fields <- rap_list_fields(pattern = "blood pressure|^participant\\.p31")

# 小规模同步提取，直接读回R
baseline <- rap_extract_pheno(
  variables = c("sex", "age", "bmi", "sbp_auto_1"),
  strip_entity_prefix = TRUE
)

# 大规模提取建议提交table-exporter云端任务
job <- rap_submit_extract(
  field_id = c(31, 53, 21022, 21001, 4080, 4079),
  file = "baseline_core"
)
job$job_id
```

`field_id`会提取这个UKB field下的全部instances/arrays；`variables`会按照`get_variable_info()`里的精确列名提取，更适合常见基线变量。`inst/python/ukb_data_loader.py`和`inst/python/protein_loader.py`仍然保留作为legacy/helper脚本，适合已经习惯Python命令行流程的用户。

其次，对于疾病诊断和前瞻性队列的生存时间计算，可以参考使用`build_survival_dataset()`函数，不同的sources对应疾病诊断来源，比如ICD、自我报告、算法定义、死亡、以及手术/操作编码OPCS4这些。`primary_disease`就是你的主疾病。参数`disease_definitions`就是一个疾病定义的格式，具体可以看文档（如果自己要定义的话）。如果一个疾病定义里没有设置`opcs4_pattern`，那么即使sources里写了`OPCS4`，也不会把手术信息纳入诊断。函数返回会有`xx_history`和`xx_incident`的列，history就是baseline及其前患病，如果是前瞻性队列就要去掉；incident对应前瞻性的事件发生，可以用来做cox回归。主疾病经过这个函数有两列是直接可以用于cox的，分别是`outcome_status`和`outcome_surv_time`，一个是事件是否发生0/1，另一个是生存时间。为什么要设置outcome sources和primary_sources？原因是有时候我们希望把基线患特定疾病作为协变量进行调整，这个函数就方便计算了。

后续的研究就是个性化的分析，常规的回归分析、生存分析、亚组分析、统计检验、机器学习等模块我们都有纳入。这个看自己情况来做即可。这个包的函数是很灵活的，希望大家好好挖掘，有建议可以提issue或者PR。

## Supplementary Materials
Here we provide some learning materials for UK Biobank in which you may be interested:
- [UK Biobank database browser](https://biobank.ndph.ox.ac.uk/ukb/index.cgi)
- [UK Biobank RAP platform](https://ukbiobank.dnanexus.com/landing)
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
