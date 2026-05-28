# Literature Search Strategies for UK Biobank Topic Scoping

Use these query patterns when `$ukb-research` needs a focused literature scan.
Adapt the terms to the user's disease, exposure, omics layer, and study design.

## Core PubMed query blocks

### Exposure-outcome epidemiology

```text
("UK Biobank"[Title/Abstract] OR "United Kingdom Biobank"[Title/Abstract])
AND (<exposure terms>)
AND (<outcome terms>)
AND (cohort OR prospective OR longitudinal OR incident)
```

### Phenotype and codelist methods

```text
("UK Biobank"[Title/Abstract])
AND (<disease terms>)
AND (ICD OR "self-report" OR "Hospital Episode Statistics" OR codelist OR phenotype)
```

### Omics profile or biomarker discovery

```text
("UK Biobank"[Title/Abstract])
AND (<disease terms>)
AND (proteomics OR metabolomics OR biomarkers OR Olink OR NMR)
AND (incident OR prospective OR prediction OR profile)
```

### Machine-learning prediction

```text
("UK Biobank"[Title/Abstract])
AND (<outcome terms>)
AND ("machine learning" OR prediction OR XGBoost OR "risk score")
AND (validation OR discrimination OR calibration OR SHAP)
```

### Mediation or causal pathway

```text
("UK Biobank"[Title/Abstract])
AND (<exposure terms>)
AND (<outcome terms>)
AND (mediation OR mediator OR pathway)
```

### Mendelian randomization / GWAS extension

```text
("UK Biobank"[Title/Abstract] OR GWAS)
AND (<exposure terms>)
AND (<outcome terms>)
AND ("Mendelian randomization" OR genetic instrument OR polygenic)
```

## Search expansion rules

- Use disease synonyms and abbreviations: e.g. COPD, chronic obstructive
  pulmonary disease; AF, atrial fibrillation.
- Search both exposure and platform names: e.g. air pollution, PM2.5, NO2,
  PM10, nitrogen oxides; proteomics, Olink, plasma proteins.
- Combine UKB-specific papers with one broader non-UKB search when biological
  plausibility or methods need support.
- Prefer papers with accessible supplementary methods, codelists, or code.

## Screening fields to extract

For each relevant paper, extract:

| Field | What to record |
|---|---|
| Citation | PMID, DOI, year, journal |
| Cohort | UKB only, multi-cohort, external validation |
| Population | sample size, age range, exclusions |
| Exposure/features | field IDs, assay, linkage, transformation |
| Outcome | disease source, ICD/self-report/death/cancer/GP, incident/prevalent |
| Model | Cox, logistic, linear, ML, MR, mediation, enrichment |
| Main finding | direction, magnitude, reproducibility |
| Limitations | what the new study can improve |
| Reusable elements | codelists, covariates, analysis scripts, software |

## Evidence ranking

Use this hierarchy when multiple references support the same claim:

1. UKB-specific prospective studies with clear phenotype definitions.
2. Multi-cohort studies that include UKB and report validation.
3. Method papers or codelist papers with reusable algorithms.
4. Mechanistic or clinical studies supporting biological plausibility.
5. Reviews, only for background context and not as sole support for a specific
   UKB design claim.
