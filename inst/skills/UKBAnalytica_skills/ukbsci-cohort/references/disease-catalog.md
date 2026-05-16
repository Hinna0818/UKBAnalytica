# Pre-validated disease catalog (`get_predefined_diseases()`)

**Canonical authority:** always call `get_predefined_diseases()` in the live
R session for authoritative code patterns and source field IDs. This catalog
is a quick-glance map only — code patterns and field lists can change between
package versions.

```r
defs <- get_predefined_diseases()
names(defs)           # see all 76 keys
str(defs$COPD)        # inspect any definition
```

---

## How to read the source-types column

| Abbreviation | Meaning |
|-------------|---------|
| ICD-10/9 | Hospital inpatient (HES) codes |
| SR | Self-reported illness (touchscreen) |
| Death | Linked death registry ICD-10 |
| FO | First Occurrence date fields (Category 1712) |
| Algo | UKB algorithmically-defined outcome (Category 42) |
| OPCS4 | Operative procedure codes |
| CR | Cancer registry (Fields 40005/40006/40011/40012) |

The table shows the source types used in the curated definition. If a source
type is listed, `create_disease_definition()` includes the relevant fields;
inspect the live object for exact patterns.

---

## 1. Cardiovascular and coronary disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `CVD` | Cardiovascular disease | ICD-10/9, SR, FO |
| `MI` | Myocardial infarction | ICD-10/9, SR, FO, Algo |
| `STEMI` | ST-elevation MI | Algo |
| `NSTEMI` | Non-ST MI | Algo |
| `HF` | Heart failure | ICD-10/9, SR |
| `Angina` | Angina pectoris | ICD-10/9, SR |

## 2. Cerebrovascular disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Stroke` | Stroke (all subtypes) | ICD-10/9, SR, FO, Algo |
| `Ischaemic_Stroke` | Ischaemic stroke | ICD-10/9, Algo |
| `Intracerebral_Haemorrhage` | ICH | ICD-10/9 |
| `Subarachnoid_Haemorrhage` | SAH | ICD-10/9 |
| `Stroke_TIA` | Stroke or TIA | ICD-10/9 |

## 3. Cardiac rhythm and conduction

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Arrhythmia` | Any arrhythmia | ICD-10/9 |
| `Atrial_Fibrillation` | Atrial fibrillation | ICD-10/9 |
| `Ventricular_Arrhythmia` | Ventricular arrhythmia | ICD-10/9 |
| `AV_Block` | AV block | ICD-10/9 |
| `Intraventricular_Block` | Intraventricular block | ICD-10/9 |
| `SVT` | Supraventricular tachycardia | ICD-10/9 |

## 4. Vascular and thromboembolic disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `AA` | Aortic aneurysm | ICD-10/9, FO |
| `TAA` | Thoracic aortic aneurysm | ICD-10/9 |
| `AAA` | Abdominal aortic aneurysm | ICD-10/9 |
| `Vascular_Disease` | Composite vascular disease | ICD-10/9 |
| `PAD` | Peripheral artery disease | ICD-10/9 |
| `VTE` | Venous thromboembolism | ICD-10/9 |
| `Hypertension` | Hypertension | ICD-10/9, SR |

## 5. Metabolic and endocrine disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Diabetes` | Diabetes (any type) | ICD-10/9, SR |
| `T1DM` | Type 1 diabetes | ICD-10/9, SR |
| `T2DM` | Type 2 diabetes | ICD-10/9, SR |
| `Hyperlipidemia` | Hyperlipidaemia | ICD-10/9, SR |
| `Thyroid_Disorders` | Thyroid disorders | ICD-10/9, SR |

## 6. Respiratory disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Asthma` | Asthma | ICD-10/9, SR, Algo |
| `COPD` | Chronic obstructive pulmonary disease | ICD-10/9, SR, Algo |
| `Bronchiectasis` | Bronchiectasis | ICD-10/9 |

## 7. Renal disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `CKD` | Chronic kidney disease | ICD-10/9 |
| `ESRD` | End-stage renal disease | ICD-10/9 |

## 8. Cancer

For cancer phenotypes, the cancer registry (`CR`) is the preferred primary
ascertainment source when it is available and the research question focuses on
incident malignancy. ICD-10 hospital codes and death registry entries serve
as complementary sources depending on the scientific question.

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Lung_Cancer` | Lung cancer | ICD-10, CR |
| `Breast_Cancer` | Breast cancer | ICD-10, CR |
| `Prostate_Cancer` | Prostate cancer | ICD-10, CR |
| `Colorectal_Cancer` | Colorectal cancer | ICD-10, CR |
| `Melanoma` | Melanoma | ICD-10, CR |
| `Non_Melanoma_Skin_Cancer` | Non-melanoma skin cancer | ICD-10, CR |
| `Ovarian_Cancer` | Ovarian cancer | ICD-10, CR |
| `Oesophageal_Cancer` | Oesophageal cancer | ICD-10, CR |
| `Stomach_Cancer` | Gastric cancer | ICD-10, CR |

### Cancer histology and behaviour codes (for custom `create_disease_definition()`)

Use `cancer_behaviour = 3L` (malignant) for almost all incident cancer
studies. Histology codes follow ICD-O-3 morphology (UKB Field 40011).

| Code range | Morphology |
|------------|-----------|
| 8070–8079 | Squamous cell carcinoma |
| 8140–8149 | Adenocarcinoma |
| 8240–8249 | Carcinoid / neuroendocrine |
| 8500–8509 | Ductal carcinoma (breast) |
| 8520–8529 | Lobular carcinoma (breast) |
| 8720–8790 | Melanoma family |

## 9. Neurological disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Parkinsons` | Parkinson's disease | ICD-10/9, SR, Algo |
| `Parkinsonism` | Parkinsonism (any cause) | ICD-10/9 |
| `Progressive_Supranuclear_Palsy` | PSP | ICD-10/9 |
| `Multiple_System_Atrophy` | MSA | ICD-10/9 |
| `Dementia` | Dementia (all-cause) | ICD-10/9, SR, Algo |
| `Alzheimers_Disease` | Alzheimer's disease | ICD-10/9, Algo |
| `Vascular_Dementia` | Vascular dementia | ICD-10/9 |
| `Frontotemporal_Dementia` | FTD | ICD-10/9 |
| `Motor_Neurone_Disease` | MND / ALS | ICD-10/9 |
| `Multiple_Sclerosis` | Multiple sclerosis | ICD-10/9, SR |
| `Epilepsy` | Epilepsy | ICD-10/9, SR |
| `Migraine` | Migraine | ICD-10/9, SR |

## 10. Mental health and substance use

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Depression` | Depression | ICD-10/9, SR |
| `Anxiety` | Anxiety disorders | ICD-10/9, SR |
| `Schizophrenia_Bipolar` | Schizophrenia / bipolar | ICD-10/9 |
| `Alcohol_Use_Disorder` | Alcohol use disorder | ICD-10/9 |
| `Substance_Use_Disorder` | Substance use disorder | ICD-10/9 |

## 11. Gastrointestinal and liver disease

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Dyspepsia` | Dyspepsia | ICD-10/9, SR |
| `Irritable_Bowel_Syndrome` | IBS | ICD-10/9, SR |
| `Inflammatory_Bowel_Disease` | IBD | ICD-10/9, SR |
| `Diverticular_Disease` | Diverticular disease | ICD-10/9 |
| `Treated_Constipation` | Treated constipation | ICD-10/9, SR |
| `Chronic_Liver_Disease` | Chronic liver disease | ICD-10/9 |

## 12. Musculoskeletal and rheumatological

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Osteoarthritis` | Osteoarthritis | ICD-10/9, SR |
| `Rheumatoid_Arthritis` | Rheumatoid arthritis | ICD-10/9, SR |
| `Fracture` | Fracture | ICD-10/9, OPCS4 |

## 13. Sensory and ophthalmological

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Menieres_Disease` | Ménière's disease | ICD-10/9 |
| `Glaucoma` | Glaucoma | ICD-10/9, SR |
| `Cataract` | Cataract | ICD-10/9, SR |
| `AMD` | Age-related macular degeneration | ICD-10/9 |

## 14. Dermatological and other chronic conditions

| Key | Name | Primary source types |
|-----|------|----------------------|
| `Pernicious_Anaemia` | Pernicious anaemia | ICD-10/9, SR |
| `Psoriasis_Eczema` | Psoriasis / eczema | ICD-10/9, SR |
| `Prostate_Disorders` | Prostate disorders | ICD-10/9, SR |

---

## First Occurrence category field IDs (quick reference)

First Occurrence fields are in UKB Category 1712. Common examples:

| Region | Example date field IDs |
|--------|------------------------|
| IHD (I20–I25) | 131296, 131298, 131300, 131302, 131304, 131306 |
| Cerebrovascular (I60–I69) | 131360, 131362, 131364, 131366, 131368 |
| Aorta / PAD (I70–I79) | 131380, 131382, 131384, 131386 |
| Diabetes (E10–E14) | 130706, 130708, 130710, 130712, 130714 |
| COPD (J40–J47) | 131490, 131492, 131494, 131496 |
| Dementia (F00–F03, G30) | 130836, 130838, 130840, 131036 |
| Parkinson's (G20–G22) | 131022, 131024, 131026 |
| CKD (N17–N19) | 132036, 132038, 132040 |

The source field ID is `date_field_id + 1` unless overridden via
`first_occurrence_source_fields` in `create_disease_definition()`.

---

## Sentinel date values excluded by the parsers

The First Occurrence and Algorithm parsers automatically drop these
UKB sentinel dates — they are not real event dates:

| Sentinel date | Meaning |
|---------------|---------|
| `1900-01-01` | Unknown / unverifiable date (Algorithm source) |
| `1901-01-01` | Before available linkage start (First Occurrence) |
| `1902-02-02`, `1903-03-03` | System-level placeholders |
| `2037-07-07` | After available linkage end |

Do not re-add these values to analysis data.

---

## Catalog maintenance note

`get_predefined_diseases()` is updated when:
- UKB Showcase releases new Category 1712 First Occurrence fields
- UKB releases new Category 42 algorithmically-defined outcomes
- Phenotyping consortia (CALIBER, OpenSAFELY) publish updated mappings

When the user's study protocol differs from a curated definition, inspect
the live definition (`str(get_predefined_diseases()[[key]])`) and build a
custom definition with `create_disease_definition()` rather than assuming
the curated version is correct.
