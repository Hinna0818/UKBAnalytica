---
name: ukbsci-metabolomics
description: >
  Metabolomics over-representation analysis (ORA) for UK Biobank Nightingale
  NMR metabolite data using UKBAnalytica. Covers metabolite classification
  (classify_metabolites), UK Biobank–to–MetaboAnalyst name mapping
  (metabolite_to_metaboanalyst_name), panel loading (load_ukb_metabolite_panel),
  ORA with a custom pathway library or MetaboAnalystR backend
  (run_metabolite_ora), and ggplot2 visualization of ORA results
  (plot_metabolite_ora_dotplot, plot_metabolite_ora_barplot). Use this skill
  when the user wants to interpret a list of significant metabolites from a
  regression or GWAS screen, identify enriched metabolic pathways, or
  visualize ORA results. Triggers: metabolomics, metabolite ORA, pathway
  enrichment, Nightingale, MetaboAnalyst, 代谢组, 代谢物富集,
  /ukbsci-metabolomics. Hard rule: local agents must not read or inspect real
  UKB RAP participant-level data; generate scripts for RAP execution and
  interpret aggregate outputs only.
---

# ukbsci-metabolomics — Metabolomics ORA for UK Biobank Nightingale data

## 0. RAP guardrails

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, or logs containing row-level values.
Generate scripts for the user to run inside RAP; only aggregate results
(ORA tables, rendered figures) may be shared back with the agent. See
`../references/agent-privacy-boundary.md`.

Input: a character vector of metabolite names (typically the significant
hits from a prior regression or association screen). Output: ORA result
table (pathway / hits / fold_enrichment / pvalue / p_adjust) and ggplot2
figures, which may be shared with the local agent.

---

## 1. When to load

- Interpret a significant metabolite list from regression, GWAS, or
  multi-omic analysis.
- Map UK Biobank Nightingale metabolite labels to MetaboAnalyst-compatible names.
- Classify metabolites (small molecules vs. lipoprotein-lipid measures vs.
  proteins) before pathway analysis.
- Run ORA against a custom pathway library or via MetaboAnalystR.
- Visualize enriched pathways as a dot plot or bar plot.

## 2. When NOT to load

- Regression of metabolites as exposures or outcomes → `ukbsci-regression`.
- Proteomics ORA → `ukbsci-proteomics`.
- Forest plots or other generic visualizations → `ukbsci-plot`.

---

## 3. Prerequisites

```r
library(UKBAnalytica)

# Optional: MetaboAnalystR backend (for smpdb_pathway, KEGG, etc.)
# remotes::install_github("xia-lab/MetaboAnalystR", ...)
```

The bundled metabolite panel is accessed automatically via `load_ukb_metabolite_panel()`.
No external files are required for the `backend = "custom"` workflow.

---

## 4. Pipeline

### Phase 1 — Inspect and classify the hit list

```r
# Load the bundled UK Biobank Nightingale metabolite panel (metadata only)
panel <- load_ukb_metabolite_panel()
head(panel)  # Description, UKB_ID, meta_ID columns

# Classify metabolites in your hit list
hits <- c("Alanine", "Glutamine", "Glycine", "Lactate", "Pyruvate",
          "LDL Cholesterol", "HDL Cholesterol", "Apolipoprotein B")

clf <- classify_metabolites(hits)
clf
# metabolite | category         | metaboanalyst_name | mapping_source
# Alanine    | small_molecule   | L-Alanine          | built_in
# LDL Cho... | lipoprotein_lipid| NA                 | NA
# Apolipopr..| protein          | NA                 | NA
```

**Categories:**
- `small_molecule` — amino acids, organic acids, ketones, etc.; eligible for ORA.
- `lipoprotein_lipid` — cholesterol fractions, triglycerides, lipoprotein particle
  counts, fatty acid composites; **excluded from small-molecule ORA by default**.
- `protein` — albumin, apolipoproteins; excluded.
- `unknown` — metabolite not matched to any built-in rule.

Only `small_molecule` metabolites are passed to `run_metabolite_ora()`.

---

### Phase 2 — Map metabolite names (optional)

```r
# Check and extend the built-in UKB → MetaboAnalyst name mapping
mapped <- metabolite_to_metaboanalyst_name(
  metabolites = hits,
  drop_unmapped = FALSE    # set TRUE to retain only successfully mapped names
)
# metabolite | metaboanalyst_name | mapping_source
```

To supply a custom mapping for non-standard names:

```r
custom_map <- data.frame(
  metabolite         = c("3-Hydroxyisobutyrate", "2-Hydroxyglutarate"),
  metaboanalyst_name = c("3-Hydroxyisobutyric acid", "2-Hydroxyglutaric acid")
)

mapped_ext <- metabolite_to_metaboanalyst_name(
  metabolites    = hits,
  mapping_table  = custom_map,
  drop_unmapped  = FALSE
)
```

Custom entries override the built-in map; `mapping_source` is set to
`"custom_table"` for those rows.

---

### Phase 3 — Run ORA

#### 3a. Custom pathway library (recommended)

Use this when you have a curated two-column pathway–metabolite library or
want full control over the background universe.

```r
pathway_db <- data.frame(
  pathway    = c(rep("Amino acid metabolism", 4),
                 rep("Energy metabolism",      3),
                 rep("Ketone body metabolism", 2)),
  metabolite = c("L-Alanine", "L-Glutamine", "Glycine", "L-Valine",
                 "Lactic acid", "Pyruvic acid", "D-Glucose",
                 "Acetoacetic acid", "3-Hydroxybutyric acid")
)

ora <- run_metabolite_ora(
  metabolites     = hits,
  pathway_db      = pathway_db,
  backend         = "custom",
  p_adjust_method = "BH",
  min_metabolites = 2      # minimum mapped hits required to test a pathway
)

# Inspect results
ora$matched      # mapped small-molecule names passed to ORA
ora$unmatched    # metabolites not mapped or not small-molecule
ora$ora_result   # data.frame: pathway, hits, pathway_size, fold_enrichment,
                 #             pvalue, p_adjust, neg_log10_p, hit_names

# Filter significant pathways
sig <- ora$ora_result[ora$ora_result$p_adjust < 0.05, ]
```

To specify a custom background universe (instead of all pathway-library metabolites):

```r
universe_metabolites <- panel$Description   # all 249 non-ratio metabolites

ora_custom_universe <- run_metabolite_ora(
  metabolites = hits,
  pathway_db  = pathway_db,
  universe    = universe_metabolites,
  backend     = "custom"
)
```

#### 3b. MetaboAnalystR backend (optional)

Requires MetaboAnalystR installed and internet access (downloads `.qs` library files).

```r
ora_ma <- run_metabolite_ora(
  metabolites    = hits,
  backend        = "metaboanalyst",
  library        = "smpdb_pathway",  # or "kegg_pathway", "hmdb_metabolites", etc.
  id_type        = "name",
  run_subprocess = TRUE   # runs MetaboAnalystR in a clean Rscript subprocess
                          # to avoid global-state side effects
)
```

`run_subprocess = TRUE` is the default and the safest option. Set `FALSE` only
if running inside an environment where subprocess calls are not available (e.g.
some HPC job submissions).

---

### Phase 4 — Visualize

```r
# Dot plot: fold enrichment (x-axis), pathways (y-axis),
#           point size = hits, colour = -log10(p)
plot_metabolite_ora_dotplot(
  ora,
  top_n = 15,
  p_col = "pvalue"
)

# Bar plot: -log10(p) bars, dashed line at p = 0.05
plot_metabolite_ora_barplot(
  ora,
  top_n = 15,
  fill_color = "#2F6FA3"
)
```

Both functions accept either a `ukb_metabolite_ora` object (the `run_metabolite_ora()`
output) or a plain data.frame (the `$ora_result` slot).

```r
# Save for manuscript
ggsave("/mnt/project/<area>/04-results/metabolomics_dotplot.pdf",
       plot_metabolite_ora_dotplot(ora), width = 6, height = 4)
```

---

## 5. Common pitfalls

1. **Lipoprotein measures are excluded from ORA.** `run_metabolite_ora()`
   restricts the query to `category == "small_molecule"`. If the user
   expects lipoprotein-related pathways, explain that these are not small
   molecules and are excluded by design.
2. **Unknown metabolites.** Non-standard or UKB-specific names may not match
   the built-in map. Use `metabolite_to_metaboanalyst_name()` with a custom
   mapping table, or check `ora$unmatched` and curate manually.
3. **`min_metabolites` too strict.** The default (`3`) drops pathways with fewer
   than 3 matched hits. Reduce to `2` for small hit lists; note this increases
   the risk of spurious enrichment.
4. **Universe choice.** The default custom-ORA universe is all metabolites in
   `pathway_db`. If the hit list was selected from a larger panel, pass the full
   panel as `universe` for a conservative, correctly calibrated background.
5. **MetaboAnalystR global state.** `run_subprocess = TRUE` (default) avoids
   side effects from MetaboAnalystR's use of global variables. Do not set
   `FALSE` unless in a constrained environment.
6. **Internet required for MetaboAnalystR backend.** Library `.qs` files are
   downloaded from MetaboAnalyst servers at runtime. Ensure outbound HTTP is
   available, or use `backend = "custom"` in air-gapped environments.
7. **Name normalisation.** The built-in name map normalises keys by lower-casing
   and removing non-alphanumeric characters, so minor spacing or capitalisation
   differences are tolerated. However, completely non-standard abbreviations
   (e.g. `"3OHB"` instead of `"3-Hydroxybutyrate"`) will not match.

---

## 6. Key functions

| Function | Returns |
|----------|---------|
| `load_ukb_metabolite_panel(file = NULL, file_encoding = "UTF-16LE")` | data.frame with `Description`, `UKB_ID`, `meta_ID` |
| `classify_metabolites(metabolites)` | data.frame: `metabolite`, `category`, `metaboanalyst_name`, `mapping_source` |
| `metabolite_to_metaboanalyst_name(metabolites, mapping_table = NULL, drop_unmapped = FALSE)` | data.frame: `metabolite`, `metaboanalyst_name`, `mapping_source` |
| `run_metabolite_ora(metabolites, pathway_db, universe, backend, library, mapping_table, min_metabolites, p_adjust_method, run_subprocess)` | `ukb_metabolite_ora` list: `input`, `mapping`, `matched`, `unmatched`, `ora_result`, `backend`, `library` |
| `plot_metabolite_ora_dotplot(x, top_n, p_col, size_col, pathway_col, color_low, color_high)` | ggplot2 object |
| `plot_metabolite_ora_barplot(x, top_n, p_col, pathway_col, fill_color)` | ggplot2 object |

**`ora_result` columns (custom backend):**

| Column | Description |
|--------|-------------|
| `pathway` | Pathway name |
| `hits` | Number of query metabolites in pathway |
| `pathway_size` | Total metabolites in pathway (within universe) |
| `query_size` | Total query metabolites used for ORA |
| `universe_size` | Background universe size |
| `expected` | Expected hits under the null |
| `fold_enrichment` | `hits / expected` |
| `pvalue` | Hypergeometric p-value |
| `p_adjust` | BH-adjusted p-value (or chosen method) |
| `neg_log10_p` | −log₁₀(pvalue) for plotting |
| `hit_names` | Semicolon-delimited list of matching metabolite names |

---

## 7. Related skills

| Skill | When |
|-------|------|
| `ukbsci-regression` | Regression of metabolites as exposures or outcomes. |
| `ukbsci-proteomics` | ORA for protein / omic-level data. |
| `ukbsci-plot` | Generic forest / volcano plots. |
| `ukbsci-workflow` | End-to-end orchestrator. |

---

## 8. References

- [`references/functions.md`](references/functions.md) — full argument tables
- [`references/examples.md`](references/examples.md) — end-to-end example
