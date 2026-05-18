# Function reference — ukbsci-metabolomics

## `load_ukb_metabolite_panel()`

```r
load_ukb_metabolite_panel(
  file          = NULL,           # path to custom panel file; NULL = bundled file
  file_encoding = "UTF-16LE"      # encoding of the panel file
)
```

Returns a data.frame with columns: `Description`, `UKB_ID`, `meta_ID` (and
any others present in the panel file). Used to inspect available metabolite names
and field IDs for the UK Biobank Nightingale non-ratio panel.

---

## `classify_metabolites()`

```r
classify_metabolites(metabolites)
# metabolites: character vector of metabolite names
```

Returns a data.frame:

| Column | Description |
|--------|-------------|
| `metabolite` | Cleaned input name |
| `category` | `"small_molecule"`, `"lipoprotein_lipid"`, `"protein"`, or `"unknown"` |
| `metaboanalyst_name` | Mapped MetaboAnalyst name (NA if no match) |
| `mapping_source` | `"built_in"`, `"custom_table"`, or NA |

Classification logic: checks built-in name map first (→ `small_molecule`);
then pattern-matches on `Albumin`, `Apolipoprotein` (→ `protein`); then
lipoprotein/lipid keyword list (→ `lipoprotein_lipid`); otherwise `"unknown"`.

---

## `metabolite_to_metaboanalyst_name()`

```r
metabolite_to_metaboanalyst_name(
  metabolites,
  mapping_table          = NULL,   # custom data.frame with 2+ columns
  mapping_metabolite_col = "metabolite",
  mapping_name_col       = "metaboanalyst_name",
  drop_unmapped          = FALSE
)
```

Returns a data.frame with `metabolite`, `metaboanalyst_name`, `mapping_source`.
Keys are normalised (lower-case, non-alphanumerics removed) for fuzzy matching.
Custom table overrides built-in map for matching keys.

---

## `run_metabolite_ora()`

```r
run_metabolite_ora(
  metabolites,
  pathway_db      = NULL,         # required for backend = "custom"
  universe        = NULL,         # background vector; NULL = pathway-library metabolites
  backend         = c("custom", "metaboanalyst"),
  id_type         = "name",       # MetaboAnalystR identifier type
  library         = "smpdb_pathway",  # MetaboAnalystR library name
  mapping_table   = NULL,         # custom mapping table (passed to name mapper)
  pathway_col     = "pathway",    # column in pathway_db
  metabolite_col  = "metabolite", # column in pathway_db
  min_metabolites = 3,            # minimum mapped hits to run ORA
  p_adjust_method = "BH",         # method for stats::p.adjust()
  run_subprocess  = TRUE          # MetaboAnalystR: clean subprocess
)
```

Returns an object of class `ukb_metabolite_ora` (a list):

| Slot | Content |
|------|---------|
| `input` | Cleaned input metabolite names |
| `mapping` | data.frame from `metabolite_to_metaboanalyst_name()` with `category` |
| `matched` | Unique mapped small-molecule names used as query |
| `unmatched` | Metabolites not mapped or not small-molecule |
| `ora_result` | ORA result data.frame (see Phase 3 column table in SKILL.md) |
| `backend` | `"custom"` or `"metaboanalyst"` |
| `library` | Library name used (`"custom"` for custom backend) |
| `p_adjust_method` | Adjustment method used |

---

## `plot_metabolite_ora_dotplot()`

```r
plot_metabolite_ora_dotplot(
  x,                        # ukb_metabolite_ora object or ora_result data.frame
  top_n       = 15,
  p_col       = "pvalue",
  size_col    = "hits",
  pathway_col = "pathway",
  color_low   = "#2F6FA3",
  color_high  = "#C74732"
)
```

Returns a ggplot2 object. X-axis = fold enrichment, y-axis = top pathways ordered
by p-value; point size = hits, colour gradient = −log₁₀(p).

---

## `plot_metabolite_ora_barplot()`

```r
plot_metabolite_ora_barplot(
  x,
  top_n       = 15,
  p_col       = "pvalue",
  pathway_col = "pathway",
  fill_color  = "#2F6FA3"
)
```

Returns a ggplot2 object. Y-axis = top pathways, x-axis = −log₁₀(pvalue).
Dashed line at −log₁₀(0.05).
