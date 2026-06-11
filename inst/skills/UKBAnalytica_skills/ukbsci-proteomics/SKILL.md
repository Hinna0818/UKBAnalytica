---
name: ukbsci-proteomics
description: >
  UK Biobank Olink / UKB-PPP proteomics post-hit analyses on the Research
  Analysis Platform (RAP) using UKBAnalytica. Covers protein-to-gene-symbol
  mapping (protein_to_gene_symbol), STRING PPI retrieval (get_protein_ppi),
  network metrics + node ranking (compute_protein_ppi_metrics,
  rank_protein_ppi_nodes), subsetting by confidence (subset_protein_ppi),
  community detection (run_protein_ppi_louvain, run_protein_ppi_fastgreedy,
  run_protein_ppi_mcl, run_protein_ppi_mcode + get_protein_mcode_res),
  cluster scoring (score_protein_ppi_clusters), robustness to target
  knockout (run_protein_ppi_robustness), GO and KEGG over-representation
  (run_protein_ora, run_protein_kegg_ora), and enrichment / network
  visualization (plot_go_ora_bar, plot_enrichment_lollipop). Use this skill
  when the user provides a list of significant proteins from UKB Olink hits
  and asks for PPI network analysis, gene-set enrichment, clustering, or
  network-robustness assessment. Triggers: UKB proteomics, Olink hits,
  UKB-PPP, STRING PPI, GO ORA, KEGG ORA, MCODE, Louvain, gene-set
  enrichment, 蛋白组分析, 通路富集, /ukbsci-proteomics. Hard rule: local agents must not read or inspect real UKB RAP participant-level data; generate scripts for RAP execution and interpret aggregate outputs only.
---

# ukbsci-proteomics — UKB Olink / UKB-PPP downstream analyses

## 0. RAP guardrails

Strict local-agent boundary: this skill is for script generation,
workflow planning, package guidance, and interpretation of aggregate outputs.
The agent must not read, inspect, summarize, or process real UK Biobank RAP
participant-level data, including de-identified row-level tables, raw RAP
fields, exact dates, per-row predictions, row-level SHAP matrices, screenshots,
tracebacks, or logs containing row-level values. Generate scripts for the user
to run inside RAP; only aggregate results or rendered figures may be shared
back with the agent. See `../references/agent-privacy-boundary.md`.

- Participant-level Olink intensities, protein matrices, and any tables that
  retain participant identifiers or row-level values stay in RAP.
- A **list of significant protein IDs / gene symbols** (no participant
  rows) is non-identifying and may be exported when needed (e.g. for
  cross-cohort replication, or for running STRING enrichment outside RAP
  if the user prefers).
- STRING / KEGG queries hit external servers — confirm RAP node has
  outbound internet **and** that the user is comfortable sending their
  gene-symbol hits.

---

## 1. When to load

- User has a list of significant proteins (e.g. from a Cox screen across
  ~3,000 Olink targets) and wants PPI / pathway downstream.
- User wants GO BP / MF / CC ORA or KEGG ORA.
- User wants to cluster the PPI network and characterize modules.
- User wants a "robustness" check: would the network collapse if specific
  proteins were knocked out?

## 2. When NOT to load

- Running the **upstream** protein-association screen (Cox / logistic for
  ~3,000 targets) → `ukbsci-regression`.
- Producing the volcano plot of effect estimates →
  `ukbsci-plot::plot_regression_volcano()`.

---

## 3. Prerequisites

```r
library(UKBAnalytica)
library(clusterProfiler)        # ORA + bitr
library(org.Hs.eg.db)           # human gene-symbol mapping
library(igraph)
library(TCMDATA)                # PPI helpers
library(data.table)
```

The protein input is typically the *significant subset* of a UKBAnalytica
Cox / logistic batch run — i.e. a `data.frame` with at least one column of
Olink protein identifiers and (optionally) effect sizes / adjusted p-values.

---

## 4. Pipeline (4 phases)

### Phase 1 — Map protein IDs to gene symbols

```r
hits <- fread("/mnt/project/<area>/04-results/protein_screen_sig.csv")
# columns: protein (Olink ID), HR, lower95, upper95, pvalue, padj

mapped <- protein_to_gene_symbol(
  hits,
  protein_col       = "protein",
  from_type         = "SYMBOL",            # Olink IDs are typically gene symbols
  drop_unmapped     = TRUE,
  organism_db       = "org.Hs.eg.db"
)
# Returns: data.frame(protein, gene_symbol, mapping_source)
genes <- unique(mapped$gene_symbol)
```

### Phase 2 — Enrichment (ORA)

**GO over-representation:**

```r
go_res <- run_protein_ora(
  proteins      = genes,
  from_type     = "SYMBOL",
  universe      = NULL,                    # NULL = whole-genome background
  ont           = "BP",                    # "BP", "MF", "CC", or "ALL"
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize     = 10, maxGSSize = 500,
  readable      = TRUE
)
go_res$ora_result                            # enrichResult

p_go <- plot_go_ora_bar(go_res)              # ggplot bar
```

**KEGG over-representation:**

```r
kegg_res <- run_protein_kegg_ora(
  proteins   = genes,
  organism   = "hsa",
  pvalueCutoff = 0.05,
  use_internal_data = FALSE                  # FALSE = live KEGG fetch
)
p_kegg <- plot_enrichment_lollipop(kegg_res) # ggplot lollipop
```

### Phase 3 — PPI network analysis

```r
ppi <- get_protein_ppi(
  proteins        = genes,
  taxID           = 9606,                    # Homo sapiens
  required_score  = 700,                     # STRING confidence × 1000
  network_type    = "functional",            # or "physical"
  output          = "igraph"
)
# ppi$ppi is an igraph object

# 3a. Confidence-based subset
ppi_top <- subset_protein_ppi(ppi$ppi, score_cutoff = 0.7, rm_isolates = TRUE)

# 3b. Compute network metrics
ppi_m <- compute_protein_ppi_metrics(ppi_top, weight_attr = "score",
                                      normalize = FALSE, seed = 42)

# 3c. Rank nodes by integrated centrality
ranked <- rank_protein_ppi_nodes(
  ppi_m,
  metrics = c("degree","betweenness","closeness",
              "eccentricity","radiality","Stress",
              "MCC","MNC","DMNC","BN","EPC"),
  use_weight = TRUE
)
ranked$table[1:20, ]
```

### Phase 4 — Community detection + cluster scoring

```r
# Pick ONE algorithm based on network size / preference
ppi_fg <- run_protein_ppi_clustering(
  ppi_m,
  method = "fastgreedy",
  n_clusters = 4,
  largest_component = TRUE
)

# Score clusters (size, density, average centrality)
scores <- score_protein_ppi_clusters(
  ppi_fg,
  cluster_attr = "fast_greedy_cluster",
  min_size     = 3
)
```

### Phase 5 — Robustness to target knockout

```r
robust <- run_protein_ppi_robustness(
  ppi          = ppi_m,
  targets      = c("APOB", "LPA"),         # known-causal candidates
  n_perm       = 100L,
  rewire_niter = 10L,
  seed         = 42L
)
robust$robustness          # TCMDATA result with permutation summary
```

### Phase 6 — Export

```r
fwrite(go_res$ora_result@result, "/mnt/project/<area>/04-results/12-go_ora.csv")
fwrite(kegg_res$ora_result@result, "/mnt/project/<area>/04-results/13-kegg_ora.csv")
fwrite(ranked$table,                "/mnt/project/<area>/04-results/14-ppi_node_rank.csv")
fwrite(scores,                       "/mnt/project/<area>/04-results/15-ppi_cluster_scores.csv")

ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig11-go_bar.pdf",
                p_go, width = 8, height = 6)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig12-kegg_lollipop.pdf",
                p_kegg, width = 8, height = 6)
```

---

## 5. Common pitfalls

1. **Olink ID format.** Most UKB Olink panel names are HGNC symbols (e.g.
   `IL6`, `CRP`), but a handful are not (`FAM3A`, `Q8NEW0`). When
   `drop_unmapped = TRUE`, those are silently removed; check
   `attr(mapped, "n_dropped")` if available.
2. **`universe = NULL` background.** Whole-genome background inflates
   enrichment. Use the union of *measured* Olink panel symbols as universe
   for honest results.
3. **`required_score = 700` is the STRING "high confidence" cutoff.**
   `400` is "medium". Don't go below `400`.
4. **Network-type choice.** `"physical"` returns only experimentally
   determined interactions; `"functional"` is broader (co-expression,
   text-mining). For mechanism-focused stories use `"physical"`.
5. **Largest connected component.** Fastgreedy requires it; the other
   algorithms tolerate disconnected components but typically cluster only
   the giant one. Note this in the manuscript.
6. **MCODE parameters.** `vwp` (vertex weight percentage) governs
   stringency; lower values yield more / smaller clusters.
7. **External-network calls.** `get_protein_ppi()` and
   `run_protein_kegg_ora(use_internal_data = FALSE)` make outbound HTTP
   calls. If the RAP node is network-restricted, set
   `use_internal_data = TRUE` for KEGG or supply your own STRING dump.
8. **`TCMDATA` dependency.** Many helpers (subset, rank, MCODE,
   robustness) thinly wrap `TCMDATA`. If the package is unavailable on the
   RAP node, fall back to manual `igraph` + `clusterProfiler` workflows.

---

## 6. Key functions

| Step | Function |
|------|----------|
| ID mapping | `protein_to_gene_symbol(proteins, protein_col, from_type, mapping_table, organism_db, drop_unmapped)` |
| GO ORA | `run_protein_ora(proteins, universe, ont, pvalueCutoff, qvalueCutoff, minGSSize, maxGSSize)` |
| KEGG ORA | `run_protein_kegg_ora(proteins, organism = "hsa", pvalueCutoff, use_internal_data)` |
| Get PPI | `get_protein_ppi(proteins, taxID, required_score, network_type, output)` |
| Subset | `subset_protein_ppi(ppi, n, score_cutoff, rm_isolates)` |
| Metrics | `compute_protein_ppi_metrics(ppi, weight_attr, normalize, seed)` |
| Rank | `rank_protein_ppi_nodes(ppi, metrics, use_weight, na_rm)` |
| Cluster | `run_protein_ppi_clustering(ppi, method = c("fastgreedy","louvain","mcode","mcl"), ...)` |
| Cluster score | `score_protein_ppi_clusters(ppi, cluster_attr, min_size)` |
| Robustness | `run_protein_ppi_robustness(ppi, targets, target_col, from_type, mapping_table, n_perm, rewire_niter, weight_attr, seed)` |
| Plot — GO bar | `plot_go_ora_bar(x)` |
| Plot — lollipop | `plot_enrichment_lollipop(x)` |

---

## 7. Related skills

| Skill | When |
|-------|------|
| `ukbsci-regression` | Run the protein-vs-outcome screen first. |
| `ukbsci-imputation` | Pool protein effects across imputations. |
| `ukbsci-plot::plot_regression_volcano` | Headline volcano of the screen. |
| `ukbsci-ml` | If the user wants ML on a protein feature set. |

---

## 8. References

- [`references/functions.md`](references/functions.md)
- [`references/examples.md`](references/examples.md)
