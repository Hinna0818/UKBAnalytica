# `ukbsci-proteomics` — function reference

All signatures verbatim from `R/Proteomics.R`.

---

## ID mapping

```r
protein_to_gene_symbol(proteins, protein_col = NULL,
                       from_type = "SYMBOL",
                       mapping_table = NULL,
                       mapping_protein_col = "protein",
                       mapping_symbol_col  = "gene_symbol",
                       organism_db = "org.Hs.eg.db",
                       drop_unmapped = TRUE)
```

**Returns:** `data.frame(protein, gene_symbol, mapping_source)`.

---

## ORA

```r
run_protein_ora(proteins, protein_col = NULL,
                from_type = "SYMBOL",
                mapping_table = NULL,
                mapping_protein_col = "protein",
                mapping_symbol_col  = "gene_symbol",
                universe = NULL,
                organism_db = "org.Hs.eg.db",
                ont          = "BP",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                pAdjustMethod = "BH",
                minGSSize = 10, maxGSSize = 500,
                readable = TRUE)
```

**Returns:** `list(gene_symbols, mapping, universe_symbols, ora_result)`.
`ora_result` is a `clusterProfiler::enrichResult`.

```r
run_protein_kegg_ora(proteins, protein_col = NULL,
                     from_type = "SYMBOL", mapping_table = NULL,
                     mapping_protein_col = "protein",
                     mapping_symbol_col  = "gene_symbol",
                     universe = NULL, organism_db = "org.Hs.eg.db",
                     organism = "hsa",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     pAdjustMethod = "BH",
                     minGSSize = 10, maxGSSize = 500,
                     readable = TRUE,
                     use_internal_data = FALSE)
```

**Returns:** `list(gene_symbols, entrez_ids, mapping, universe_symbols,
universe_entrez_ids, ora_result)`.

---

## PPI

```r
get_protein_ppi(proteins, protein_col = NULL,
                from_type = "SYMBOL",
                mapping_table = NULL,
                mapping_protein_col = "protein",
                mapping_symbol_col  = "gene_symbol",
                organism_db = "org.Hs.eg.db",
                taxID = 9606,
                required_score = NULL,
                network_type = "functional",
                add_nodes = 0,
                show_query_node_labels = 0,
                output = c("igraph", "data.frame"))
```

**Returns:** `list(source = "STRING", gene_symbols, mapping, ppi)`.

```r
subset_protein_ppi(ppi, n = NULL, score_cutoff = 0.7,
                   edge_attr = "score", rm_isolates = TRUE)

compute_protein_ppi_metrics(ppi, weight_attr = "score",
                            normalize = FALSE, seed = 42)

rank_protein_ppi_nodes(ppi,
                       metrics = c("degree","betweenness","closeness",
                                   "eccentricity","radiality","Stress",
                                   "MCC","MNC","DMNC","BN","EPC"),
                       weights = NULL, use_weight = TRUE, na_rm = TRUE)
```

Wraps `TCMDATA::ppi_subset`, `TCMDATA::compute_nodeinfo`,
`TCMDATA::rank_ppi_nodes`.

---

## Community detection

```r
run_protein_ppi_clustering(
  ppi,
  method = c("fastgreedy", "louvain", "mcode", "mcl"),
  ...
)
score_protein_ppi_clusters(ppi, cluster_attr = "fast_greedy_cluster",
                            min_size = 3)
```

---

## Robustness

```r
run_protein_ppi_robustness(ppi, targets, target_col = NULL,
                           from_type = "SYMBOL",
                           mapping_table = NULL,
                           mapping_protein_col = "protein",
                           mapping_symbol_col  = "gene_symbol",
                           organism_db = "org.Hs.eg.db",
                           n_perm = 100L,
                           weight_attr = "score",
                           rewire_niter = 10L,
                           seed = 42L)
```

**Returns:** `list(targets, mapping, robustness)`.

---

## Visualization

```r
plot_go_ora_bar(x, ...)
plot_enrichment_lollipop(x, ...)
```

Both accept either an `enrichResult` directly or a list with an `ora_result`
element.

---

## Suggests dependencies

`clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `igraph`, `TCMDATA`,
`xml2`.

---

## Privacy note

Gene symbols and protein IDs are non-identifying — they can be exported off
RAP. Per-participant Olink intensities (long format with `eid`) cannot.
