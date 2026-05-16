# `ukbsci-proteomics` — examples

```r
library(UKBAnalytica); library(clusterProfiler); library(org.Hs.eg.db)
library(igraph); library(TCMDATA); library(data.table)

hits <- fread("/mnt/project/<area>/04-results/protein_screen_sig.csv")
```

## A. GO + KEGG ORA from significant proteins

```r
mapped <- protein_to_gene_symbol(hits, protein_col = "protein",
                                  from_type = "SYMBOL")
genes <- unique(mapped$gene_symbol)

go_res <- run_protein_ora(genes, ont = "BP", pvalueCutoff = 0.05,
                          minGSSize = 10, maxGSSize = 500, readable = TRUE)
kegg_res <- run_protein_kegg_ora(genes, organism = "hsa", pvalueCutoff = 0.05)

ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig11-go_bar.pdf",
                plot_go_ora_bar(go_res), width = 8, height = 6)
ggplot2::ggsave("/mnt/project/<area>/05-figs/Fig12-kegg_lollipop.pdf",
                plot_enrichment_lollipop(kegg_res), width = 8, height = 6)

fwrite(go_res$ora_result@result,
       "/mnt/project/<area>/04-results/12-go_ora.csv")
fwrite(kegg_res$ora_result@result,
       "/mnt/project/<area>/04-results/13-kegg_ora.csv")
```

## B. PPI network + Louvain modules

```r
ppi <- get_protein_ppi(genes, taxID = 9606,
                       required_score = 700, network_type = "functional",
                       output = "igraph")
ppi_top <- subset_protein_ppi(ppi$ppi, score_cutoff = 0.7, rm_isolates = TRUE)
ppi_m   <- compute_protein_ppi_metrics(ppi_top, weight_attr = "score")
ranked  <- rank_protein_ppi_nodes(ppi_m,
            metrics = c("degree","betweenness","closeness","Stress",
                        "MCC","MNC","DMNC","BN","EPC"))
ppi_louv <- run_protein_ppi_louvain(ppi_m, resolution = 1.0)
scores   <- score_protein_ppi_clusters(ppi_louv,
              cluster_attr = "louvain_cluster", min_size = 3)
fwrite(ranked$table, "/mnt/project/<area>/04-results/14-ppi_node_rank.csv")
fwrite(scores,      "/mnt/project/<area>/04-results/15-ppi_cluster_scores.csv")
```

## C. MCODE clusters + robustness to knockout

```r
ppi_mc <- run_protein_ppi_mcode(ppi_m, vwp = 0.2, degree_cutoff = 2)
mc_df  <- get_protein_mcode_res(ppi_mc, only_clusters = FALSE)

robust <- run_protein_ppi_robustness(
  ppi = ppi_m, targets = c("APOB","LPA"),
  n_perm = 100L, rewire_niter = 10L, seed = 42L
)
```
