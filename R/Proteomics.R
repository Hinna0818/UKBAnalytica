#' @importFrom igraph simplify as_undirected components induced_subgraph vcount
#'   ecount cluster_fast_greedy cut_at set_vertex_attr
NULL

#' Convert protein identifiers to gene symbols
#'
#' @description
#' Convert a vector of protein identifiers into HGNC gene symbols for downstream
#' enrichment analysis. When a custom mapping table is supplied, it is used
#' first. Remaining unmatched identifiers can then be mapped with
#' `clusterProfiler::bitr()`. Inputs in UK Biobank Olink coding 143 format,
#' such as `"IL6;Interleukin-6"`, and RAP-exported Olink column names such as
#' `"olink_instance_0.eno2"` are parsed automatically. Multi-target Olink
#' symbols such as `"IL12A_IL12B"` are expanded into one row per gene symbol.
#'
#' @param proteins A character vector of protein identifiers, or a data.frame
#'   containing a protein identifier column.
#' @param protein_col Optional column name when `proteins` is a data.frame.
#' @param from_type Character string. Identifier type used by Bioconductor when
#'   `mapping_table` does not fully resolve the input. Default is `"SYMBOL"`.
#' @param input_format Input naming convention. `"auto"` parses only recognized
#'   RAP Olink column names and Olink semicolon coding; `"rap_olink"` explicitly
#'   extracts the suffix after the final dot; `"symbol"` validates gene-symbol
#'   style identifiers; and `"raw"` preserves identifiers apart from trimming
#'   whitespace and converting to upper case.
#' @param mapping_table Optional data.frame containing custom protein-to-symbol
#'   mappings.
#' @param mapping_protein_col Column name in `mapping_table` containing protein
#'   identifiers. Default is `"protein"`.
#' @param mapping_symbol_col Column name in `mapping_table` containing gene
#'   symbols. Default is `"gene_symbol"`.
#' @param organism_db Character string naming the OrgDb package. Default is
#'   `"org.Hs.eg.db"`.
#' @param drop_unmapped Logical. If `TRUE`, drop rows without a mapped gene
#'   symbol. Default is `TRUE`.
#'
#' @return A data.frame with columns `protein`, `gene_symbol`, and
#'   `mapping_source`. Its `mapping_summary` attribute reports input, mapped,
#'   and unmapped protein counts.
#'
#' @export
protein_to_gene_symbol <- function(proteins,
                                   protein_col = NULL,
                                   from_type = "SYMBOL",
                                   mapping_table = NULL,
                                   mapping_protein_col = "protein",
                                   mapping_symbol_col = "gene_symbol",
                                   organism_db = "org.Hs.eg.db",
                                   drop_unmapped = TRUE,
                                   input_format = c("auto", "rap_olink", "symbol", "raw")) {
  from_type <- toupper(from_type)
  input_format <- match.arg(input_format)
  protein_vec <- .extract_protein_identifiers(proteins, protein_col)
  protein_vec <- .normalize_identifiers(protein_vec, input_format = input_format)

  if (!length(protein_vec)) {
    stop("No valid protein identifiers were provided.")
  }

  result <- data.frame(
    protein = protein_vec,
    gene_symbol = NA_character_,
    mapping_source = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(mapping_table)) {
    if (!is.data.frame(mapping_table)) {
      stop("`mapping_table` must be a data.frame.")
    }

    required_cols <- c(mapping_protein_col, mapping_symbol_col)
    if (!all(required_cols %in% names(mapping_table))) {
      stop(
        sprintf(
          "`mapping_table` must contain columns: %s.",
          paste(required_cols, collapse = ", ")
        )
      )
    }

    map_df <- data.frame(
      protein = .normalize_identifier_rows(
        mapping_table[[mapping_protein_col]],
        input_format = input_format
      ),
      gene_symbol = .normalize_gene_symbols(mapping_table[[mapping_symbol_col]]),
      stringsAsFactors = FALSE
    )
    map_df <- map_df[
      !is.na(map_df$protein) & nzchar(map_df$protein) &
        !is.na(map_df$gene_symbol) & nzchar(map_df$gene_symbol),
      ,
      drop = FALSE
    ]
    map_df <- unique(map_df)

    if (nrow(map_df)) {
      input_order <- match(map_df$protein, result$protein)
      custom_hits <- map_df[!is.na(input_order), , drop = FALSE]
      custom_hits$input_order <- input_order[!is.na(input_order)]
      custom_hits$mapping_source <- "custom_table"
      custom_hits <- custom_hits[
        order(custom_hits$input_order),
        c("protein", "gene_symbol", "mapping_source"),
        drop = FALSE
      ]

      unmatched <- result[!result$protein %in% custom_hits$protein, , drop = FALSE]
      result <- rbind(custom_hits, unmatched)
    }
  }

  remaining <- is.na(result$gene_symbol) & !is.na(result$protein)
  if (any(remaining)) {
    result <- .map_remaining_to_symbol(
      result = result,
      from_type = from_type,
      organism_db = organism_db
    )
  }

  result$gene_symbol <- .normalize_gene_symbols(result$gene_symbol)
  result <- .expand_multi_target_symbols(result)

  if (drop_unmapped) {
    result <- result[!is.na(result$gene_symbol) & nzchar(result$gene_symbol), , drop = FALSE]
  }

  result <- unique(result)
  rownames(result) <- NULL
  attr(result, "mapping_summary") <- list(
    input = length(protein_vec),
    mapped = length(unique(result$protein[!is.na(result$gene_symbol)])),
    unmapped = length(setdiff(protein_vec, result$protein[!is.na(result$gene_symbol)]))
  )
  result
}


#' Run GO ORA enrichment for proteomics hits
#'
#' @description
#' Convert protein identifiers to gene symbols and run over-representation
#' analysis (ORA) with `clusterProfiler::enrichGO()`. This is the GO-specific
#' interface for proteomics hits extracted from UK Biobank RAP Olink data.
#'
#' @param proteins A character vector of protein identifiers, or a data.frame
#'   containing a protein identifier column.
#' @param protein_col Optional column name when `proteins` is a data.frame.
#' @param from_type Character string describing the input identifier type for
#'   Bioconductor-based mapping. Default is `"SYMBOL"`.
#' @param input_format Input naming convention passed to
#'   [protein_to_gene_symbol()].
#' @param mapping_table Optional data.frame containing custom protein-to-symbol
#'   mappings.
#' @param mapping_protein_col Column name in `mapping_table` containing protein
#'   identifiers. Default is `"protein"`.
#' @param mapping_symbol_col Column name in `mapping_table` containing gene
#'   symbols. Default is `"gene_symbol"`.
#' @param universe Optional character vector of background protein identifiers.
#'   These identifiers are converted with the same rules as `proteins`.
#' @param organism_db Character string naming the OrgDb package. Default is
#'   `"org.Hs.eg.db"`.
#' @param ont One of `"BP"`, `"MF"`, `"CC"`, or `"ALL"`. Passed to
#'   `clusterProfiler::enrichGO()`.
#' @param pvalueCutoff Numeric p-value cutoff for ORA. Default is `0.05`.
#' @param qvalueCutoff Numeric q-value cutoff for ORA. Default is `0.2`.
#' @param pAdjustMethod Character string for multiple-testing correction.
#'   Default is `"BH"`.
#' @param minGSSize Minimum gene set size. Default is `10`.
#' @param maxGSSize Maximum gene set size. Default is `500`.
#' @param readable Logical. Passed to `clusterProfiler::enrichGO()`. Default is
#'   `TRUE`.
#'
#' @return A list with components `gene_symbols`, `mapping`, `mapping_summary`,
#'   `universe_symbols`, and `ora_result`.
#'
#' @export
run_protein_ora <- function(proteins,
                            protein_col = NULL,
                            from_type = "SYMBOL",
                            mapping_table = NULL,
                            mapping_protein_col = "protein",
                            mapping_symbol_col = "gene_symbol",
                            universe = NULL,
                            organism_db = "org.Hs.eg.db",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2,
                            pAdjustMethod = "BH",
                            minGSSize = 10,
                            maxGSSize = 500,
                            readable = TRUE,
                            input_format = c("auto", "rap_olink", "symbol", "raw")) {
  .require_package("clusterProfiler")
  input_format <- match.arg(input_format)

  enrich_input <- .prepare_protein_enrichment_input(
    proteins = proteins,
    protein_col = protein_col,
    from_type = from_type,
    input_format = input_format,
    mapping_table = mapping_table,
    mapping_protein_col = mapping_protein_col,
    mapping_symbol_col = mapping_symbol_col,
    universe = universe,
    organism_db = organism_db,
    target_type = "SYMBOL"
  )

  orgdb_obj <- .get_orgdb_object(organism_db)
  ora_result <- clusterProfiler::enrichGO(
    gene = enrich_input$gene_symbols,
    universe = enrich_input$universe_symbols,
    OrgDb = orgdb_obj,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    readable = readable
  )

  list(
    database = "GO",
    gene_symbols = enrich_input$gene_symbols,
    mapping = enrich_input$mapping,
    mapping_summary = enrich_input$mapping_summary,
    universe_symbols = enrich_input$universe_symbols,
    ora_result = ora_result
  )
}


#' Run KEGG ORA enrichment for proteomics hits
#'
#' @description
#' Convert protein identifiers to gene symbols, then to Entrez IDs, and run
#' over-representation analysis (ORA) with `clusterProfiler::enrichKEGG()`.
#'
#' @inheritParams run_protein_ora
#' @param organism Character string for KEGG organism code. Default is `"hsa"`.
#' @param use_internal_data Logical. Passed to `clusterProfiler::enrichKEGG()`.
#'   Default is `FALSE`.
#'
#' @return A list with components `gene_symbols`, `entrez_ids`, `mapping`,
#'   `mapping_summary`, `universe_symbols`, `universe_entrez_ids`, and
#'   `ora_result`.
#'
#' @export
run_protein_kegg_ora <- function(proteins,
                                 protein_col = NULL,
                                 from_type = "SYMBOL",
                                 mapping_table = NULL,
                                 mapping_protein_col = "protein",
                                 mapping_symbol_col = "gene_symbol",
                                 universe = NULL,
                                 organism_db = "org.Hs.eg.db",
                                 organism = "hsa",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.2,
                                 pAdjustMethod = "BH",
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 readable = TRUE,
                                 use_internal_data = FALSE,
                                 input_format = c("auto", "rap_olink", "symbol", "raw")) {
  .require_package("clusterProfiler")
  input_format <- match.arg(input_format)

  enrich_input <- .prepare_protein_enrichment_input(
    proteins = proteins,
    protein_col = protein_col,
    from_type = from_type,
    input_format = input_format,
    mapping_table = mapping_table,
    mapping_protein_col = mapping_protein_col,
    mapping_symbol_col = mapping_symbol_col,
    universe = universe,
    organism_db = organism_db,
    target_type = "ENTREZID"
  )

  ora_result <- clusterProfiler::enrichKEGG(
    gene = enrich_input$target_ids,
    universe = enrich_input$universe_target_ids,
    organism = organism,
    keyType = "ncbi-geneid",
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    use_internal_data = use_internal_data
  )

  if (isTRUE(readable) && .has_enrichment_rows(ora_result)) {
    ora_result <- clusterProfiler::setReadable(
      ora_result,
      OrgDb = .get_orgdb_object(organism_db),
      keyType = "ENTREZID"
    )
  }

  list(
    database = "KEGG",
    gene_symbols = enrich_input$gene_symbols,
    entrez_ids = enrich_input$target_ids,
    mapping = enrich_input$mapping,
    mapping_summary = enrich_input$mapping_summary,
    universe_symbols = enrich_input$universe_symbols,
    universe_entrez_ids = enrich_input$universe_target_ids,
    ora_result = ora_result
  )
}


#' Retrieve a STRING PPI network for proteomics hits
#'
#' @description
#' Convert protein identifiers to gene symbols and retrieve a protein-protein
#' interaction network from STRING via `clusterProfiler::getPPI()`.
#'
#' @param proteins A character vector of protein identifiers, or a data.frame
#'   containing a protein identifier column.
#' @param protein_col Optional column name when `proteins` is a data.frame.
#' @param from_type Character string describing the input identifier type for
#'   Bioconductor-based mapping. Default is `"SYMBOL"`.
#' @param input_format Input naming convention passed to
#'   [protein_to_gene_symbol()].
#' @param mapping_table Optional data.frame containing custom protein-to-symbol
#'   mappings.
#' @param mapping_protein_col Column name in `mapping_table` containing protein
#'   identifiers. Default is `"protein"`.
#' @param mapping_symbol_col Column name in `mapping_table` containing gene
#'   symbols. Default is `"gene_symbol"`.
#' @param organism_db Character string naming the OrgDb package. Default is
#'   `"org.Hs.eg.db"`.
#' @param taxID NCBI taxon identifier passed to `clusterProfiler::getPPI()`.
#'   Default is `9606`.
#' @param required_score Optional STRING score cutoff passed to
#'   `clusterProfiler::getPPI()`.
#' @param network_type STRING network type. One of `"functional"` or
#'   `"physical"`. Default is `"functional"`.
#' @param add_nodes Number of partner nodes to add in STRING. Default is `0`.
#' @param show_query_node_labels Passed to `clusterProfiler::getPPI()`. Default
#'   is `0`.
#' @param output One of `"igraph"` or `"data.frame"`. Default is `"igraph"`.
#'
#' @return A list with components `source`, `gene_symbols`, `mapping`, and
#'   `ppi`.
#'
#' @export
get_protein_ppi <- function(proteins,
                            protein_col = NULL,
                            from_type = "SYMBOL",
                            mapping_table = NULL,
                            mapping_protein_col = "protein",
                            mapping_symbol_col = "gene_symbol",
                            organism_db = "org.Hs.eg.db",
                            taxID = 9606,
                            required_score = NULL,
                            network_type = "functional",
                            add_nodes = 0,
                            show_query_node_labels = 0,
                            output = c("igraph", "data.frame"),
                            input_format = c("auto", "rap_olink", "symbol", "raw")) {
  .require_package("clusterProfiler")
  output <- match.arg(output)
  input_format <- match.arg(input_format)

  mapping <- protein_to_gene_symbol(
    proteins = proteins,
    protein_col = protein_col,
    from_type = from_type,
    input_format = input_format,
    mapping_table = mapping_table,
    mapping_protein_col = mapping_protein_col,
    mapping_symbol_col = mapping_symbol_col,
    organism_db = organism_db,
    drop_unmapped = TRUE
  )

  gene_symbols <- unique(mapping$gene_symbol)
  if (!length(gene_symbols)) {
    stop("No gene symbols remained after protein-to-symbol conversion.")
  }

  ppi <- tryCatch(
    clusterProfiler::getPPI(
      x = gene_symbols,
      taxID = taxID,
      required_score = required_score,
      network_type = network_type,
      add_nodes = add_nodes,
      show_query_node_labels = show_query_node_labels,
      output = output
    ),
    error = function(e) {
      stop(sprintf("Failed to retrieve STRING PPI network: %s", e$message))
    }
  )

  list(
    source = "STRING",
    gene_symbols = gene_symbols,
    mapping = mapping,
    ppi = ppi
  )
}


#' Filter a STRING PPI network via TCMDATA
#'
#' @description
#' A thin wrapper around `TCMDATA::ppi_subset()` for STRING-derived PPI
#' networks.
#'
#' @param ppi An `igraph` object, a list returned by `get_protein_ppi()`, or a
#'   list containing a `graph` element.
#' @param n Integer. Number of top-degree nodes to keep. If `NULL`, no degree
#'   filtering is applied.
#' @param score_cutoff Numeric. Minimum STRING confidence score to retain.
#'   Default is `0.7`.
#' @param edge_attr Character. Edge attribute containing the confidence score.
#'   Default is `"score"`.
#' @param rm_isolates Logical. Remove isolated nodes after filtering? Default is
#'   `TRUE`.
#'
#' @return An `igraph` object.
#'
#' @export
subset_protein_ppi <- function(ppi,
                               n = NULL,
                               score_cutoff = 0.7,
                               edge_attr = "score",
                               rm_isolates = TRUE) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "ppi_subset")(
    ppi_obj = graph,
    n = n,
    score_cutoff = score_cutoff,
    edge_attr = edge_attr,
    rm_isolates = rm_isolates
  )
}


#' Compute topological metrics for a PPI network
#'
#' @description
#' A thin wrapper around `TCMDATA::compute_nodeinfo()`.
#'
#' @inheritParams subset_protein_ppi
#' @param weight_attr Character. Edge attribute used as weight. Default is
#'   `"score"`.
#' @param normalize Logical. Whether to normalize betweenness and closeness.
#'   Default is `FALSE`.
#' @param seed Numeric random seed used by the EPC calculation. Default is
#'   `42`.
#'
#' @return An `igraph` object with additional vertex attributes.
#'
#' @export
compute_protein_ppi_metrics <- function(ppi,
                                        weight_attr = "score",
                                        normalize = FALSE,
                                        seed = 42) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "compute_nodeinfo")(
    g = graph,
    weight_attr = weight_attr,
    normalize = normalize,
    seed = seed
  )
}


#' Rank nodes in a PPI network by integrated centrality
#'
#' @description
#' A thin wrapper around `TCMDATA::rank_ppi_nodes()`.
#'
#' @inheritParams subset_protein_ppi
#' @param metrics Character vector of node metrics used for integrated ranking.
#' @param weights Optional numeric weights for `metrics`.
#' @param use_weight Logical. Whether to prefer weighted betweenness and
#'   closeness metrics. Default is `TRUE`.
#' @param na_rm Logical. Whether to ignore missing values during normalization.
#'   Default is `TRUE`.
#'
#' @return A list with components `graph` and `table`.
#'
#' @export
rank_protein_ppi_nodes <- function(ppi,
                                   metrics = c(
                                     "degree",
                                     "betweenness",
                                     "closeness",
                                     "eccentricity",
                                     "radiality",
                                     "Stress",
                                     "MCC",
                                     "MNC",
                                     "DMNC",
                                     "BN",
                                     "EPC"
                                   ),
                                   weights = NULL,
                                   use_weight = TRUE,
                                   na_rm = TRUE) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "rank_ppi_nodes")(
    g = graph,
    metrics = metrics,
    weights = weights,
    use_weight = use_weight,
    na_rm = na_rm
  )
}


#' Cluster a protein-protein interaction network
#'
#' @description
#' Unified interface for community detection in STRING-derived PPI networks.
#' New analyses should call this function and choose the algorithm with
#' \code{method}. Method-specific helper functions are retained internally.
#'
#' @inheritParams subset_protein_ppi
#' @param method Clustering algorithm. Options are \code{"fastgreedy"},
#'   \code{"louvain"}, \code{"mcode"}, and \code{"mcl"}.
#' @param ... Method-specific arguments passed to the selected clustering
#'   helper, such as \code{n_clusters} for \code{"fastgreedy"},
#'   \code{resolution} for \code{"louvain"}, \code{vwp} for \code{"mcode"},
#'   or \code{inflation} for \code{"mcl"}.
#'
#' @return An \code{igraph} object with method-specific cluster attributes.
#' @export
#'
run_protein_ppi_clustering <- function(ppi,
                                       method = c("fastgreedy", "louvain", "mcode", "mcl"),
                                       ...) {
  method <- match.arg(method)
  switch(
    method,
    fastgreedy = run_protein_ppi_fastgreedy(ppi, ...),
    louvain = run_protein_ppi_louvain(ppi, ...),
    mcode = run_protein_ppi_mcode(ppi, ...),
    mcl = run_protein_ppi_mcl(ppi, ...)
  )
}


#' Run Louvain clustering on a PPI network
#'
#' @description
#' A thin wrapper around `TCMDATA::run_louvain()`.
#'
#' @inheritParams subset_protein_ppi
#' @param resolution Numeric resolution parameter for Louvain clustering.
#'   Default is `1.0`.
#' @param weights Optional edge weights.
#'
#' @return An `igraph` object with a `louvain_cluster` vertex attribute.
#'
#' @noRd
run_protein_ppi_louvain <- function(ppi,
                                    resolution = 1.0,
                                    weights = NULL) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "run_louvain")(
    g = graph,
    resolution = resolution,
    weights = weights
  )
}


#' Run fast greedy clustering on a PPI network
#'
#' @description
#' Cluster a PPI network with igraph's fast greedy community detection. By
#' default the function uses the largest connected component and cuts the
#' resulting dendrogram into a user-specified number of clusters. This is useful
#' when sparse STRING networks contain many small disconnected components but
#' the case study needs interpretable modules in the main connected network.
#'
#' @inheritParams subset_protein_ppi
#' @param n_clusters Integer number of clusters to return. Default is `4`.
#' @param largest_component Logical. If `TRUE`, run clustering on the largest
#'   connected component. If `FALSE`, run on the full graph. Default is `TRUE`.
#' @param cluster_attr Character vertex attribute used to store cluster labels.
#'   Default is `"fast_greedy_cluster"`.
#' @param prefix Character prefix for cluster labels. Default is `"FG"`.
#'
#' @return An `igraph` object with a fast greedy cluster vertex attribute.
#'
#' @noRd
run_protein_ppi_fastgreedy <- function(ppi,
                                       n_clusters = 4L,
                                       largest_component = TRUE,
                                       cluster_attr = "fast_greedy_cluster",
                                       prefix = "FG") {
  .require_package("igraph")
  graph <- .extract_ppi_graph(ppi)

  if (!is.numeric(n_clusters) || length(n_clusters) != 1L || n_clusters < 1) {
    stop("`n_clusters` must be a positive integer.")
  }
  n_clusters <- as.integer(n_clusters)

  graph <- simplify(
    as_undirected(graph),
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = "mean"
  )

  if (isTRUE(largest_component)) {
    comp <- components(graph)
    largest <- which.max(comp$csize)
    keep_nodes <- names(comp$membership)[comp$membership == largest]
    graph <- induced_subgraph(graph, vids = keep_nodes)
  }

  if (vcount(graph) == 0L) {
    stop("The graph has no vertices after preprocessing.")
  }
  if (ecount(graph) == 0L) {
    stop("Fast greedy clustering requires at least one edge.")
  }

  n_clusters <- min(n_clusters, vcount(graph))
  communities <- cluster_fast_greedy(graph)
  membership <- cut_at(communities, no = n_clusters)
  labels <- paste0(prefix, membership)
  graph <- set_vertex_attr(graph, cluster_attr, value = labels)
  graph <- set_vertex_attr(graph, "fast_greedy_membership", value = membership)
  graph <- set_vertex_attr(graph, "fast_greedy_n_clusters", value = n_clusters)
  graph
}


#' Run MCODE clustering on a PPI network
#'
#' @description
#' A thin wrapper around `TCMDATA::run_mcode()`.
#'
#' @inheritParams subset_protein_ppi
#' @param vwp Numeric node score cutoff. Default is `0.2`.
#' @param degree_cutoff Numeric degree cutoff used in node weighting. Default
#'   is `2`.
#' @param k_core_threshold Numeric minimum k-core threshold for retained
#'   modules. Default is `2`.
#' @param haircut Logical. Whether to prune singly connected nodes. Default is
#'   `TRUE`.
#' @param fluff Logical. Whether to expand modules with dense neighbors. Default
#'   is `FALSE`.
#' @param fdt Numeric fluff density cutoff. Default is `0.1`.
#' @param loops Logical. Whether to include self-loops. Default is `FALSE`.
#' @param max_depth Numeric maximum recursion depth. Default is `100`.
#'
#' @return An `igraph` object with MCODE-related vertex attributes.
#'
#' @noRd
run_protein_ppi_mcode <- function(ppi,
                                  vwp = 0.2,
                                  degree_cutoff = 2,
                                  k_core_threshold = 2,
                                  haircut = TRUE,
                                  fluff = FALSE,
                                  fdt = 0.1,
                                  loops = FALSE,
                                  max_depth = 100) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "run_mcode")(
    g = graph,
    vwp = vwp,
    degree_cutoff = degree_cutoff,
    k_core_threshold = k_core_threshold,
    haircut = haircut,
    fluff = fluff,
    fdt = fdt,
    loops = loops,
    max_depth = max_depth
  )
}


#' Extract MCODE clustering results from a PPI network
#'
#' @description
#' A thin wrapper around `TCMDATA::get_mcode_res()`.
#'
#' @inheritParams subset_protein_ppi
#' @param only_clusters Logical. If `TRUE`, return only nodes assigned to an
#'   MCODE cluster. Default is `FALSE`.
#'
#' @return A data.frame containing node-level MCODE results.
#'
#' @noRd
get_protein_mcode_res <- function(ppi, only_clusters = FALSE) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "get_mcode_res")(g = graph, only_clusters = only_clusters)
}


#' Run Markov clustering on a PPI network
#'
#' @description
#' A thin wrapper around `TCMDATA::run_MCL()`.
#'
#' @inheritParams subset_protein_ppi
#' @param inflation Numeric inflation parameter controlling cluster granularity.
#'   Default is `2.5`.
#' @param max_iter Integer maximum number of MCL iterations. Default is `100`.
#' @param pruning Numeric pruning threshold. Default is `1e-5`.
#' @param allow1 Logical. Whether singleton clusters are allowed. Default is
#'   `FALSE`.
#'
#' @return An `igraph` object with an `MCL_cluster` vertex attribute.
#'
#' @noRd
run_protein_ppi_mcl <- function(ppi,
                                inflation = 2.5,
                                max_iter = 100,
                                pruning = 1e-5,
                                allow1 = FALSE) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "run_MCL")(
    g = graph,
    inflation = inflation,
    max_iter = max_iter,
    pruning = pruning,
    allow1 = allow1
  )
}


#' Score network clusters in a PPI graph
#'
#' @description
#' A thin wrapper around `TCMDATA::add_cluster_score()`.
#'
#' @inheritParams subset_protein_ppi
#' @param cluster_attr Character. Vertex attribute containing cluster labels.
#'   Default is `"louvain_cluster"`.
#' @param min_size Integer. Minimum cluster size to retain. Default is `3`.
#'
#' @return A data.frame containing cluster scores.
#'
#' @export
score_protein_ppi_clusters <- function(ppi,
                                       cluster_attr = "louvain_cluster",
                                       min_size = 3) {
  .require_package("TCMDATA")
  graph <- .extract_ppi_graph(ppi)
  .get_exported_fun("TCMDATA", "add_cluster_score")(
    g = graph,
    cluster_attr = cluster_attr,
    min_size = min_size
  )
}


#' Evaluate PPI network robustness for selected protein targets
#'
#' @description
#' Convert target protein identifiers to gene symbols and run STRING-network
#' robustness analysis via `TCMDATA::ppi_knock()`.
#'
#' @inheritParams subset_protein_ppi
#' @param targets A character vector of target protein identifiers, or a
#'   data.frame containing a target identifier column.
#' @param target_col Optional column name when `targets` is a data.frame.
#' @param from_type Character string describing the input identifier type for
#'   Bioconductor-based mapping. Default is `"SYMBOL"`.
#' @param input_format Input naming convention passed to
#'   [protein_to_gene_symbol()].
#' @param mapping_table Optional data.frame containing custom protein-to-symbol
#'   mappings.
#' @param mapping_protein_col Column name in `mapping_table` containing protein
#'   identifiers. Default is `"protein"`.
#' @param mapping_symbol_col Column name in `mapping_table` containing gene
#'   symbols. Default is `"gene_symbol"`.
#' @param organism_db Character string naming the OrgDb package. Default is
#'   `"org.Hs.eg.db"`.
#' @param n_perm Integer. Number of permutation iterations. Default is `100`.
#' @param weight_attr Character. Edge attribute containing the confidence score.
#'   Default is `"score"`.
#' @param rewire_niter Integer. Rewiring multiplier used in the null model.
#'   Default is `10`.
#' @param seed Integer random seed. Default is `42`.
#'
#' @return A list with components `targets`, `mapping`, and `robustness`.
#'
#' @export
run_protein_ppi_robustness <- function(ppi,
                                       targets,
                                       target_col = NULL,
                                       from_type = "SYMBOL",
                                       mapping_table = NULL,
                                       mapping_protein_col = "protein",
                                       mapping_symbol_col = "gene_symbol",
                                       organism_db = "org.Hs.eg.db",
                                       n_perm = 100L,
                                       weight_attr = "score",
                                       rewire_niter = 10L,
                                       seed = 42L,
                                       input_format = c("auto", "rap_olink", "symbol", "raw")) {
  .require_package("TCMDATA")
  input_format <- match.arg(input_format)
  graph <- .extract_ppi_graph(ppi)

  mapping <- protein_to_gene_symbol(
    proteins = targets,
    protein_col = target_col,
    from_type = from_type,
    input_format = input_format,
    mapping_table = mapping_table,
    mapping_protein_col = mapping_protein_col,
    mapping_symbol_col = mapping_symbol_col,
    organism_db = organism_db,
    drop_unmapped = TRUE
  )

  target_symbols <- unique(mapping$gene_symbol)
  if (!length(target_symbols)) {
    stop("No target gene symbols remained after protein-to-symbol conversion.")
  }

  robustness <- .get_exported_fun("TCMDATA", "ppi_knock")(
    g = graph,
    targets = target_symbols,
    n_perm = n_perm,
    weight_attr = weight_attr,
    rewire_niter = rewire_niter,
    seed = seed
  )

  list(
    targets = target_symbols,
    mapping = mapping,
    robustness = robustness
  )
}


#' Plot GO ORA results as a bar chart via TCMDATA
#'
#' @description
#' A thin wrapper around `TCMDATA::go_barplot()` for GO enrichment results.
#' This function accepts either a raw `enrichResult` object or a list returned
#' by `run_protein_ora()`.
#'
#' @param x An `enrichResult` object, or a list containing `ora_result`.
#' @param ... Additional arguments passed to `TCMDATA::go_barplot()`.
#'
#' @return A `ggplot2` object.
#'
#' @export
plot_go_ora_bar <- function(x, ...) {
  .require_package("TCMDATA")
  enrich_obj <- .extract_enrichment_result(x)
  .get_exported_fun("TCMDATA", "go_barplot")(enrich_obj, ...)
}


#' Plot enrichment results as a lollipop chart via TCMDATA
#'
#' @description
#' A thin wrapper around `TCMDATA::gglollipop()` for enrichment results. This
#' function accepts either a raw `enrichResult` object or a list returned by one
#' of the proteomics ORA helpers in this package.
#'
#' @param x An `enrichResult` object, or a list containing `ora_result`.
#' @param ... Additional arguments passed to `TCMDATA::gglollipop()`.
#'
#' @return A `ggplot2` object.
#'
#' @export
plot_enrichment_lollipop <- function(x, ...) {
  .require_package("TCMDATA")
  enrich_obj <- .extract_enrichment_result(x)
  .get_exported_fun("TCMDATA", "gglollipop")(enrich_obj, ...)
}


.extract_ppi_graph <- function(x) {
  if (inherits(x, "igraph")) {
    return(x)
  }

  if (is.list(x)) {
    if ("ppi" %in% names(x) && inherits(x$ppi, "igraph")) {
      return(x$ppi)
    }

    if ("graph" %in% names(x) && inherits(x$graph, "igraph")) {
      return(x$graph)
    }
  }

  stop("`ppi` must be an igraph object or a list containing `ppi` or `graph`.")
}


.extract_protein_identifiers <- function(proteins, protein_col = NULL) {
  if (is.data.frame(proteins)) {
    if (is.null(protein_col) || !nzchar(protein_col)) {
      stop("When `proteins` is a data.frame, `protein_col` must be provided.")
    }
    if (!protein_col %in% names(proteins)) {
      stop(sprintf("Column `%s` was not found in `proteins`.", protein_col))
    }
    proteins <- proteins[[protein_col]]
  }

  if (!is.character(proteins)) {
    stop("`proteins` must be a character vector or a data.frame containing one.")
  }

  proteins
}


.normalize_identifier_rows <- function(x,
                                       input_format = c("auto", "rap_olink", "symbol", "raw")) {
  input_format <- match.arg(input_format)
  x <- trimws(as.character(x))
  if (input_format == "auto") {
    rap_olink <- grepl(
      "(^|\\.)olink_instance_[0-9]+\\.",
      x,
      ignore.case = TRUE
    )
    x[rap_olink] <- sub(
      "^.*olink_instance_[0-9]+\\.",
      "",
      x[rap_olink],
      ignore.case = TRUE
    )
    x <- sub(";.*$", "", x)
    x <- gsub("[^A-Za-z0-9_.\\-]", "", x)
  } else if (input_format == "rap_olink") {
    x <- sub("^.*\\.", "", x)
    x <- sub(";.*$", "", x)
    x <- gsub("[^A-Za-z0-9_.\\-]", "", x)
  } else if (input_format == "symbol") {
    invalid <- !is.na(x) & nzchar(x) & !grepl("^[A-Za-z0-9_.\\-]+$", x)
    if (any(invalid)) {
      stop(
        "`input_format = \"symbol\"` received invalid identifier(s): ",
        paste(utils::head(x[invalid], 5L), collapse = ", "),
        call. = FALSE
      )
    }
  }
  x[x == ""] <- NA_character_
  toupper(x)
}


.normalize_identifiers <- function(x,
                                   input_format = c("auto", "rap_olink", "symbol", "raw")) {
  input_format <- match.arg(input_format)
  x <- .normalize_identifier_rows(x, input_format = input_format)
  unique(x[!is.na(x)])
}


.normalize_gene_symbols <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x <- toupper(x)
  x
}


.expand_multi_target_symbols <- function(result) {
  if (!nrow(result)) {
    return(result)
  }

  split_symbols <- strsplit(result$gene_symbol, "_", fixed = TRUE)
  row_index <- rep(seq_len(nrow(result)), lengths(split_symbols))
  expanded <- result[row_index, , drop = FALSE]
  expanded$gene_symbol <- unlist(split_symbols, use.names = FALSE)
  expanded$gene_symbol <- .normalize_gene_symbols(expanded$gene_symbol)
  expanded <- expanded[!is.na(expanded$gene_symbol) & nzchar(expanded$gene_symbol), , drop = FALSE]
  rownames(expanded) <- NULL
  expanded
}


.map_remaining_to_symbol <- function(result, from_type, organism_db) {
  remaining_ids <- unique(result$protein[is.na(result$gene_symbol)])
  if (!length(remaining_ids)) {
    return(result)
  }

  if (identical(toupper(from_type), "SYMBOL")) {
    result$gene_symbol[is.na(result$gene_symbol)] <- remaining_ids[
      match(result$protein[is.na(result$gene_symbol)], remaining_ids)
    ]
    result$mapping_source[is.na(result$mapping_source)] <- "input_symbol"
    return(result)
  }

  .require_package("clusterProfiler")
  .require_package("AnnotationDbi")
  orgdb_obj <- .get_orgdb_object(organism_db)

  mapped <- tryCatch(
    clusterProfiler::bitr(
      remaining_ids,
      fromType = from_type,
      toType = "SYMBOL",
      OrgDb = orgdb_obj
    ),
    error = function(e) {
      stop(sprintf("Failed to map identifiers from `%s` to `SYMBOL`: %s", from_type, e$message))
    }
  )

  if (!nrow(mapped)) {
    return(result)
  }

  mapped <- mapped[!duplicated(mapped[[from_type]]), , drop = FALSE]
  matched_idx <- match(result$protein, mapped[[from_type]])
  has_match <- !is.na(matched_idx) & is.na(result$gene_symbol)
  if (any(has_match)) {
    result$gene_symbol[has_match] <- mapped$SYMBOL[matched_idx[has_match]]
    result$mapping_source[has_match] <- paste0("bitr:", from_type)
  }

  result
}


.prepare_protein_enrichment_input <- function(proteins,
                                              protein_col = NULL,
                                              from_type = "SYMBOL",
                                              input_format = c("auto", "rap_olink", "symbol", "raw"),
                                              mapping_table = NULL,
                                              mapping_protein_col = "protein",
                                              mapping_symbol_col = "gene_symbol",
                                              universe = NULL,
                                              organism_db = "org.Hs.eg.db",
                                              target_type = c("SYMBOL", "ENTREZID")) {
  input_format <- match.arg(input_format)
  target_type <- match.arg(target_type)
  mapping <- protein_to_gene_symbol(
    proteins = proteins,
    protein_col = protein_col,
    from_type = from_type,
    input_format = input_format,
    mapping_table = mapping_table,
    mapping_protein_col = mapping_protein_col,
    mapping_symbol_col = mapping_symbol_col,
    organism_db = organism_db,
    drop_unmapped = TRUE
  )

  gene_symbols <- unique(mapping$gene_symbol)
  if (!length(gene_symbols)) {
    stop("No gene symbols remained after protein-to-symbol conversion.")
  }

  universe_symbols <- NULL
  if (!is.null(universe)) {
    universe_mapping <- protein_to_gene_symbol(
      proteins = universe,
      from_type = from_type,
      input_format = input_format,
      mapping_table = mapping_table,
      mapping_protein_col = mapping_protein_col,
      mapping_symbol_col = mapping_symbol_col,
      organism_db = organism_db,
      drop_unmapped = TRUE
    )
    universe_symbols <- unique(universe_mapping$gene_symbol)
    if (!length(universe_symbols)) {
      stop("`universe` was provided, but no background gene symbols were mapped.")
    }
    missing_from_universe <- setdiff(gene_symbols, universe_symbols)
    if (length(missing_from_universe)) {
      stop(
        "All mapped query proteins must be present in `universe`; missing gene symbols: ",
        paste(utils::head(missing_from_universe, 10L), collapse = ", ")
      )
    }
  }

  target_ids <- gene_symbols
  universe_target_ids <- universe_symbols
  if (identical(target_type, "ENTREZID")) {
    target_ids <- .map_symbols_to_target(
      symbols = gene_symbols,
      target_type = target_type,
      organism_db = organism_db
    )

    if (!is.null(universe_symbols)) {
      universe_target_ids <- .map_symbols_to_target(
        symbols = universe_symbols,
        target_type = target_type,
        organism_db = organism_db
      )
      missing_target_ids <- setdiff(target_ids, universe_target_ids)
      if (length(missing_target_ids)) {
        stop("Mapped query identifiers are not a subset of the mapped ORA universe.")
      }
    }
  }

  list(
    mapping = mapping,
    mapping_summary = list(
      query = attr(mapping, "mapping_summary"),
      universe = if (exists("universe_mapping", inherits = FALSE)) {
        attr(universe_mapping, "mapping_summary")
      } else {
        NULL
      }
    ),
    gene_symbols = gene_symbols,
    universe_symbols = universe_symbols,
    target_ids = target_ids,
    universe_target_ids = universe_target_ids
  )
}


.map_symbols_to_target <- function(symbols,
                                   target_type = "ENTREZID",
                                   organism_db = "org.Hs.eg.db") {
  .require_package("clusterProfiler")
  .require_package("AnnotationDbi")
  orgdb_obj <- .get_orgdb_object(organism_db)
  symbols <- unique(.normalize_gene_symbols(symbols))

  mapped <- tryCatch(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = target_type,
      OrgDb = orgdb_obj
    ),
    error = function(e) {
      stop(sprintf("Failed to map gene symbols to `%s`: %s", target_type, e$message))
    }
  )

  if (!nrow(mapped)) {
    stop(sprintf("No `%s` identifiers were mapped from the supplied gene symbols.", target_type))
  }

  unique(as.character(mapped[[target_type]]))
}


.extract_enrichment_result <- function(x) {
  if (inherits(x, "enrichResult")) {
    return(x)
  }

  if (is.list(x) && "ora_result" %in% names(x)) {
    return(x$ora_result)
  }

  stop("`x` must be an enrichResult object or a list returned by the proteomics ORA helpers.")
}


.has_enrichment_rows <- function(x) {
  inherits(x, "enrichResult") && !is.null(x@result) && nrow(x@result) > 0
}


.get_orgdb_object <- function(organism_db) {
  .require_package(organism_db)
  get(organism_db, envir = asNamespace(organism_db))
}


.require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package `%s` is required for this function. Please install it first.",
        pkg
      ),
      call. = FALSE
    )
  }
}


.get_exported_fun <- function(pkg, fun) {
  .require_package(pkg)
  getExportedValue(pkg, fun)
}
