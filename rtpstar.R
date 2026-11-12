# rtpstar.R -- R port of RTP-STAR (Regression Tree Pipeline for Spatial, Temporal, And Replicate data),
# the RTP-STAR tab of TuxNet (Natalie M. Clark, 2019; RTPSTAR_MAIN.m, regression_tree_pipeline.m,
# run_regressiontree.m, clustering.m, biograph_to_text.m, CheckInputsRTPSTAR.m).
#
# Algorithm (same steps and defaults as the MATLAB version):
#   1. keep the DE genes; optionally cluster them (spatial: k-means, temporal: hierarchical;
#      number of clusters by silhouette) or use clusters you provide;
#   2. per cluster, infer a network with GENIE3 (random forests, regulators = TFs only, 10000 trees),
#      keep the strongest edges (a number that depends on the cluster size and its fraction of TFs),
#      and, if a time course is given, label each edge "activates"/"inhibits"/"regulates";
#   3. optionally connect the clusters through their hubs (nodes with most outgoing edges);
#   4. optionally repeat everything `numiters` times (clustering is random) and keep the edges that
#      appear in at least `maxprop` of the runs;
#   5. write "A activates B" lines for Cytoscape.
#
# Differences from the MATLAB code (deliberate):
#   * genes are matched across files by name (MATLAB matched the clustering/time-course rows to the
#     expression rows by position, which silently mislabels genes if the files are ordered differently);
#   * GENIE3 is re-implemented with `ranger` (same method as the published GENIE3: unit-variance targets,
#     TF inputs, total variance-reduction importance) instead of the GENIE3 MATLAB MEX code;
#   * k-means uses several random starts, constant rows are z-scored to 0 instead of NaN.

# ---- GENIE3 (random forest, K = sqrt) -------------------------------------------------------------
# expr: genes x samples matrix. Returns an n x n matrix W[regulator, target] (zero for non-TF rows).
genie3_weights <- function(expr, tf_idx, ntrees = 10000, cores = 1, seed = NULL) {
  n <- nrow(expr); X <- t(expr)                       # samples x genes
  W <- matrix(0, n, n, dimnames = list(rownames(expr), rownames(expr)))
  for (j in seq_len(n)) {
    inputs <- setdiff(tf_idx, j)
    if (!length(inputs)) next
    y <- X[, j]; s <- sd(y)
    if (!is.finite(s) || s == 0) next                  # constant target: no information
    y <- y / s                                         # unit-variance targets, as in GENIE3
    xin <- X[, inputs, drop = FALSE]
    if (length(inputs) == 1L) {                        # one regulator: its total variance reduction on a stump forest
      fit <- ranger::ranger(x = data.frame(v = xin[, 1]), y = y, num.trees = ntrees, mtry = 1,
                            importance = "impurity", num.threads = cores, seed = seed)
    } else {
      fit <- ranger::ranger(x = as.data.frame(xin, optional = TRUE), y = y, num.trees = ntrees,
                            mtry = max(1, round(sqrt(length(inputs)))), importance = "impurity",
                            num.threads = cores, seed = seed)
    }
    W[inputs, j] <- fit$variable.importance
  }
  W
}

# ---- gene labels (symbols, made unique) -----------------------------------------------------------
rtp_labels <- function(ids, symbol_df) {
  lab <- ids
  if (!is.null(symbol_df) && nrow(symbol_df) && ncol(symbol_df) >= 2) {
    sy <- map_symbols(ids, symbol_df); used <- character()
    for (i in seq_along(ids)) if (sy[i] != ids[i]) {
      if (sy[i] %in% used) next                        # symbol already taken by an earlier gene: keep the id
      lab[i] <- sy[i]; used <- c(used, sy[i])
    }
  }
  dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
  lab[dup] <- ids[dup]
  lab
}

# ---- clustering (clustering.m) ----------------------------------------------------------------------
# data: genes x conditions numeric matrix. Returns integer cluster labels (named by gene).
rtp_cluster <- function(data, type = c("S", "T"), seed = NULL, maxclusters = NULL) {
  type <- match.arg(type)
  if (is.null(seed)) set.seed(NULL) else set.seed(seed)
  n <- nrow(data)
  if (is.null(maxclusters)) maxclusters <- floor(n / 10 + 5)
  minclusters <- if (floor(n / 10) <= 5) floor(n / 10) else floor(n / 10 - 5)
  if (minclusters > maxclusters) minclusters <- max(1, maxclusters - 10)
  Z <- (data - rowMeans(data)) / apply(data, 1, sd)
  Z[!is.finite(Z)] <- 0                                # constant rows
  ks <- seq(max(2, minclusters), maxclusters); ks <- ks[ks < n]
  if (!length(ks) || n < 3) return(setNames(rep(1L, n), rownames(data)))
  d <- dist(Z)
  if (type == "S") {
    sil <- vapply(ks, function(k) mean(cluster::silhouette(kmeans(Z, k, nstart = 5)$cluster, d)[, "sil_width"]), numeric(1))
    k <- ks[which.max(sil)]
    cl <- kmeans(Z, k, nstart = 5)$cluster
  } else {
    hc_eval <- hclust(d, "ward.D2")
    sil <- vapply(ks, function(k) mean(cluster::silhouette(cutree(hc_eval, k), d)[, "sil_width"]), numeric(1))
    k <- ks[which.max(sil)]
    cl <- cutree(hclust(d, "single"), k)               # (sic) MATLAB's linkage() default is single linkage
  }
  setNames(as.integer(cl), rownames(data))
}

# ---- one network (run_regressiontree.m) -------------------------------------------------------------
# expr: genes x samples; is_tf: logical; tc: genes x time points (or NULL). labels: display names.
# Returns list(edges = data.frame(from, to, dir), hub = character()) or NULL when nothing can be inferred.
rtp_infer <- function(expr, is_tf, tc, labels, time_threshold = 1.25, edge_number = c(0.5, 1.5, 2),
                      ntrees = 10000, cores = 1, seed = NULL) {
  n <- nrow(expr)
  if (n <= 1 || sum(is_tf) == 0) return(NULL)
  W <- genie3_weights(expr, which(is_tf), ntrees, cores, seed)
  nz <- which(W != 0, arr.ind = TRUE)
  if (!nrow(nz)) return(list(edges = data.frame(from = character(), to = character(), dir = character()),
                             hub = character()))
  w_sorted <- sort(W[nz], decreasing = TRUE)
  frac <- sum(is_tf) / n
  ntk <- if (n <= 2) 1 else if (n == 3) 5 else
    if (frac < 0.25) floor(edge_number[1] * n) else if (frac < 0.5) floor(edge_number[2] * n) else floor(edge_number[3] * n)
  thr <- if (ntk > length(w_sorted)) 0 else w_sorted[max(1, ntk)]
  keep <- which(W != 0 & W >= thr, arr.ind = TRUE)
  from <- keep[, 1]; to <- keep[, 2]

  dir <- rep("regulates", length(from))
  if (!is.null(tc) && ncol(tc) >= 3) {
    Tn <- ncol(tc); kk <- Tn - 2
    for (e in seq_along(from)) {
      a <- tc[from[e], 2:(Tn - 1)] / tc[from[e], 1:(Tn - 2)]       # regulator change at step k
      b <- tc[to[e], 3:Tn] / tc[to[e], 2:(Tn - 1)]                  # target change at step k+1
      ga <- .isT(a > time_threshold); la <- .isT(a < time_threshold)
      gb <- .isT(b > time_threshold); lb <- .isT(b < time_threshold)
      rep_n <- sum((ga & lb) | (la & gb)); act_n <- sum((ga & gb) | (la & lb))
      dir[e] <- if (rep_n > floor(kk / 2)) "inhibits" else if (act_n > floor(kk / 2)) "activates" else "regulates"
    }
  }
  cnt <- table(from)
  hub <- labels[as.integer(names(cnt)[cnt == max(cnt)])]
  list(edges = data.frame(from = labels[from], to = labels[to], dir = dir, stringsAsFactors = FALSE), hub = hub)
}
# ---- one full pass: clustering -> per-cluster networks -> hubs (regression_tree_pipeline.m) ------------
# expr_df: data.frame(gene, TF, expression...) already filtered to the DE genes.
rtp_pass <- function(ids, is_tf, expr, tc, clu, labels, connect_hubs, seed, clustering_type, preset, cluster_file,
                     time_threshold, ntrees, cores, rf_seed, progress_base = 0, progress_span = 100, symbol_df = NULL) {
  fmt <- function(r) if (is.null(r) || !nrow(r$edges)) character() else paste(r$edges$from, r$edges$dir, r$edges$to)

  cl <- NULL
  if (!is.null(preset)) {
    cl <- setNames(as.integer(preset[[3]]), preset[[1]])
  } else if (!is.null(clu)) {
    cl <- rtp_cluster(clu, clustering_type, seed)
    tab <- data.frame(AGI = names(cl), gene = rtp_labels(names(cl), symbol_df), cluster = unname(cl), stringsAsFactors = FALSE)
    utils::write.csv(tab, cluster_file, row.names = FALSE)
  }

  if (is.null(cl)) {
    emit_progress(progress_base + progress_span * 0.1, "Inferring network")
    r <- rtp_infer(expr, is_tf, tc, labels, time_threshold, ntrees = ntrees, cores = cores, seed = rf_seed)
    return(c_unique_sort(fmt(r)))
  }
  K <- max(cl); lines <- character(); hubs <- character()
  for (i in seq_len(K)) {
    emit_progress(progress_base + progress_span * 0.9 * (i - 1) / K, sprintf("Cluster %d of %d", i, K))
    sel <- which(ids %in% names(cl)[cl == i])
    r <- rtp_infer(expr[sel, , drop = FALSE], is_tf[sel], if (is.null(tc)) NULL else tc[sel, , drop = FALSE],
                   labels[sel], time_threshold, ntrees = ntrees, cores = cores, seed = rf_seed)
    if (!is.null(r)) { lines <- c(lines, fmt(r)); hubs <- c(hubs, r$hub) }
  }
  if (connect_hubs && length(hubs) > 1) {
    emit_progress(progress_base + progress_span * 0.92, "Connecting cluster hubs")
    sel <- which(labels %in% hubs)
    if (length(sel) > 1) {
      r <- rtp_infer(expr[sel, , drop = FALSE], is_tf[sel], if (is.null(tc)) NULL else tc[sel, , drop = FALSE],
                     labels[sel], time_threshold, ntrees = ntrees, cores = cores, seed = rf_seed)
      lines <- c(lines, fmt(r))
    }
  }
  c_unique_sort(lines)
}

# ---- input checks (CheckInputsRTPSTAR) ------------------------------------------------------------------
rtpstar_check_inputs <- function(genes_file, expression_file, clustering_file = "", time_course_file = "",
                                 symbol_file = "", tf_file = "", output_folder = "", output_file = "biograph.txt",
                                 preset_clusters_file = "") {
  p <- character()
  req <- function(path, what) if (!nzchar(path)) p <<- c(p, paste("Please provide", what, ".")) else
    if (!file.exists(path)) p <<- c(p, paste("Provided", what, "cannot be found."))
  req(genes_file, "a genes file"); req(expression_file, "a gene expression file")
  for (x in list(list(clustering_file, "clustering file"), list(time_course_file, "time course file"),
                 list(symbol_file, "gene name file"), list(tf_file, "TF file"),
                 list(preset_clusters_file, "clusters file")))
    if (nzchar(x[[1]]) && !file.exists(x[[1]])) p <- c(p, paste("Provided", x[[2]], "cannot be found."))
  if (!nzchar(output_folder)) p <- c(p, "Please provide an outputs folder.")
  if (!nzchar(output_file)) p <- c(p, "Please provide a results file name.")
  if (nzchar(output_folder) && nzchar(output_file)) {                      # the results must never replace one of the researcher's own files
    ins <- c(`genes file` = genes_file, `expression file` = expression_file, `clustering data file` = clustering_file, `time course file` = time_course_file,
             `gene name file` = symbol_file, `TF file` = tf_file, `clusters file` = preset_clusters_file)
    own <- c(output_file, "expressionWithTFs.csv", "clusters.csv", sprintf("clusters%d.csv", 1:20), sprintf("biograph%d.txt", 1:20))
    p <- c(p, output_overwrites_input(ins[!is.na(ins) & nzchar(ins)], file.path(output_folder, own)))
  }
  p
}

# ---- main (RTPSTAR_MAIN.m + CheckInputsRTPSTAR.m) -------------------------------------------------------
#' @param genes_file       DE genes, first column (Excel/CSV/TXT, first row = header)
#' @param expression_file  gene, [TF indicator], expression columns ...
#' @param tf_file          optional list of TFs (first column); used if the expression file has no "TF" column.
#'                         With neither, every gene is treated as a potential regulator.
#' @param numiters         number of runs to combine (default 1); maxprop: edge proportion kept (default 0.33)
#' @param clustering_type  "S" spatial (k-means) or "T" temporal (hierarchical)
#' @param ntrees           trees per GENIE3 forest (MATLAB: 10000)
rtpstar_run <- function(genes_file, expression_file, output_folder, output_file = "biograph.txt",
                        clustering_file = NULL, time_course_file = NULL, symbol_file = NULL, tf_file = NULL,
                        preset_clusters_file = NULL, numiters = 1, maxprop = 0.33, connect_hubs = TRUE,
                        clustering_seed = NULL, clustering_type = "S", ntrees = 10000, cores = 1, rf_seed = NULL) {
  blank <- function(x) is.null(x) || is.na(x) || !nzchar(x)
  nz <- function(x) if (blank(x)) "" else x
  problems <- rtpstar_check_inputs(genes_file, expression_file, nz(clustering_file), nz(time_course_file),
                                   nz(symbol_file), nz(tf_file), output_folder, output_file, nz(preset_clusters_file))
  if (length(problems)) stop(paste(problems, collapse = " "), call. = FALSE)
  if (numiters < 1 || numiters != round(numiters)) stop("numiters must be an integer >= 1", call. = FALSE)
  if (!clustering_type %in% c("S", "T")) stop("clustering_type must be \"S\" or \"T\"", call. = FALSE)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  out <- function(f) file.path(output_folder, f)
  if (!blank(clustering_seed) && clustering_seed == 0) clustering_seed <- NULL

  emit_progress(2, "Reading input files")
  de <- unique(read_id_column(genes_file))
  ex <- read_table_any(expression_file)
  ex <- ex[ex[[1]] %in% de & !is.na(ex[[1]]) & !duplicated(ex[[1]]), , drop = FALSE]
  if (nrow(ex) < 2) stop("Fewer than 2 genes of the genes file were found in the expression file.", call. = FALSE)

  # TF indicator column (CheckInputsRTPSTAR): use a "TF" column if present, else build one
  if ("TF" %in% names(ex)) {
    tfcol <- ex[["TF"]]
    ex <- ex[c(1, setdiff(seq_along(ex), c(1, which(names(ex) == "TF"))))]
    tv <- tolower(trimws(as.character(tfcol))); tn <- suppressWarnings(as.numeric(tv))                    # 1/0, TRUE/FALSE, yes/no; an empty cell is "not a TF"
    is_tf <- (tv %in% c("true", "t", "yes", "y", "tf")) | (!is.na(tn) & tn == 1)
  } else if (!blank(tf_file)) {
    is_tf <- ex[[1]] %in% read_id_column(tf_file)
    tf_df <- cbind(ex[1], TF = is_tf, ex[-1]); utils::write.csv(tf_df, out("expressionWithTFs.csv"), row.names = FALSE)
  } else {
    is_tf <- rep(TRUE, nrow(ex))
    message("Note: no TF file or TF column given; every gene is treated as a potential regulator.")
  }
  ids <- ex[[1]]
  expr <- as.matrix(sapply(ex[-1], function(v) suppressWarnings(as.numeric(v))))
  if (!is.matrix(expr)) expr <- matrix(expr, nrow = length(ids))
  rownames(expr) <- ids
  if (anyNA(expr)) stop("The expression file contains non-numeric values in its data columns.", call. = FALSE)
  threshold <- 2 * sum(is_tf)                                    # minimum number of edges kept over iterations

  align <- function(path) {                                       # rows by gene name, in the expression order
    if (blank(path)) return(NULL)
    m <- note_dropped_columns(table_to_matrix(read_table_any(path)), "file")
    m <- m[match(ids, rownames(m)), , drop = FALSE]
    rownames(m) <- ids
    m
  }
  tc <- align(time_course_file)
  if (!is.null(tc) && ncol(tc) < 3) stop("The time course file needs at least 3 time points.", call. = FALSE)
  clu <- align(clustering_file)
  if (!is.null(clu)) clu <- clu[rowSums(is.na(clu)) == 0, , drop = FALSE]    # only genes present in the clustering file
  sym <- if (blank(symbol_file)) NULL else read_table_any(symbol_file, header = FALSE)   # name files often have no header row
  labels <- rtp_labels(ids, sym)
  preset <- if (blank(preset_clusters_file)) NULL else read_table_any(preset_clusters_file)
  if (!is.null(preset) && ncol(preset) < 3) stop("The clusters file needs the cluster number in column 3.", call. = FALSE)
  if (!is.null(clu) && nrow(clu) < 3) stop("Fewer than 3 genes of the genes file were found in the clustering file.", call. = FALSE)
  use_clusters <- !is.null(clu) || !is.null(preset)
  message(sprintf("%d genes (%d regulators), %s, %s", length(ids), sum(is_tf),
                  if (is.null(tc)) "no time course (edges are undirected in sign)" else sprintf("time course with %d points", ncol(tc)),
                  if (!is.null(preset)) "using your clusters" else if (!is.null(clu)) paste(if (clustering_type == "S") "spatial" else "temporal", "clustering") else "no clustering"))

  one <- function(i, cluster_file) rtp_pass(ids, is_tf, expr, tc, clu, labels, connect_hubs, clustering_seed,
                                           clustering_type, preset, cluster_file, 1.25, ntrees, cores, rf_seed,
                                           progress_base = 5 + 90 * (i - 1) / numiters, progress_span = 90 / numiters,
                                           symbol_df = sym)
  final_path <- out(output_file)
  if (numiters == 1) {
    emit_progress(5, "Running RTP-STAR")
    lines <- one(1, out("clusters.csv"))
    writeLines(lines, final_path)
  } else {
    runs <- vector("list", numiters)
    for (i in seq_len(numiters)) {
      emit_progress(5 + 90 * (i - 1) / numiters, sprintf("Iteration %d of %d", i, numiters))
      runs[[i]] <- one(i, out(sprintf("clusters%d.csv", i)))
      writeLines(runs[[i]], out(sprintf("biograph%d.txt", i)))
    }
    edges <- c_unique_sort(unlist(runs))
    number <- vapply(edges, function(e) sum(vapply(runs, function(r) e %in% r, logical(1))), numeric(1))
    ord <- order(-number)                                         # stable: ties keep alphabetical order
    edges <- edges[ord]; number <- number[ord]
    if (length(edges)) {
      kth <- if (threshold >= 1) number[min(threshold, length(number))] else Inf
      keep <- if (kth < maxprop * numiters) number >= kth else number >= maxprop * numiters
      writeLines(sprintf("%s %s", edges[keep], format(number[keep] / numiters, digits = 15)), final_path)
    } else writeLines(character(), final_path)
  }
  emit_progress(100, "Done")
  invisible(final_path)
}
