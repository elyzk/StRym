# analysis/genist_network.R -- GENIST regulatory network, started automatically from the results of the entropy steps.
#
# GENIST needs two files the researcher normally builds by hand: a GENES file (which genes go into the network) and a TIME COURSE file
# (their expression at every time point, earliest first). Here StRym builds both for every treated-vs-control pair of the experiment and
# runs GENIST on them:
#   genes        the main-set genes the entropy step(s) of this run found for that treated condition, strongest first (at most
#                `genist_max_genes`, 500 by default: the run time grows steeply with the number of genes)
#   time course  the mean expression of the treated condition at each time point, from the same expression data the other steps use
#                (a Tux run, or your own FPKM / counts table); gene ids are replaced by readable names where the data has them
# No transcription-factor list is needed: without one every gene is a possible regulator (give a TF list in the GENIST tab to restrict it).
# Output : <step folder>/<treated condition>/inputs/  genist_genes.csv, genist_time_course.csv   (usable in the GENIST tab too)
#          <step folder>/<treated condition>/results/ <project>_<comparison>_GENIST_network_for_Cytoscape.txt (the network table),
#                                                     a *_network.png picture, clusters.xlsx and GENIST's other files
# Input  : bulk time course (>= 3 time points) with entropy results (Cuffdiff q-values) -- the entropy steps must be run in the same go.

.genist_gene_labels <- function(ids, gene_names) {                              # readable name when known and unique, else the id
  nm <- if (is.null(gene_names)) rep(NA_character_, length(ids)) else unname(gene_names[ids])
  nm <- trimws(sub("[,;].*$", "", nm)); nm[is.na(nm) | !nzchar(nm) | nm %in% c("-", "NA")] <- NA
  lab <- ifelse(is.na(nm), ids, nm)
  dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
  ifelse(dup, ids, lab)                                                         # two loci with one name: keep both, by id
}

strym_step_genist_network <- function(ctx) {
  p <- strym_normalize_params(ctx$design$params)
  if (!exists("genist_run", mode = "function")) source(file.path(ctx$app_dir, "genist.R"), local = FALSE)
  data <- strym_drop_samples(strym_load_expression(ctx), p)
  S <- data$samples
  tps <- strym_order_timepoints(S$time_point, ctx$design$timepoint_names)
  if (length(tps) < 3) stop(sprintf("GENIST needs a time course of at least 3 time points; the sample sheet has %d.", length(tps)), call. = FALSE)
  pairs <- strym_pairs(S, p); keys <- make.unique(gsub("[^A-Za-z0-9._-]+", "_", .pair_keys(pairs)))        # folder names: no slashes, colons ...
  old_opt <- options(strym.quiet_progress = TRUE); on.exit(options(old_opt), add = TRUE)                  # GENIST's own progress markers would overwrite the design page's bar

  # the genes the entropy step(s) of this run found (their "main set"), per comparison
  found <- list(); n_files <- 0L
  for (i in names(ctx$prereq)) for (ff in list.files(ctx$prereq[[i]], pattern = "_main_genes\\.csv$", full.names = TRUE)) {
    x <- tryCatch(utils::read.csv(ff, stringsAsFactors = FALSE, colClasses = c(gene = "character")), error = function(e) NULL)      # ids as text: 00123 stays 00123
    if (!is.null(x) && all(c("gene", "comparison") %in% names(x))) { n_files <- n_files + 1L; if (nrow(x)) found[[length(found) + 1]] <- x }
  }
  if (!n_files) stop("No genes to build the network from: tick Shannon entropy or Modified Shannon entropy (their main-set genes are used).", call. = FALSE)
  if (!length(found)) return(list(status = "skipped", files = character(), message = sprintf("The entropy step(s) found no genes that change at >= %g time points, so there is nothing to build a network from. Lower \u201cchanged at least this many time points\u201d, the q-value or the fold-change limit in Analysis settings.", p$min_sig_tp)))
  found <- do.call(rbind, lapply(found, function(d) { if (!"n_sig_timepoints" %in% names(d)) d$n_sig_timepoints <- 0; d[, c("gene", "comparison", "n_sig_timepoints")] }))

  M <- strym_wide(data, na_bad = TRUE); genes_all <- rownames(M)
  # gene labels in the files: readable names when the data has them, else ids. With a transcription-factor list, whichever of the two
  # matches the list's ids better is used (so "AT1G01030" in your list finds the gene whether the table's key column is the id or a symbol).
  tf_file <- p$genist_tf_file
  if (nzchar(tf_file) && !file.exists(tf_file)) stop("The transcription factor list was not found: ", tf_file, call. = FALSE)
  lab <- .genist_gene_labels(genes_all, data$gene_names); names(lab) <- genes_all
  if (nzchar(tf_file)) {
    tf_ids <- read_id_column(tf_file)
    if (sum(genes_all %in% tf_ids) >= sum(lab %in% tf_ids)) { lab <- genes_all; names(lab) <- genes_all }
    if (!any(lab %in% tf_ids)) stop("None of the genes in the data are in the transcription factor list (do the gene ids match? e.g. ", paste(utils::head(unique(lab), 3), collapse = ", "),
                                    " vs ", paste(utils::head(tf_ids, 3), collapse = ", "), ").", call. = FALSE)
  }
  top <- ctx$out_dir; files <- character(); notes <- character(); nrun <- 0L; failed_any <- FALSE
  for (i in seq_len(nrow(pairs))) {
    tr <- pairs$treated[i]; key <- keys[i]; pk <- paste0(tr, "_vs_", pairs$control[i])
    ctx$sub_progress((i - 1) / nrow(pairs), sprintf("GENIST for %s (%d of %d)", tr, i, nrow(pairs)))
    x <- found[tolower(found$comparison) == tolower(pk) & found$gene %in% genes_all, , drop = FALSE]
    x <- x[order(-x$n_sig_timepoints, x$gene, method = "radix"), , drop = FALSE]; g <- unique(x$gene)
    n_found <- length(g)
    if (p$genist_max_genes > 0 && length(g) > p$genist_max_genes) g <- g[seq_len(p$genist_max_genes)]
    if (length(g) < 2) { notes <- c(notes, sprintf("%s: only %d gene(s) in the entropy main set, not enough for a network.", tr, length(g))); next }
    tps_i <- tps[vapply(tps, function(t) any(S$condition == tr & S$time_point == t), TRUE)]       # the time points this condition has
    if (length(tps_i) < 3) { notes <- c(notes, sprintf("%s: only %d time point(s) with samples, GENIST needs 3.", tr, length(tps_i))); next }
    if (length(tps_i) < length(tps)) notes <- c(notes, sprintf("%s: no samples at %s (left out of its time course)", tr, paste(setdiff(tps, tps_i), collapse = ", ")))
    tc <- do.call(cbind, lapply(tps_i, function(t) {
      cols <- S$sample[S$condition == tr & S$time_point == t]
      v <- suppressWarnings(rowMeans(M[, cols, drop = FALSE], na.rm = TRUE)); v[is.na(v)] <- 0; v }))
    colnames(tc) <- tps_i; rownames(tc) <- genes_all
    tc <- tc[apply(tc, 1, max) > 0 | genes_all %in% g, , drop = FALSE]; tc <- tc[order(lab[rownames(tc)], method = "radix"), , drop = FALSE]
    d_in <- file.path(top, key, "inputs"); d_out <- file.path(top, key, "results"); dir.create(d_in, recursive = TRUE, showWarnings = FALSE)
    gf <- file.path(d_in, "genist_genes.csv"); tf <- file.path(d_in, "genist_time_course.csv")
    utils::write.csv(data.frame(gene = c_sort(unname(lab[g]))), gf, row.names = FALSE)
    utils::write.csv(data.frame(gene = unname(lab[rownames(tc)]), as.data.frame(tc, check.names = FALSE), check.names = FALSE), tf, row.names = FALSE)
    r <- tryCatch(genist_run(gf, tf, tf_file = if (nzchar(tf_file)) tf_file else NULL, output_folder = d_out,
                             output_file = strym_network_filename(ctx$design$project_name, "GENIST", key),       # e.g. Heat_shock_Heat_vs_Control_GENIST_network_for_Cytoscape.txt
                             time_lapse = 1, fc_thresh = 1.1, time_percent = 0.1, n_levels = 3, low_conn = 0.05, cores = max(1L, min(as.integer(ctx$design$threads), 4L))),
                  error = function(e) e)
    if (inherits(r, "error")) { notes <- c(notes, sprintf("%s: GENIST failed (%s)", tr, conditionMessage(r))); failed_any <- TRUE; next }
    nrun <- nrun + 1L; files <- c(files, file.path(key, "inputs", c("genist_genes.csv", "genist_time_course.csv")), file.path(key, "results", list.files(d_out)))
    notes <- c(notes, sprintf("%s: network from %d genes%s over %d time points", tr, length(g), if (n_found > length(g)) sprintf(" (the %d most often called of the %d in the main set; change the limit in Analysis settings)", length(g), n_found) else "", length(tps_i)))
  }
  if (!nrun && !failed_any) return(list(status = "skipped", files = files, message = paste(c("No network was built:", notes), collapse = " ")))
  if (!nrun) stop(paste(c("No network could be built.", notes), collapse = " "), call. = FALSE)
  list(status = "ok", files = files, message = paste(notes, collapse = "; "))
}
