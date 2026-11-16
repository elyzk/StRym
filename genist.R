# genist.R -- R port of GENIST (Gene regulatory Network Inference from Spatio-Temporal data),
# the GENIST tab of TuxNet (M. Angels de Luis Balaguer, 2016; GENIST2_0.m, GRN_preprocessing2_0.m,
# GRN_inference2_0.m, matrix_to_cytoscape_table2_0.m).
#
# Algorithm (same steps and defaults as the MATLAB version):
#   1. keep the genes of interest; optionally cluster them on a spatial/condition dataset
#      (z-scored, Ward hierarchical clustering; number of clusters by silhouette);
#   2. inside each cluster, pick candidate regulators from the time course (genes whose
#      expression changed >= `fc_thresh` in the previous / same time point), discretize the
#      time course, and choose for each target the regulator set (1-3 genes) with the best
#      Bayesian (BDe) score; the sign comes from the direction of the changes;
#   3. connect the clusters through their hubs (genes with most outgoing edges) the same way;
#   4. write a Cytoscape-ready edge table (source, target, width, sign) and the cluster lists.
#
# Differences from the MATLAB code (deliberate):
#   * genes are matched across files by name (MATLAB assumed identical row order after sorting and
#     fails if a gene is missing from one file); genes missing from a file are dropped with a note;
#   * cluster numbering follows R's cutree (MATLAB's `cluster` may number clusters differently);
#   * everything else, including tie-breaking, follows the MATLAB code.

# ---- Step 2a: candidate regulators + discretization (GRN_preprocessing2_0) ---------------------
genist_preprocess <- function(temp, fc, pct, n_levels, tf_pos = NULL) {
  n <- nrow(temp); Tn <- ncol(temp)
  change <- matrix(FALSE, n, Tn); up <- change; dn <- change
  for (i in seq_len(Tn)[-1]) {
    a <- temp[, i]; b <- temp[, i - 1]
    change[, i] <- .isT(a > fc * b) | .isT(a < (1 / fc) * b)
    up[, i] <- .isT(a > 1.1 * b)             # used for the sign of a regulation
    dn[, i] <- .isT(a < 0.8 * b)
  }

  # reg1[r, t]: r changed in the PREVIOUS time point of a change of t;  reg2: in the SAME time point
  reg1 <- matrix(FALSE, n, n); reg2 <- reg1
  if (Tn >= 2) for (t in seq_len(n)) {
    js <- which(change[t, 2:Tn]) + 1L
    if (length(js)) {
      reg1[, t] <- rowSums(change[, js - 1L, drop = FALSE]) > 0
      reg2[, t] <- rowSums(change[, js, drop = FALSE]) > 0
    }
  }
  if (!is.null(tf_pos)) { reg1[!tf_pos, ] <- FALSE; reg2[!tf_pos, ] <- FALSE }
  diag(reg1) <- FALSE; diag(reg2) <- FALSE            # no autoregulation

  # keep only pairs whose change patterns agree in at least `pct` of the time points
  # (MATLAB counts positions where the two 0/1 vectors are EQUAL, zeros included)
  C <- change * 1; thr <- pct * Tn
  if (Tn >= 2) {
    A1 <- C[, 2:Tn, drop = FALSE]; B1 <- C[, 1:(Tn - 1), drop = FALSE]
    S1 <- A1 %*% t(B1) + (1 - A1) %*% t(1 - B1)        # S1[target, regulator]
    reg1 <- reg1 & t(S1 >= thr)
  }
  S2 <- C %*% t(C) + (1 - C) %*% t(1 - C)
  reg2 <- reg2 & t(S2 >= thr)

  # discretize
  if (n_levels == 2) {
    disc <- (temp >= rowMeans(temp, na.rm = TRUE)) * 1
  } else {
    q <- t(apply(temp, 1, quantile, probs = (1:(n_levels - 1)) / n_levels, type = 5, names = FALSE, na.rm = TRUE))      # an empty cell in the time course counts as "no information", not a crash
    disc <- matrix(0, n, Tn)
    for (g in seq_len(n)) {
      x <- temp[g, ]; qq <- q[g, ]
      d <- numeric(Tn)                                # values below the first quantile stay 0
      for (k in 2:(n_levels - 1)) d[.isT(x > qq[k - 1] & x < qq[k])] <- k - 1
      d[.isT(x > qq[n_levels - 1])] <- n_levels - 1
      disc[g, ] <- d
    }
  }
  max_reg <- max(colSums(reg1), colSums(reg2))
  list(disc = disc, reg1 = reg1, reg2 = reg2, set_max_reg = min(3, max_reg),
       up = up, dn = dn, n_levels = n_levels)
}

# ---- Step 2b: Bayesian (BDe) score of a regulator set (inside GRN_inference2_0) ----------------
# Emat: regulators x time (values 0..L-1); t: target states (0..L-1).
# log P(t | regulators) with the same equivalent sample size as the MATLAB code (1e-15).
bde_score <- function(Emat, t, L, ESS = 1e-15) {
  m <- nrow(Emat); q <- L^m
  idx <- colSums(Emat * L^(seq_len(m) - 1))          # configuration index 0..q-1 per time point
  Nj  <- tabulate(idx + 1, q)
  Njk <- tabulate(idx * L + t + 1, q * L)
  a <- ESS / q; b <- ESS / (L * q)
  sum(lgamma(a) - lgamma(Nj + a)) + sum(lgamma(Njk + b) - lgamma(b))
}

# ---- Step 2c: choose regulators for every target (GRN_inference2_0) ----------------------------
# time_lapse: 1 (previous time point), 0 (same time point) or -1 (both).
genist_infer_orig <- function(pre, time_lapse, low_conn, cores = 1) {
  disc <- pre$disc; up <- pre$up; dn <- pre$dn; L <- pre$n_levels
  n <- nrow(disc); Tn <- ncol(disc); set_max_reg <- pre$set_max_reg

  per_layer <- function(R, shift) {
    cE <- seq_len(Tn - shift); ct <- (1L + shift):Tn
    res <- parallel::mclapply(seq_len(n), function(j) {
      ind3 <- which(R[, j]); nreg <- length(ind3)
      if (nreg == 0L) return(NULL)
      tv <- disc[j, ct]; tu <- up[j, ct]; td <- dn[j, ct]
      best <- -Inf; idx <- integer(); sg <- numeric()
      consider <- function(regs) {
        s <- bde_score(disc[regs, cE, drop = FALSE], tv, L)
        best <<- max(best, s)
        if (s == best && s != -Inf) {                 # ties replace the earlier candidate (as in MATLAB)
          idx <<- regs
          sg <<- vapply(regs, function(r) {
            eu <- up[r, cE]; ed <- dn[r, cE]
            sign(sum((eu & tu) | (ed & td)) - sum((eu & td) | (ed & tu)))
          }, numeric(1))
        }
      }
      for (k in nreg:1) {
        consider(ind3[k])
        if (nreg > 1 && set_max_reg > 1) for (kk in seq_len(k - 1)) {
          consider(ind3[c(k, kk)])
          if (nreg > 2 && set_max_reg > 2) for (kkk in seq_len(kk - 1)) consider(ind3[c(k, kk, kkk)])
        }
      }
      list(best = best, idx = idx, sign = sg)
    }, mc.cores = cores)
    net <- matrix(0, n, n); sgn <- matrix(0, n, n)
    for (j in seq_len(n)) {
      r <- res[[j]]
      if (is.null(r) || !length(r$idx)) next
      net[r$idx, j] <- 10^r$best                      # (sic) the MATLAB code undoes the log with base 10
      sgn[r$idx, j] <- r$sign
    }
    list(net = net, sgn = sgn)
  }

  net <- NULL
  if (time_lapse != 0) l1 <- per_layer(pre$reg1, 1L)
  if (time_lapse != 1) l2 <- per_layer(pre$reg2, 0L)
  if (time_lapse == 0) { net <- l2$net; sgn <- l2$sgn
  } else if (time_lapse == -1) {
    net <- l1$net + l2$net; sgn <- l1$sgn + l2$sgn
    sgn[sgn == 2] <- 1; sgn[sgn == -2] <- -1          # opposite signs from the two layers cancel to 0
  } else { net <- l1$net; sgn <- l1$sgn }

  m1 <- max(abs(net)); if (m1 != 0) net <- net / m1   # relative weights
  nz <- abs(net) > 0
  if (any(nz)) net[abs(net) < low_conn * mean(abs(net[nz]))] <- 0   # drop the weakest regulations
  list(network = net, sign = sgn)
}

# ---- cluster number selection (silhouette over Ward clusterings) --------------------------------
genist_n_clusters <- function(Z, n_tf) {
  if (n_tf <= 10 || ncol(Z) == 1) return(list(k = 1L, hc = NULL))
  hc <- hclust(dist(Z), "ward.D2")                      # Ward linkage on Euclidean distance
  d2 <- as.matrix(dist(Z))^2                            # silhouette uses squared Euclidean distance
  ks <- 2:floor(n_tf / 5)
  ks <- ks[ks < nrow(Z)]
  if (!length(ks)) return(list(k = 1L, hc = hc))
  sil <- vapply(ks, function(k) mean(cluster::silhouette(cutree(hc, k), dmatrix = d2)[, "sil_width"]),
                numeric(1))
  list(k = ks[which.max(sil)], hc = hc)
}

# ---- Cytoscape edge table (matrix_to_cytoscape_table2_0) ----------------------------------------
genist_edge_table <- function(names, network, sign_network) {
  n <- nrow(network); src <- character(); tgt <- character(); w <- numeric(); sg <- numeric()
  for (i in seq_len(n)) {
    if (sum(abs(network[i, ])) + abs(sum(network[, i])) == 0) {   # isolated node: keep it visible
      src <- c(src, names[i]); tgt <- c(tgt, names[i]); w <- c(w, 0); sg <- c(sg, 1)
    }
    for (j in which(network[i, ] != 0)) {
      src <- c(src, names[i]); tgt <- c(tgt, names[j]); w <- c(w, abs(network[i, j])); sg <- c(sg, sign_network[i, j])
    }
  }
  data.frame(source = src, target = tgt, width = w, sign = sg, stringsAsFactors = FALSE)
}

# ---- input checks (CheckInputsGENIST) ------------------------------------------------------------
# Returns a character vector of problems (empty = fine).
genist_check_inputs <- function(genes_file, time_course_file, clustering_file = "", tf_file = "",
                                symbol_file = "", output_folder = "", output_file = "cytoscape_table.txt") {
  p <- character()
  need <- function(path, what) if (!nzchar(path)) p <<- c(p, paste("Please provide", what, ".")) else
    if (!file.exists(path)) p <<- c(p, paste("Provided", what, "cannot be found."))
  need(genes_file, "a genes file"); need(time_course_file, "a time course file")
  for (x in list(list(clustering_file, "clustering file"), list(tf_file, "transcription factors file"),
                 list(symbol_file, "gene name file")))
    if (nzchar(x[[1]]) && !file.exists(x[[1]])) p <- c(p, paste("Provided", x[[2]], "cannot be found."))
  if (!nzchar(output_folder)) p <- c(p, "Please provide an outputs folder.")
  if (!nzchar(output_file)) p <- c(p, "Please provide an output file name.")
  if (nzchar(output_folder) && nzchar(output_file)) {                      # the results must never replace one of the researcher's own files
    ins <- c(`genes file` = genes_file, `time course file` = time_course_file, `clustering data file` = clustering_file,
             `transcription factors file` = tf_file, `gene name file` = symbol_file)
    p <- c(p, output_overwrites_input(ins[!is.na(ins) & nzchar(ins)], file.path(output_folder, c(output_file, "clusters.xlsx", "clusters.csv"))))
  }
  p
}

# ---- Inputs built from Tux (Cuffdiff) results -------------------------------------------------------
# GENIST needs (1) a time course file: the expression of every gene at each time point, earliest first, and (2) a genes file:
# the genes to put in the network. A Tux run holds both: genes.fpkm_tracking has the FPKM of every gene in every group, and
# gene_exp.diff has the group-to-group tests. Nothing in the Tux folder is changed; the two files are written elsewhere.
.san_group <- function(x) gsub("[^A-Za-z0-9._-]+", "_", trimws(x))       # the way Tux step names a group

# Find the Tux results inside `path`: a StRym project folder, its 01_quantification or tux folder, or diff_out itself.
# Returns list(diff_dir, project_dir (NULL if not a StRym project), design (the project's design record, or NULL)).
tux_find_results <- function(path) {
  path <- sub("[/\\\\]+$", "", path)
  if (is.na(path) || !nzchar(path) || !dir.exists(path)) stop("Folder not found.", call. = FALSE)
  cand <- c(path, file.path(path, "diff_out"), file.path(path, "tux", "diff_out"), file.path(path, "01_quantification", "tux", "diff_out"))
  ok <- cand[file.exists(file.path(cand, "genes.fpkm_tracking")) & file.exists(file.path(cand, "gene_exp.diff"))]
  if (!length(ok)) stop("No Tux results found there. Choose your StRym project folder, or the Tux “diff_out” folder ",
                        "(it holds genes.fpkm_tracking and gene_exp.diff).", call. = FALSE)
  d <- normalizePath(ok[1], winslash = "/"); up <- d; design <- NULL; project <- NULL
  for (i in 1:4) {                                                              # diff_out -> tux -> 01_quantification -> project
    up <- dirname(up); f <- file.path(up, "00_design", "design.json")
    if (file.exists(f)) { design <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL); project <- up; break }
  }
  list(diff_dir = d, project_dir = project, design = design)
}

# The groups of a Tux result in time order. `series` is a list with one entry per condition: a named character vector whose
# names are the time point labels (earliest first) and whose values are the group names inside the Tux files. The order comes
# from the time point names typed on the Experimental Design page; without that record it is Tux's own order, which is
# NOT necessarily the time order (`guessed` = TRUE then).
tux_groups <- function(found) {
  hdr <- strsplit(readLines(file.path(found$diff_dir, "genes.fpkm_tracking"), n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  avail <- sub("_FPKM$", "", hdr[grepl("_FPKM$", hdr)])
  if (!length(avail)) stop("genes.fpkm_tracking has no expression columns (…_FPKM).", call. = FALSE)
  series <- list(); tps <- as.character(found$design$timepoint_names); tps <- tps[nzchar(tps)]
  if (length(tps) >= 2) {
    ts <- .san_group(tps)
    pre <- lapply(ts, function(tp) { m <- avail[endsWith(avail, paste0("_", tp))]; substr(m, 1, nchar(m) - nchar(tp) - 1L) })
    for (p in Reduce(intersect, pre)) series[[p]] <- stats::setNames(paste(p, ts, sep = "_"), tps)
  }
  guessed <- !length(series)
  if (guessed) series <- list(`all groups` = stats::setNames(avail, avail))
  ctrl <- if (!is.null(found$design$params$control_pattern)) found$design$params$control_pattern else "control|ctrl|(^|[^a-z])ctl([^a-z]|$)|mock|untreated|wild.?type|^wt$|^ck$|baseline|reference"
  nm <- names(series); not_ctrl <- nm[!grepl(ctrl, nm, ignore.case = TRUE)]
  list(available = avail, series = series, guessed = guessed, default = if (length(not_ctrl)) not_ctrl[1] else nm[1])
}

#' Write a genes file and a time course file for GENIST from a Tux run.
#' @param tux_path  StRym project folder, or the Tux diff_out folder
#' @param output_folder  where the two files are written
#' @param order  time point labels, earliest first (default: the design's time point order)
#' @param series which condition's groups form the time course (default: the treated one)
#' @param q_cutoff,min_log2fc  a gene is chosen when Tux finds it significantly different between two of the time points
#'        (q <= q_cutoff) by at least this |log2 fold change|
#' @param min_fpkm  ... and it is expressed at at least this FPKM at some time point
#' @param max_genes keep at most this many genes (the strongest changes first); 0 / NA = no limit. GENIST's run time grows steeply with
#'        the number of genes and of transcription factors among them (measured on 2 threads: 500 genes ~2 min, 1500 genes > 30 min).
genist_inputs_from_tux <- function(tux_path, output_folder, order = NULL, series = NULL,
                                   q_cutoff = 0.05, min_log2fc = 1, min_fpkm = 1, max_genes = 500) {
  found <- tux_find_results(tux_path); tg <- tux_groups(found)
  s <- if (length(series) == 1 && !is.na(series) && series %in% names(tg$series)) series else tg$default
  grp <- tg$series[[s]]
  order <- trimws(as.character(order)); order <- order[!is.na(order) & nzchar(order)]
  if (length(order)) {
    bad <- setdiff(order, names(grp))
    if (length(bad)) stop("These time points are not in the Tux results: ", paste(bad, collapse = ", "),
                          ". Available: ", paste(names(grp), collapse = ", "), ".", call. = FALSE)
    if (anyDuplicated(order)) stop("A time point is listed twice.", call. = FALSE)
    grp <- grp[order]
  }
  if (length(grp) < 3) stop("GENIST needs a time course of at least 3 time points; this selection has ", length(grp), ".", call. = FALSE)
  if (!is.numeric(q_cutoff) || is.na(q_cutoff) || q_cutoff <= 0 || q_cutoff > 1) stop("The q-value cutoff must be between 0 and 1.", call. = FALSE)
  if (is.na(min_log2fc) || min_log2fc < 0) stop("The minimum fold change must be 0 or more.", call. = FALSE)

  # expression of each gene at each time point
  tr <- utils::read.delim(file.path(found$diff_dir, "genes.fpkm_tracking"), check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  cols <- paste0(grp, "_FPKM")
  if (!"gene_short_name" %in% names(tr) || !all(cols %in% names(tr))) stop("genes.fpkm_tracking does not have the expected columns.", call. = FALSE)
  id <- trimws(as.character(tr$gene_short_name))
  fp <- sapply(tr[cols], function(v) suppressWarnings(as.numeric(v)))
  if (!is.matrix(fp)) fp <- matrix(fp, nrow = nrow(tr))
  keep <- !is.na(id) & nzchar(id) & id != "-"                                    # loci without a gene name cannot be used
  id <- id[keep]; fp <- fp[keep, , drop = FALSE]
  o <- order(-rowMeans(fp, na.rm = TRUE)); id <- id[o]; fp <- fp[o, , drop = FALSE]    # one row per gene name (the most expressed locus)
  first <- !duplicated(id); id <- id[first]; fp <- fp[first, , drop = FALSE]
  o <- order(id, method = "radix"); id <- id[o]; fp <- fp[o, , drop = FALSE]; rownames(fp) <- id
  colnames(fp) <- names(grp)
  if (length(id) < 2) stop("Fewer than 2 genes with a name in genes.fpkm_tracking.", call. = FALSE)

  # genes with a significant, large change between two of the time points
  de <- utils::read.delim(file.path(found$diff_dir, "gene_exp.diff"), check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  need <- c("gene", "sample_1", "sample_2", "status", "log2(fold_change)", "q_value")
  if (!all(need %in% names(de))) stop("gene_exp.diff does not have the expected columns.", call. = FALSE)
  lf <- suppressWarnings(as.numeric(de[["log2(fold_change)"]])); q <- suppressWarnings(as.numeric(de$q_value))
  nm <- trimws(as.character(de$gene))
  hit <- de$sample_1 %in% grp & de$sample_2 %in% grp & de$status == "OK" & .isT(q <= q_cutoff) & .isT(abs(lf) >= min_log2fc) & !is.na(nm) & nzchar(nm) & nm != "-"
  if (!any(de$sample_1 %in% grp & de$sample_2 %in% grp)) stop("gene_exp.diff has no comparison between these time points.", call. = FALSE)
  strength <- if (any(hit)) tapply(pmin(abs(lf[hit]), 20), nm[hit], max) else numeric()   # strongest change of each gene (Inf capped)
  n_sig <- length(strength)
  strength <- strength[names(strength) %in% id]
  if (length(strength)) strength <- strength[apply(fp[names(strength), , drop = FALSE], 1, function(r) any(r >= min_fpkm, na.rm = TRUE))]
  n_pass <- length(strength)
  capped <- !is.na(max_genes) && max_genes > 0 && n_pass > max_genes
  if (capped) strength <- strength[order(-strength, names(strength), method = "radix")][seq_len(max_genes)]
  genes <- if (length(strength)) c_sort(names(strength)) else character()
  if (length(genes) < 2) stop(sprintf("Only %d gene(s) pass these filters. Relax the q-value cutoff, the fold change or the expression limit.", length(genes)), call. = FALSE)

  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  genes_file <- file.path(output_folder, "genist_genes_from_tux.csv"); tc_file <- file.path(output_folder, "genist_time_course_from_tux.csv")
  utils::write.csv(data.frame(gene = genes), genes_file, row.names = FALSE)
  utils::write.csv(data.frame(gene = id, as.data.frame(fp, check.names = FALSE), check.names = FALSE), tc_file, row.names = FALSE)
  summary <- sprintf("%s genes chosen (q ≤ %s and |log2 fold change| ≥ %s between two of the time points, FPKM ≥ %s at some time point)%s. Time course: %s genes × %d time points (%s).",
                     format(length(genes), big.mark = ","), format(q_cutoff), format(min_log2fc), format(min_fpkm),
                     if (capped) sprintf("; kept the %s strongest of %s", format(max_genes, big.mark = ","), format(n_pass, big.mark = ",")) else "",
                     format(length(id), big.mark = ","), ncol(fp), paste(colnames(fp), collapse = " → "))
  list(genes_file = genes_file, time_course_file = tc_file, n_genes = length(genes), n_pass = n_pass, n_significant = n_sig, capped = capped,
       time_points = colnames(fp), groups = unname(grp), series = s, guessed = tg$guessed, summary = summary)
}

# ---- main ------------------------------------------------------------------------------------------
#' @param genes_file,time_course_file,clustering_file,tf_file,symbol_file  Excel/CSV/TXT files (first row = header,
#'        first column = gene ids). clustering/tf/symbol files are optional ("" or NULL to skip).
#' @param time_lapse  1 (regulators act at the next time point), 0 (same time point) or "both"
#' @param fc_thresh   minimum fold change between consecutive time points to count as a change (default 1.3)
#' @param time_percent minimum fraction of time points with matching changes (default 0.3)
#' @param n_levels    discretization levels, integer >= 2 (default 2)
#' @param low_conn    fraction of the mean edge weight below which edges are dropped (default 0.2)
#' @param cores       parallel workers for the inference step
genist_run <- function(genes_file, time_course_file, clustering_file = NULL, tf_file = NULL, symbol_file = NULL,
                       output_folder, output_file = "cytoscape_table.txt",
                       time_lapse = 1, fc_thresh = 1.3, time_percent = 0.3, n_levels = 2, low_conn = 0.2,
                       cores = 1) {
  blank <- function(x) is.null(x) || is.na(x) || !nzchar(x)
  tl <- if (identical(time_lapse, "both") || identical(time_lapse, c(0, 1))) -1 else as.numeric(time_lapse)
  if (!tl %in% c(0, 1, -1)) stop("time_lapse must be 0, 1 or \"both\"", call. = FALSE)
  if (fc_thresh < 1) stop("fc_thresh must be >= 1", call. = FALSE)
  if (n_levels < 2 || n_levels != round(n_levels)) stop("n_levels must be an integer >= 2", call. = FALSE)
  problems <- genist_check_inputs(genes_file, time_course_file, if (blank(clustering_file)) "" else clustering_file,
                                  if (blank(tf_file)) "" else tf_file, if (blank(symbol_file)) "" else symbol_file,
                                  output_folder, output_file)
  if (length(problems)) stop(paste(problems, collapse = " "), call. = FALSE)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

  emit_progress(2, "Reading input files")
  genes <- c_unique_sort(read_id_column(genes_file))
  tc <- note_dropped_columns(table_to_matrix(read_table_any(time_course_file)), "time course file")
  has_clust <- !blank(clustering_file)
  clu <- if (has_clust) note_dropped_columns(table_to_matrix(read_table_any(clustering_file)), "clustering data file") else NULL
  tf_list <- if (blank(tf_file)) genes else read_id_column(tf_file)
  sym_df <- if (blank(symbol_file)) NULL else read_table_any(symbol_file, header = FALSE)   # name files often have no header row

  # genes present in every needed file
  keep <- genes[genes %in% rownames(tc)]
  if (has_clust) keep <- keep[keep %in% rownames(clu)]
  dropped <- setdiff(genes, keep)
  if (length(dropped)) message("Note: ", length(dropped), " of ", length(genes),
                               " genes are missing from the time course", if (has_clust) "/clustering",
                               " file and were skipped.")
  if (length(keep) < 2) stop("Fewer than 2 genes of the genes file were found in the data files.", call. = FALSE)
  genes <- keep
  if (ncol(tc) < 3) stop("The time course file needs at least 3 time points (found ", ncol(tc), ").", call. = FALSE)
  tcm <- tc[genes, , drop = FALSE]
  if (anyNA(tcm)) message("Note: the time course contains missing values; comparisons involving them count as no change.")

  emit_progress(8, "Clustering genes")
  n_tf <- sum(genes %in% tf_list)
  if (has_clust && ncol(clu) > 1) {
    cd <- clu[genes, , drop = FALSE]
    Z <- (cd - rowMeans(cd)) / (apply(cd, 1, sd) + 1e-5)
  } else Z <- matrix(1, length(genes), 1, dimnames = list(genes, NULL))
  nc <- genist_n_clusters(Z, n_tf)
  K <- nc$k
  if (K > 1) {
    idx <- cutree(nc$hc, K)
    members <- lapply(seq_len(K), function(i) genes[idx == i])
  } else members <- list(genes)
  message(sprintf("%d genes, %d transcription factors, %d cluster%s", length(genes), n_tf, K, if (K > 1) "s" else ""))

  # cluster lists (clusters.xlsx: sheet 1 = ids, sheet 2 = symbols)
  maxn <- max(lengths(members))
  pad <- function(v) c(v, rep("", maxn - length(v)))
  ids_tab <- as.data.frame(lapply(members, pad), stringsAsFactors = FALSE, col.names = paste0("Cluster_", seq_len(K)))
  sym_tab <- as.data.frame(lapply(members, function(v) pad(map_symbols(v, sym_df))), stringsAsFactors = FALSE,
                           col.names = paste0("Cluster_", seq_len(K)))
  if (requireNamespace("writexl", quietly = TRUE))
    writexl::write_xlsx(list(Gene_IDs = ids_tab, Gene_symbols = sym_tab), file.path(output_folder, "clusters.xlsx"))
  else utils::write.csv(ids_tab, file.path(output_folder, "clusters.csv"), row.names = FALSE)

  # one GRN per cluster
  names_all <- character(); blocks <- list(); hubs <- character()
  for (ci in seq_len(K)) {
    emit_progress(10 + 70 * (ci - 1) / K, sprintf("Inferring network for cluster %d of %d (%d genes)", ci, K, length(members[[ci]])))
    g <- members[[ci]]; ng <- length(g)
    if (ng > 1) {
      pre <- genist_preprocess(tcm[g, , drop = FALSE], fc_thresh, time_percent, n_levels, tf_pos = g %in% tf_list)
      res <- genist_infer(pre, tl, low_conn, cores)
      net <- res$network; sgn <- res$sign
    } else { net <- matrix(0, 1, 1); sgn <- matrix(0, 1, 1) }
    blocks[[ci]] <- list(net = net, sgn = sgn)
    names_all <- c(names_all, g)
    nconn <- rowSums(net != 0)                      # out-connections of each gene
    hubs <- c(hubs, g[nconn == max(nconn)])         # all ties; a cluster without edges contributes every gene
  }
  N <- length(names_all); final <- matrix(0, N, N); fsgn <- matrix(0, N, N); pos <- 0L
  for (b in blocks) { r <- pos + seq_len(nrow(b$net)); final[r, r] <- b$net; fsgn[r, r] <- b$sgn; pos <- pos + nrow(b$net) }

  # connect the clusters through their hubs
  if (K > 1) {
    emit_progress(85, "Connecting clusters through their hubs")
    hub_genes <- genes[genes %in% hubs]
    if (length(hub_genes) > 1) {
      pre <- genist_preprocess(tcm[hub_genes, , drop = FALSE], fc_thresh, time_percent, n_levels)   # no TF filter, as in MATLAB
      res <- genist_infer(pre, tl, low_conn, cores)
      p <- match(hub_genes, names_all)
      final[p, p] <- final[p, p] + res$network
      fsgn[p, p] <- fsgn[p, p] + res$sign
    }
  }

  emit_progress(95, "Writing results")
  edges <- genist_edge_table(map_symbols(names_all, sym_df), final, fsgn)
  out_path <- file.path(output_folder, output_file)
  utils::write.table(edges, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  e <- edges[edges$width > 0, , drop = FALSE]                                           # a quick picture for the interactive browser (never stops the run)
  strym_network_png(data.frame(from = e$source, to = e$target, type = ifelse(e$sign > 0, "activates", ifelse(e$sign < 0, "inhibits", "regulates")), weight = e$width,
                               stringsAsFactors = FALSE),
                    sub("\\.[A-Za-z0-9]+$", "_network.png", out_path), sprintf("GENIST network: %d genes, %d regulations", length(unique(c(e$source, e$target))), nrow(e)))
  emit_progress(100, "Done")
  invisible(list(edges = edges, clusters = members, output = out_path, network = final, sign = fsgn, genes = names_all))
}

# ---- faster inference (default) ---------------------------------------------------------------------------------------
# genist_infer_fast(): behaviour-preserving faster replacement for genist_infer_orig() (written and verified in the overnight audit: identical network/sign matrices on 48 randomised cases).
# Same signature / result. Idea: the BDe score of a candidate set only needs lgamma() of small integer counts (0..Tn), so
# (1) look the terms up in a table instead of calling lgamma per cell (identical doubles), and
# (2) score ALL (kk, kkk) partners of a given first regulator k in one vectorised batch (tabulate on a matrix),
# keeping the original candidate order for tie-breaking (the LAST candidate with the maximal score wins, as in the MATLAB code).
genist_infer_fast <- function(pre, time_lapse, low_conn, cores = 1, chunk_cells = 2e6) {
  disc <- pre$disc; up <- pre$up; dn <- pre$dn; L <- pre$n_levels
  n <- nrow(disc); Tn <- ncol(disc); set_max_reg <- pre$set_max_reg
  ESS <- 1e-15
  tabs <- lapply(1:3, function(m) { q <- L^m; a <- ESS / q; b <- ESS / (L * q)
    list(q = q, A = lgamma(a) - lgamma(0:Tn + a), B = lgamma(0:Tn + b) - lgamma(b)) })

  # scores of the columns of IDX (config index 0..q-1 per time point), target states tv (length nrow(IDX)); m regulators
  score_batch <- function(IDX, m, tv) {
    tb <- tabs[[m]]; q <- tb$q; nc <- ncol(IDX); Tc <- nrow(IDX)
    out <- numeric(nc); step <- max(1L, as.integer(chunk_cells %/% (q * L)))
    for (s0 in seq(1L, nc, by = step)) {
      cols <- s0:min(nc, s0 + step - 1L); k <- length(cols); X <- IDX[, cols, drop = FALSE]
      Nj  <- tabulate(X + rep((0:(k - 1L)) * q, each = Tc) + 1, q * k)
      Njk <- tabulate(X * L + tv + rep((0:(k - 1L)) * q * L, each = Tc) + 1, q * L * k)
      out[cols] <- colSums(matrix(tb$A[Nj + 1L], q)) + colSums(matrix(tb$B[Njk + 1L], q * L))
    }
    out
  }
  seq_cache <- new.env()                                   # candidate order inside one k: for kk in 1..k-1: pair(kk), triples(kk, 1..kk-1)
  seq_for_k <- function(k) {
    key <- as.character(k); if (!is.null(h <- seq_cache[[key]])) return(h)
    kk <- if (k > 1) rep(seq_len(k - 1L), times = seq_len(k - 1L)) else integer()
    kkk <- if (k > 1) unlist(lapply(seq_len(k - 1L), function(x) 0:(x - 1L))) else integer()   # 0 = the pair itself
    assign(key, list(kk = kk, kkk = kkk), envir = seq_cache); seq_cache[[key]]
  }

  per_layer <- function(R, shift) {
    cE <- seq_len(Tn - shift); ct <- (1L + shift):Tn
    res <- parallel::mclapply(seq_len(n), function(j) {
      ind3 <- which(R[, j]); nreg <- length(ind3)
      if (nreg == 0L) return(list(best = -Inf, idx = integer(), sign = numeric()))      # no candidate regulators (NULL is reserved for a failed worker)
      tv <- disc[j, ct]; tu <- up[j, ct]; td <- dn[j, ct]
      Ed <- t(disc[ind3, cE, drop = FALSE])                # time x regulators
      best <- -Inf; sel <- NULL
      for (k in nreg:1) {
        ek <- Ed[, k]
        sc <- score_batch(matrix(ek, ncol = 1), 1L, tv)    # the single regulator comes first
        if (nreg > 1 && set_max_reg > 1 && k > 1) {
          sq <- seq_for_k(k); np <- length(sq$kk); seqsc <- numeric(np)
          isp <- sq$kkk == 0L
          seqsc[isp] <- score_batch(ek + L * Ed[, sq$kk[isp], drop = FALSE], 2L, tv)
          if (nreg > 2 && set_max_reg > 2 && any(!isp))
            seqsc[!isp] <- score_batch(ek + L * Ed[, sq$kk[!isp], drop = FALSE] + L^2 * Ed[, sq$kkk[!isp], drop = FALSE], 3L, tv)
          sc <- c(sc, seqsc)
        }
        M <- max(sc)
        if (M >= best && M != -Inf) {                      # ties replace the earlier candidate
          best <- M; w <- which(sc == M); w <- w[length(w)]
          sel <- if (w == 1L) ind3[k] else { sq <- seq_for_k(k); p <- w - 1L
            if (sq$kkk[p] == 0L) ind3[c(k, sq$kk[p])] else ind3[c(k, sq$kk[p], sq$kkk[p])] }
        }
      }
      if (is.null(sel)) return(list(best = best, idx = integer(), sign = numeric()))
      sg <- vapply(sel, function(r) { eu <- up[r, cE]; ed <- dn[r, cE]
        sign(sum((eu & tu) | (ed & td)) - sum((eu & td) | (ed & tu))) }, numeric(1))
      list(best = best, idx = sel, sign = sg)
    }, mc.cores = cores)
    if (any(vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)))
      stop("A GENIST worker process failed (probably out of memory). Run again with fewer threads.", call. = FALSE)
    net <- matrix(0, n, n); sgn <- matrix(0, n, n)
    for (j in seq_len(n)) {
      r <- res[[j]]
      if (is.null(r) || !length(r$idx)) next
      net[r$idx, j] <- 10^r$best
      sgn[r$idx, j] <- r$sign
    }
    list(net = net, sgn = sgn)
  }

  net <- NULL
  if (time_lapse != 0) l1 <- per_layer(pre$reg1, 1L)
  if (time_lapse != 1) l2 <- per_layer(pre$reg2, 0L)
  if (time_lapse == 0) { net <- l2$net; sgn <- l2$sgn
  } else if (time_lapse == -1) {
    net <- l1$net + l2$net; sgn <- l1$sgn + l2$sgn
    sgn[sgn == 2] <- 1; sgn[sgn == -2] <- -1
  } else { net <- l1$net; sgn <- l1$sgn }

  m1 <- max(abs(net)); if (m1 != 0) net <- net / m1
  nz <- abs(net) > 0
  if (any(nz)) net[abs(net) < low_conn * mean(abs(net[nz]))] <- 0
  list(network = net, sign = sgn)
}

# The fast version is used unless options(strym.genist_orig = TRUE) asks for the original loop (kept for comparison / bit-for-bit reproduction).
genist_infer <- function(pre, time_lapse, low_conn, cores = 1) {
  if (isTRUE(getOption("strym.genist_orig"))) genist_infer_orig(pre, time_lapse, low_conn, cores) else genist_infer_fast(pre, time_lapse, low_conn, cores)
}
