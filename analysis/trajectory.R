# analysis/trajectory.R -- expression trajectories over time, treated vs control (from the baseline simple_trajectories.R,
# generalized). One panel row per gene, one column per treated-vs-control pair; line = mean, points = replicates.
#
# Which genes: (1) the gene list file chosen in Analysis settings (one id or name per line), else (2) the main-set genes of the
# entropy step(s) run in this project, strongest first. Gene ids or names (e.g. AT1G01010) are both understood.
# Input  : Tux results (FPKM per replicate) or a counts table (counts per million), with >= 2 time points.
# Output : trajectories.png (the first trajectory_max genes), trajectories.pdf (all of them), trajectory_values.csv (means).
.read_gene_list <- function(path) {
  ext <- tolower(tools::file_ext(path))
  g <- if (ext %in% c("xlsx", "xls")) read_id_column(path) else {
    x <- trimws(readLines(path, warn = FALSE)); x <- x[nzchar(x)]
    x <- trimws(sub("[,\t;].*$", "", x)); gsub("^[\"']|[\"']$", "", x)
  }
  g <- g[nzchar(g)]
  if (length(g) && grepl("^(genes?|gene_?id|ids?|names?|symbols?)$", g[1], ignore.case = TRUE)) g <- g[-1]   # a header line
  unique(g)
}

# Genes to plot and where they came from. Returns list(genes, source).
.trajectory_genes <- function(ctx, p) {
  if (nzchar(p$trajectory_genes)) {
    if (!file.exists(p$trajectory_genes)) stop("The gene list file for trajectories was not found: ", p$trajectory_genes, call. = FALSE)
    g <- .read_gene_list(p$trajectory_genes)
    if (!length(g)) stop("The gene list file for trajectories is empty.", call. = FALSE)
    return(list(genes = g, source = "your gene list"))
  }
  found <- list()
  for (i in names(ctx$prereq)) {
    f <- list.files(ctx$prereq[[i]], pattern = "_main_genes\\.csv$", full.names = TRUE)
    for (ff in f) { x <- tryCatch(utils::read.csv(ff, stringsAsFactors = FALSE), error = function(e) NULL)
                    if (!is.null(x) && nrow(x) && "gene" %in% names(x)) found[[length(found) + 1]] <- x[, intersect(c("gene", "n_sig_timepoints"), names(x)), drop = FALSE] }
  }
  if (!length(found)) stop("No genes to plot: choose a gene list file in Analysis settings, or run an entropy step that finds genes (its main set is plotted).", call. = FALSE)
  x <- do.call(rbind, lapply(found, function(d) { if (!"n_sig_timepoints" %in% names(d)) d$n_sig_timepoints <- 0; d }))
  x <- x[order(-x$n_sig_timepoints), ]
  list(genes = unique(x$gene), source = "the entropy results")
}

# Gene ids in the data for the requested ids / names (NA = not found).
.match_genes <- function(req, data) {
  ids <- strym_gene_ids(data); id <- ids[match(tolower(req), tolower(ids))]
  nm <- data$gene_names
  if (!is.null(nm) && anyNA(id)) {
    lst <- strsplit(ifelse(is.na(nm), "", nm), "[,;]\\s*"); key <- tolower(unlist(lst)); owner <- rep(names(nm), lengths(lst))
    miss <- is.na(id); id[miss] <- owner[match(tolower(req[miss]), key)]
  }
  data.frame(request = req, id = id, stringsAsFactors = FALSE)
}

strym_step_trajectory <- function(ctx) {
  p <- strym_normalize_params(ctx$design$params)
  data <- strym_drop_samples(strym_load_expression(ctx), p)
  S <- data$samples
  tps <- strym_order_timepoints(S$time_point, ctx$design$timepoint_names)
  if (length(tps) < 2) stop("Trajectories need at least 2 time points in the sample sheet.", call. = FALSE)
  want <- .trajectory_genes(ctx, p)
  mt <- .match_genes(want$genes, data); mt <- mt[!duplicated(mt$id) | is.na(mt$id), ]
  missing <- mt$request[is.na(mt$id)]; mt <- mt[!is.na(mt$id), , drop = FALSE]
  if (!nrow(mt)) stop(sprintf("None of the %d genes from %s were found in the data (they look like: %s).", length(want$genes), want$source,
                              paste(utils::head(strym_gene_ids(data), 3), collapse = ", ")), call. = FALSE)
  if (length(missing)) ctx$log(sprintf("Trajectories: %d of %d genes not found in the data: %s%s", length(missing), length(want$genes),
                                       paste(utils::head(missing, 8), collapse = ", "), if (length(missing) > 8) ", ..." else ""))
  all_ids <- mt$id; mt <- utils::head(mt, 200)                                    # PDF holds up to 200 genes, the PNG the first trajectory_max
  nm <- data$gene_names
  label <- function(id) { n <- if (!is.null(nm)) unname(nm[id]) else NA_character_; ifelse(is.na(n) | n == id, id, paste0(n, " (", id, ")")) }
  r <- strym_long(data, mt$id)
  r <- merge(r, S[, c("sample", "condition", "time_point", "replicate")], by = "sample")
  r$time_point <- factor(r$time_point, levels = tps); r <- r[!is.na(r$time_point), ]
  r$label <- factor(label(r$gene), levels = unique(label(mt$id)))

  pairs <- tryCatch(strym_pairs(S, p), error = function(e) NULL)                  # no control found: just one line per condition
  if (!is.null(pairs)) {
    r <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
      x <- r[r$condition %in% c(pairs$treated[i], pairs$control[i]), ]
      if (!nrow(x)) return(NULL)
      x$panel <- paste0(pairs$treated[i], " vs ", pairs$control[i]); x$role <- ifelse(x$condition == pairs$control[i], "Control", x$condition); x
    }))
    cols <- c(Control = "#808080", setNames(rep(c("#4156A1", "#6F7D1C", "#008473", "#344891", "#586800", "#427E93"), length.out = length(unique(pairs$treated))), unique(pairs$treated)))
  } else {
    r$panel <- "all conditions"; r$role <- r$condition
    cols <- setNames(rep(c("#4156A1", "#6F7D1C", "#808080", "#008473", "#344891", "#586800"), length.out = length(unique(r$condition))), unique(r$condition))
  }
  r$role <- factor(r$role, levels = names(cols)[names(cols) %in% r$role]); r$panel <- factor(r$panel, levels = unique(r$panel))
  agg <- stats::aggregate(value ~ gene + label + panel + role + time_point, r, mean)
  utils::write.csv(agg[order(agg$label, agg$panel, agg$role, agg$time_point), ], file.path(ctx$out_dir, "trajectory_values.csv"), row.names = FALSE)

  plot_it <- function(df, n) {
    ggplot2::ggplot(df, ggplot2::aes(time_point, value, colour = role, group = role)) +
      ggplot2::geom_point(alpha = 0.4) + ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 1) +
      ggplot2::facet_grid(label ~ panel, scales = "free_y") + ggplot2::scale_colour_manual(values = cols) +
      ggplot2::expand_limits(y = 0) + ggplot2::labs(x = "Time point", y = data$unit, colour = NULL) + ggplot2::theme_bw() +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, size = 7))
  }
  top <- utils::head(levels(r$label), p$trajectory_max)
  png_df <- r[r$label %in% top, ]
  ggplot2::ggsave(file.path(ctx$out_dir, "trajectories.png"), plot_it(png_df), width = 3 + 2.6 * nlevels(r$panel), height = 1.5 * length(top) + 1, dpi = 120, limitsize = FALSE)
  ggplot2::ggsave(file.path(ctx$out_dir, "trajectories.pdf"), plot_it(r), width = 3 + 2.6 * nlevels(r$panel), height = 1.5 * nlevels(r$label) + 1, limitsize = FALSE)
  msg <- sprintf("Plotted %d of %d genes from %s%s.", nlevels(r$label), length(want$genes), want$source,
                 if (nlevels(r$label) > length(top)) sprintf(" (picture: first %d; the PDF has all)", length(top)) else "")
  ctx$log(paste("Trajectories:", msg))
  list(status = "ok", files = c("trajectories.png", "trajectories.pdf", "trajectory_values.csv"), message = msg)
}
