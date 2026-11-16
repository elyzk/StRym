# analysis/qc_overview.R -- a quick look at the DATA before (or alongside) any analysis: is the experiment complete, are the samples comparable?
#
# Real experiments are untidy (a failed sample, a time point only some conditions have, one library ten times smaller than the rest). This
# step draws four small pictures and writes the problems it sees in plain words, so they are noticed before they become a surprising result:
#   design_grid.png         samples per condition x time point (empty cells in red: the experiment has a hole there; 1 = no replicate)
#   sample_depth.png        reads (counts) per sample, with the samples far below the typical depth marked
#   genes_detected.png      genes with any signal, per sample
#   sample_correlation.png  how alike the samples are (top variable genes), with condition / time point beside each sample: replicates should
#                           sit together; a sample that does not is an outlier you may want to leave out (Analysis settings)
#   data_overview_notes.csv every observation, one per line (issue, where, what to do)
# Input  : an expression table (counts / FPKM), long-read counts, or a Tux run. Never stops a run: if something cannot be drawn it says so.
strym_step_qc_overview <- function(ctx)
  tryCatch(.qc_overview_run(ctx), error = function(e) list(status = "skipped", files = character(),
                                                           message = paste("The data overview could not be made:", conditionMessage(e))))

.qc_overview_run <- function(ctx) {
  p <- strym_normalize_params(ctx$design$params); d <- ctx$design; out <- ctx$out_dir
  data <- strym_drop_samples(strym_load_expression(ctx), p)
  S <- data$samples; M <- strym_wide(data); files <- character(); notes <- list()
  note <- function(level, what, where, advice) notes[[length(notes) + 1]] <<- data.frame(level = level, issue = what, where = where, what_to_do = advice, stringsAsFactors = FALSE)
  S$time_point[is.na(S$time_point) | !nzchar(S$time_point)] <- "all"
  tps <- strym_order_timepoints(S$time_point, d$timepoint_names); conds <- unique(S$condition)
  S$time_point <- factor(S$time_point, levels = tps); S$condition <- factor(S$condition, levels = conds)
  S <- S[order(S$condition, S$time_point, S$replicate), , drop = FALSE]
  theme <- ggplot2::theme_bw(base_size = 11)
  save <- function(g, f, w, h) { ggplot2::ggsave(file.path(out, f), g, width = min(w, 40), height = min(h, 40), dpi = 110, limitsize = FALSE); files <<- c(files, f) }
  section <- function(what, expr) tryCatch(expr, error = function(e) note("info", sprintf("%s could not be drawn (%s)", what, conditionMessage(e)), "overview", "The other pictures and the notes are not affected."))
  dr <- data$dropped                                                                    # what the table reader left out: a failed library must not vanish silently
  if (length(dr$empty)) note("warning", sprintf("Sample(s) with no counts at all were left out: %s", paste(utils::head(dr$empty, 6), collapse = ", ")), paste(utils::head(dr$empty, 6), collapse = ", "), "A library with no reads is a failed sample: check the sequencing or the export of the table.")
  if (length(dr$unlisted)) note("warning", sprintf("%d column(s) of the table are not in the sample sheet and were left out: %s", length(dr$unlisted), paste(utils::head(dr$unlisted, 6), collapse = ", ")), paste(utils::head(dr$unlisted, 6), collapse = ", "), "Add them to the sample sheet if they are samples.")
  if (length(dr$missing)) note("info", sprintf("%d sample(s) of the sample sheet are not in the table: %s", length(dr$missing), paste(utils::head(dr$missing, 6), collapse = ", ")), "sample sheet", "Check the spelling of the names.")
  # what the table reader repaired on the way in: said here once, in the notes file and the summary page, not only in the log
  if (isTRUE(dr$rows_blank > 0) || isTRUE(dr$rows_dup > 0)) note("info", sprintf("%d row(s) of the table were left out: %d without a gene id, %d with a gene id that was already listed (the first is kept)", dr$rows_blank + dr$rows_dup, dr$rows_blank, dr$rows_dup), "table", "Check the gene id column of your table.")
  if (isTRUE(dr$na_cells > 0)) note("warning", sprintf("%s empty or non-numeric cell(s) were read as 0", format(dr$na_cells, big.mark = ",")), "table", "Check the table for text, '-' or 'NA' where a number belongs.")
  if (isTRUE(dr$inf_cells > 0)) note("warning", sprintf("%d cell(s) held Inf / -Inf and were read as 0", dr$inf_cells), "table", "Inf is not a count: check where it came from (a division by zero in a spreadsheet?).")
  if (isTRUE(dr$decimal_cells > 0)) note("info", sprintf("%d cell(s) of a table of whole numbers held decimals and were rounded", dr$decimal_cells), "table", "Check that no value was typed with a decimal comma by mistake.")
  if (isTRUE(data$whole_numbers) && is.null(data$counts) && identical(d$data_kind, "counts"))
    note("info", "Every value is a whole number although the table is set as already normalized (FPKM / TPM)", "table", "If these are raw counts, choose “Raw counts” as the unit: DESeq2 then becomes available and the depth pictures show reads.")
  vv <- as.numeric(M[M > 0]); if (length(vv) > 100) { q999 <- stats::quantile(vv, 0.999); mx <- max(vv)
    if (is.finite(q999) && q999 > 0 && mx > 1000 * q999) { w <- which(M == mx, arr.ind = TRUE)[1, ]
      note("warning", sprintf("Unusually large value: %s in sample %s is %s (over 1000 times the 99.9th percentile of the table)", rownames(M)[w[1]], colnames(M)[w[2]], format(mx, big.mark = ",", scientific = FALSE)), colnames(M)[w[2]],
           "One huge value dominates its sample (it shrinks every other gene in the per-million scaling and distorts PCA). Check for a typo such as extra digits.") } }

  # 1. the design grid
  grid <- as.data.frame(table(condition = S$condition, time_point = S$time_point)); names(grid)[3] <- "samples"
  utils::write.csv(grid, file.path(out, "design_grid.csv"), row.names = FALSE); files <- c(files, "design_grid.csv")
  section("The design grid", {
  g <- ggplot2::ggplot(grid, ggplot2::aes(time_point, condition, fill = factor(pmin(samples, 3), levels = 0:3, labels = c("none", "1", "2", "3 or more")))) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) + ggplot2::geom_text(ggplot2::aes(label = samples), size = 4.2) +
    ggplot2::scale_fill_manual(values = c(none = "#E8A09A", `1` = "#F3D9A4", `2` = "#B5D6B2", `3 or more` = "#7FBF7B"), name = "samples", drop = FALSE) +
    ggplot2::labs(x = "Time point", y = NULL, title = "Samples per condition and time point") + theme + ggplot2::theme(panel.grid = ggplot2::element_blank())
  save(g, "design_grid.png", max(4.5, 1.1 * length(tps) + 2.5), max(2.5, 0.55 * length(conds) + 1.6)) })
  if (length(conds) > 1 && length(tps) > 1) {
    holes <- grid[grid$samples == 0, , drop = FALSE]
    for (i in seq_len(nrow(holes))) if (any(grid$samples[grid$time_point == holes$time_point[i]] > 0))
      note("warning", sprintf("No samples for %s at %s", holes$condition[i], holes$time_point[i]), "design", "Time-course comparisons leave this time point out. If the sample failed, check the sample sheet; if it was never taken, this is expected.")
  }
  two <- grid[grid$samples == 2, , drop = FALSE]
  if (nrow(two) && !nrow(grid[grid$samples >= 3, , drop = FALSE])) note("info", sprintf("Every condition/time point combination with samples has only 2 replicates (%d of them)", nrow(two)), "design", "Two replicates are the minimum for the statistics: they work, but only large changes will be found. Three or more replicates are better.")
  one <- grid[grid$samples == 1, , drop = FALSE]
  if (nrow(one)) note("info", sprintf("%d condition/time point combination(s) have a single sample (e.g. %s)", nrow(one), if (identical(as.character(one$time_point[1]), "all")) as.character(one$condition[1]) else sprintf("%s at %s", one$condition[1], one$time_point[1])), "design", "Without replicates there is no estimate of variability: the statistics are weak there.")

  # 2. depth and genes detected, from raw counts when there are some
  raw <- data$counts
  if (is.null(raw) && identical(d$data_kind, "fastq")) { cp <- strym_counts_path(d)
    if (file.exists(cp)) raw <- tryCatch(.counts_raw(cp, "")$m, error = function(e) NULL) }
  src <- if (!is.null(raw)) raw[, intersect(colnames(raw), S$sample), drop = FALSE] else M[, S$sample, drop = FALSE]
  lab <- if (!is.null(raw)) "Reads counted per sample (millions)" else sprintf("Total %s per sample", data$unit)
  tot <- colSums(src); scale <- if (!is.null(raw)) 1e6 else 1
  section("The depth picture", {
  dep <- data.frame(sample = factor(names(tot), levels = S$sample[S$sample %in% names(tot)]), total = tot / scale, S[match(names(tot), S$sample), c("condition", "time_point")])
  g <- ggplot2::ggplot(dep, ggplot2::aes(sample, total, fill = condition)) + ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.2 * stats::median(dep$total), linetype = 2, colour = "#B03A2E") +
    ggplot2::labs(x = NULL, y = lab, title = "Sequencing depth per sample", subtitle = "dashed line = 20% of the median sample") + theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = if (nrow(dep) > 60) 5 else 8))
  save(g, "sample_depth.png", max(5.5, 0.22 * nrow(dep) + 2), 4) })
  dep <- data.frame(sample = factor(names(tot), levels = S$sample[S$sample %in% names(tot)]), total = tot / scale, S[match(names(tot), S$sample), c("condition", "time_point")])
  low <- names(tot)[tot < 0.2 * stats::median(tot)]
  if (length(low)) note("warning", sprintf("Low depth: %s (under 20%% of the median sample)", paste(utils::head(low, 5), collapse = ", ")), paste(utils::head(low, 5), collapse = ", "),
                        "Normalization copes with moderate differences, but a very small library is noisy: consider leaving it out (Analysis settings: samples to leave out).")
  det <- data.frame(sample = factor(colnames(src), levels = levels(dep$sample)), genes = colSums(src > 0), condition = S$condition[match(colnames(src), S$sample)], stringsAsFactors = FALSE)
  section("The genes-detected picture", {
  g <- ggplot2::ggplot(det, ggplot2::aes(sample, genes, fill = condition)) + ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = "Genes with any signal", title = "Genes detected per sample") + theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = if (nrow(det) > 60) 5 else 8))
  save(g, "genes_detected.png", max(5.5, 0.22 * nrow(det) + 2), 4) })
  lowdet <- as.character(det$sample[det$genes < 0.6 * stats::median(det$genes)])
  if (length(lowdet)) note("warning", sprintf("Few genes detected: %s (under 60%% of the median)", paste(utils::head(lowdet, 5), collapse = ", ")), paste(utils::head(lowdet, 5), collapse = ", "), "A sample with far fewer detected genes than its replicates is usually a failed library.")
  dead <- sum(rowSums(src) == 0); if (dead) note("info", sprintf("%s of %s genes have no signal in any sample", format(dead, big.mark = ","), format(nrow(src), big.mark = ",")), "genes", "Harmless: they are left out of the statistics.")

  # 3. how alike are the samples?
  if (ncol(M) >= 3 && requireNamespace("pheatmap", quietly = TRUE)) {
    L <- log2(M[, S$sample, drop = FALSE] + 1); L <- L[rowSums(L) > 0, , drop = FALSE]
    v <- rowSums((L - rowMeans(L))^2) / max(1, ncol(L) - 1); top <- utils::head(order(-v), min(2000, length(v))); cm <- stats::cor(L[top, , drop = FALSE])
    ann <- data.frame(condition = S$condition, row.names = S$sample); if (length(tps) > 1) ann$time_point <- S$time_point
    n <- ncol(cm); f <- file.path(out, "sample_correlation.png"); sz <- max(5.5, min(14, 0.3 * n + 3))
    brk <- if (is.finite(diff(range(cm, na.rm = TRUE))) && diff(range(cm, na.rm = TRUE)) > 1e-8) seq(min(cm, na.rm = TRUE), 1, length.out = 101) else NA      # identical samples: no colour scale to cut
    section("The sample correlation picture", suppressWarnings(pheatmap::pheatmap(cm, annotation_col = ann, annotation_row = ann, show_rownames = n <= 60, show_colnames = n <= 60, fontsize = if (n > 40) 6 else 9,
                                        color = grDevices::colorRampPalette(c("#2C5F8A", "#F7F7F7", "#B2182B"))(100), breaks = brk,
                                        main = "Sample-to-sample correlation (top variable genes)", filename = f, width = sz, height = sz * 0.9)))
    if (file.exists(f)) files <- c(files, "sample_correlation.png")
    off <- vapply(seq_len(n), function(i) mean(cm[i, -i]), 0); z <- (off - mean(off)) / (stats::sd(off) + 1e-9)
    odd <- colnames(cm)[z < -2.5 & off < stats::median(off) - 0.03]
    if (length(odd)) note("warning", sprintf("Sample(s) that look unlike the rest: %s", paste(utils::head(odd, 5), collapse = ", ")), paste(utils::head(odd, 5), collapse = ", "),
                          "Check the correlation picture and the PCA: if it sits away from its replicates, consider leaving it out (Analysis settings: samples to leave out).")
    utils::write.csv(data.frame(sample = colnames(cm), mean_correlation_with_others = round(off, 4)), file.path(out, "sample_correlation_summary.csv"), row.names = FALSE); files <- c(files, "sample_correlation_summary.csv")
  } else note("info", "The sample correlation picture was not drawn", "samples", "It needs at least 3 samples and the pheatmap package.")

  nt <- if (length(notes)) do.call(rbind, notes) else data.frame(level = character(), issue = character(), where = character(), what_to_do = character())
  if (!nrow(nt)) nt <- data.frame(level = "ok", issue = "No problems found in the design or the samples", where = "", what_to_do = "", stringsAsFactors = FALSE)
  utils::write.csv(nt, file.path(out, "data_overview_notes.csv"), row.names = FALSE); files <- c(files, "data_overview_notes.csv")
  nw <- sum(nt$level == "warning")
  list(status = "ok", files = files, message = sprintf("%d samples, %s genes: %s%s.", nrow(S), format(nrow(M), big.mark = ","),
                                                        if (nw) sprintf("%d thing%s to look at (see data_overview_notes.csv)", nw, if (nw == 1) "" else "s") else "no problems found",
                                                        if (any(nt$level == "info")) sprintf("; %d note%s", sum(nt$level == "info"), if (sum(nt$level == "info") == 1) "" else "s") else ""))
}
