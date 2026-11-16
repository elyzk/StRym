# analysis/pca.R -- PCA of the samples (from the baseline pca_plot_script.R, generalized).
#
# Input  : Tux results (Cuffdiff genes.read_group_tracking, FPKM per replicate) or, for a counts table, counts per million.
#          Samples to leave out (e.g. outliers seen in an earlier PCA): Analysis settings -> "Samples to leave out".
# Level  : Analysis settings -> "PCA level" (pca_level):
#            replicate  every sample (replicate) is one point: shows how the replicates group and which one is an outlier   [default]
#            average    the replicates of each condition x time point are averaged first: one point per condition and time point
#            both       the two plots
# Method : log2(value + 1); genes with mean >= pca_min_mean; the pca_genes (500) most variable genes; prcomp on samples (or averages).
# Output : replicate level: pca_scores.csv, pca_variance.csv, pca_plot.png (PC1 vs PC2; label = sample name, to spot outliers)
#          average level  : pca_averaged_scores.csv, pca_averaged_variance.csv, pca_averaged_plot.png (label = condition_time point)

# One PCA. `m` is genes x points (already in the working unit), `info` has one row per column of m
# (sample, condition, time_point, replicate or n_replicates). `tag` is "" or "averaged_" (file prefix); returns the files written.
.pca_one <- function(m, info, out_dir, tag, tps, data_unit, p, point_label, title_what) {
  m <- log2(m[rowMeans(m) >= p$pca_min_mean, , drop = FALSE] + 1)                  # expressed genes, compressed scale
  if (nrow(m) < 3) stop(sprintf("Fewer than 3 genes reach a mean %s of %g.", data_unit, p$pca_min_mean), call. = FALSE)
  top <- order(apply(m, 1, var), decreasing = TRUE)[seq_len(min(p$pca_genes, nrow(m)))]    # the most variable genes inform the PCA
  pca <- stats::prcomp(t(m[top, , drop = FALSE]))                                  # samples (or averages) as rows
  pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  npc <- min(10, ncol(pca$x))
  scores <- data.frame(info[match(rownames(pca$x), info$sample), , drop = FALSE], pca$x[, seq_len(npc), drop = FALSE], row.names = NULL)
  f_scores <- sprintf("pca_%sscores.csv", tag); f_var <- sprintf("pca_%svariance.csv", tag); f_plot <- sprintf("pca_%splot.png", tag)
  utils::write.csv(scores, file.path(out_dir, f_scores), row.names = FALSE)
  utils::write.csv(data.frame(PC = paste0("PC", seq_along(pve)), percent_variance = round(pve, 2), cumulative_percent = round(cumsum(pve), 2)),
                   file.path(out_dir, f_var), row.names = FALSE)

  multi <- length(tps) > 1
  scores$colour_by <- if (multi) factor(scores$time_point, levels = tps) else factor(scores$condition)
  scores$condition <- factor(scores$condition, levels = unique(info$condition))
  pcx <- if (ncol(pca$x) >= 2) "PC2" else "PC1"
  gg <- ggplot2::ggplot(scores, ggplot2::aes(PC1, .data[[pcx]])) +
    ggplot2::geom_point(ggplot2::aes(colour = colour_by, shape = condition), size = 3, stroke = 1.2) +
    ggplot2::scale_shape_manual(values = rep(c(16, 17, 15, 18, 1, 2, 0, 5), length.out = nlevels(scores$condition))) +
    ggplot2::scale_colour_viridis_d(option = "turbo") +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("%s (%.1f%%)", pcx, pve[min(2, length(pve))]),
                  colour = if (multi) "Time point" else "Condition", shape = "Condition",
                  title = sprintf("PCA of %d %s (%d most variable genes, log2(%s + 1))", nrow(scores), title_what, length(top), data_unit)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +                    # room for the sample labels at the edges (they were cut off)
    ggplot2::theme_bw()
  if (nrow(scores) <= 40) gg <- gg + ggplot2::geom_text(ggplot2::aes(label = .data[[point_label]]), vjust = -1, size = 2.6, check_overlap = TRUE, show.legend = FALSE)
  ggplot2::ggsave(file.path(out_dir, f_plot), gg, width = 7, height = 5, dpi = 300)
  list(files = c(f_scores, f_var, f_plot), pve = pve, n = nrow(scores), genes = length(top))
}

strym_step_pca <- function(ctx) {
  p <- strym_normalize_params(ctx$design$params)
  data <- strym_drop_samples(strym_load_expression(ctx), p)
  S <- data$samples
  if (nrow(S) < 3) stop(sprintf("PCA needs at least 3 samples; found %d.", nrow(S)), call. = FALSE)
  m <- strym_wide(data, na_bad = FALSE)[, S$sample, drop = FALSE]; m[is.na(m)] <- 0
  tps <- strym_order_timepoints(S$time_point, ctx$design$timepoint_names)
  level <- p$pca_level
  files <- character(); msg <- character()

  if (level %in% c("replicate", "both")) {
    r <- .pca_one(m, S[, c("sample", "condition", "time_point", "replicate")], ctx$out_dir, "", tps, data$unit, p, "sample", "samples")
    files <- c(files, r$files)
    ctx$log(sprintf("PCA: %d samples, %d genes; PC1 %.1f%%, PC2 %.1f%%", r$n, r$genes, r$pve[1], r$pve[min(2, length(r$pve))]))
    msg <- c(msg, sprintf("%d samples; PC1 explains %.1f%% and PC2 %.1f%% of the variance (%d most variable genes, %s).", r$n, r$pve[1], r$pve[min(2, length(r$pve))], r$genes, data$unit))
  }
  if (level %in% c("average", "both")) {
    grp <- paste(S$condition, S$time_point, sep = "_")                             # one point per condition x time point
    ug <- unique(grp)
    if (length(ug) < 3) {
      if (level == "average") stop(sprintf("The averaged PCA needs at least 3 condition x time point groups; this experiment has %d. Use the replicate level instead.", length(ug)), call. = FALSE)
      ctx$log(sprintf("Averaged PCA skipped: only %d condition x time point groups (it needs 3).", length(ug)))
      msg <- c(msg, sprintf("The averaged PCA was skipped: only %d condition x time point groups (it needs 3).", length(ug)))
    } else {
      avg <- do.call(cbind, lapply(ug, function(g) rowMeans(m[, S$sample[grp == g], drop = FALSE])))   # mean of the replicates (in the working unit, before the log)
      colnames(avg) <- ug
      first <- match(ug, grp)
      info <- data.frame(sample = ug, condition = S$condition[first], time_point = S$time_point[first],
                         n_replicates = as.integer(table(grp)[ug]), stringsAsFactors = FALSE)
      info$point <- ug
      r <- .pca_one(avg, info, ctx$out_dir, "averaged_", tps, data$unit, p, "point", "condition x time point averages")
      files <- c(files, r$files)
      ctx$log(sprintf("PCA (averaged): %d groups, %d genes; PC1 %.1f%%, PC2 %.1f%%", r$n, r$genes, r$pve[1], r$pve[min(2, length(r$pve))]))
      msg <- c(msg, sprintf("Averaged over replicates: %d condition x time point groups; PC1 %.1f%% and PC2 %.1f%%.", r$n, r$pve[1], r$pve[min(2, length(r$pve))]))
    }
  }
  list(status = "ok", files = files, message = paste(msg, collapse = " "))
}
