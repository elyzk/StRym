# analysis/pca.R -- replicate-level PCA (from the baseline pca_plot_script.R, generalized).
#
# Input  : Tux results (Cuffdiff genes.read_group_tracking, FPKM per replicate) or, for a counts table, counts per million.
#          Samples to leave out (e.g. outliers seen in an earlier PCA): Analysis settings -> "Samples to leave out".
# Method : log2(value + 1); genes with mean >= pca_min_mean; the pca_genes (500) most variable genes; prcomp on samples.
# Output : pca_scores.csv, pca_variance.csv, pca_plot.png (PC1 vs PC2; label = sample name, to spot outliers).
strym_step_pca <- function(ctx) {
  p <- strym_normalize_params(ctx$design$params)
  data <- strym_drop_samples(strym_load_expression(ctx), p)
  S <- data$samples
  if (nrow(S) < 3) stop(sprintf("PCA needs at least 3 samples; found %d.", nrow(S)), call. = FALSE)
  m <- strym_wide(data, na_bad = FALSE)[, S$sample, drop = FALSE]; m[is.na(m)] <- 0
  m <- log2(m[rowMeans(m) >= p$pca_min_mean, , drop = FALSE] + 1)                 # expressed genes, compressed scale
  if (nrow(m) < 3) stop(sprintf("Fewer than 3 genes reach a mean %s of %g.", data$unit, p$pca_min_mean), call. = FALSE)
  top <- order(apply(m, 1, var), decreasing = TRUE)[seq_len(min(p$pca_genes, nrow(m)))]    # the most variable genes inform the PCA
  pca <- stats::prcomp(t(m[top, , drop = FALSE]))                                  # samples as rows
  pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  npc <- min(10, ncol(pca$x))
  scores <- data.frame(S[match(rownames(pca$x), S$sample), c("sample", "condition", "time_point", "replicate")], pca$x[, seq_len(npc), drop = FALSE], row.names = NULL)
  utils::write.csv(scores, file.path(ctx$out_dir, "pca_scores.csv"), row.names = FALSE)
  utils::write.csv(data.frame(PC = paste0("PC", seq_along(pve)), percent_variance = round(pve, 2), cumulative_percent = round(cumsum(pve), 2)),
                   file.path(ctx$out_dir, "pca_variance.csv"), row.names = FALSE)

  tps <- strym_order_timepoints(S$time_point, ctx$design$timepoint_names)
  multi <- length(tps) > 1
  scores$colour_by <- if (multi) factor(scores$time_point, levels = tps) else factor(scores$condition)
  scores$condition <- factor(scores$condition, levels = unique(S$condition))
  pcx <- if (ncol(pca$x) >= 2) "PC2" else "PC1"
  gg <- ggplot2::ggplot(scores, ggplot2::aes(PC1, .data[[pcx]])) +
    ggplot2::geom_point(ggplot2::aes(colour = colour_by, shape = condition), size = 3, stroke = 1.2) +
    ggplot2::scale_shape_manual(values = rep(c(16, 17, 15, 18, 1, 2, 0, 5), length.out = nlevels(scores$condition))) +
    ggplot2::scale_colour_viridis_d(option = "turbo") +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("%s (%.1f%%)", pcx, pve[min(2, length(pve))]),
                  colour = if (multi) "Time point" else "Condition", shape = "Condition",
                  title = sprintf("PCA of %d samples (%d most variable genes, log2(%s + 1))", nrow(scores), length(top), data$unit)) +
    ggplot2::theme_bw()
  if (nrow(scores) <= 40) gg <- gg + ggplot2::geom_text(ggplot2::aes(label = sample), vjust = -1, size = 2.6, check_overlap = TRUE, show.legend = FALSE)
  ggplot2::ggsave(file.path(ctx$out_dir, "pca_plot.png"), gg, width = 7, height = 5, dpi = 300)
  ctx$log(sprintf("PCA: %d samples, %d genes; PC1 %.1f%%, PC2 %.1f%%", nrow(scores), length(top), pve[1], pve[min(2, length(pve))]))
  list(status = "ok", files = c("pca_scores.csv", "pca_variance.csv", "pca_plot.png"),
       message = sprintf("%d samples; PC1 explains %.1f%% and PC2 %.1f%% of the variance (%d most variable genes, %s).", nrow(scores), pve[1], pve[min(2, length(pve))], length(top), data$unit))
}
