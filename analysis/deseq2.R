# analysis/deseq2.R -- differential expression with DESeq2 (structure of the baseline "easy_deseq2" script, generalized).
#
# easy_deseq2 (github.com/zhangdongqin/easy_deseq2): count matrix + sample table -> DESeqDataSetFromMatrix -> DESeq (Wald test) ->
# normalized counts, results tables, up/down gene lists, volcano plot, heatmap of the DE genes. Changes made here:
#   * the comparison is explicit: every treated condition against its control (found from the condition names or set in Analysis
#     settings, as for the entropy steps), at every time point (design ~ group, group = condition_time point). The original calls
#     results() without a contrast, so its direction is whatever DESeq2 happens to pick;
#   * the thresholds are StRym's shared ones (adjusted p < q, |log2 FC| >= min_log2fc; defaults 0.05 and 1), not a mix of 0.1 / 0.01;
#   * the volcano title and axes follow the data (the original hard-codes "control vs treat" and fixed limits);
#   * the human GO / KEGG enrichment is left out (it only works for human gene symbols);
#   * genes with fewer than 10 counts in total are dropped before testing (saves time; DESeq2's own filtering would remove them anyway).
# Input  : RAW counts (whole numbers, e.g. the Tux counts.csv, or your counts table) + the sample sheet (condition, time_point).
#          DESeq2 needs replicates: at least one group with 2 or more samples.
# Output : deseq2_<comparison>.csv (every tested gene), deseq2_DEGs.csv (the significant genes, all comparisons), normalized_counts.csv,
#          size_factors.csv, deseq2_comparisons_used.csv, deseq2_volcano.png, deseq2_DEG_heatmap.png.
strym_step_deseq2 <- function(ctx) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("The DESeq2 package is not installed in StRym's R environment. Run setup_modules.sh (see the README) to add it.", call. = FALSE)
  p <- strym_normalize_params(ctx$design$params); d <- ctx$design; out <- ctx$out_dir; log <- ctx$log
  raw <- .counts_raw(ctx$data_path, ctx$sample_sheet)
  drop <- strym_split_list(p$drop_reps); S <- raw$samples[!raw$samples$sample %in% drop, , drop = FALSE]
  if (length(setdiff(drop, raw$samples$sample))) log(sprintf("Not a sample name, ignored: %s.", paste(setdiff(drop, raw$samples$sample), collapse = ", ")))
  m <- raw$m[, S$sample, drop = FALSE]
  if (any(m < 0)) stop("The counts table has negative values: DESeq2 needs raw counts.", call. = FALSE)
  if (any(abs(m - round(m)) > 1e-6)) stop("DESeq2 needs RAW counts (whole numbers), but this table has fractional values: it looks normalized (FPKM / TPM / CPM). Use the raw counts instead.", call. = FALSE)
  m <- round(m); storage.mode(m) <- "integer"
  tps <- strym_order_timepoints(S$time_point, d$timepoint_names); multi <- length(tps) > 1
  S$group <- if (multi) paste(S$condition, S$time_point, sep = "_") else S$condition
  if (!any(table(S$group) >= 2)) stop("DESeq2 needs replicates: every group (condition at a time point) has a single sample.", call. = FALSE)
  if (length(unique(S$condition)) < 2) stop("DESeq2 compares conditions, but the sample sheet has only one.", call. = FALSE)

  # which comparisons: treated vs control; with no control to be found and exactly two conditions: the 2nd vs the 1st
  pairs <- tryCatch(strym_pairs(S, p), error = function(e) {
    conds <- unique(S$condition)
    if (length(conds) != 2) stop(conditionMessage(e), call. = FALSE)
    log(sprintf("No control condition found: comparing %s with %s (positive log2 fold change = higher in %s).", conds[2], conds[1], conds[2]))
    data.frame(treated = conds[2], control = conds[1], stringsAsFactors = FALSE) })
  groups <- if (multi) tps else ""
  comps <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) do.call(rbind, lapply(groups, function(t) {
    g1 <- if (multi) paste(pairs$treated[i], t, sep = "_") else pairs$treated[i]; g2 <- if (multi) paste(pairs$control[i], t, sep = "_") else pairs$control[i]
    data.frame(treated = pairs$treated[i], control = pairs$control[i], time_point = t, g_treated = g1, g_control = g2, stringsAsFactors = FALSE) }))))
  lost <- comps[!(comps$g_treated %in% S$group & comps$g_control %in% S$group), , drop = FALSE]
  if (nrow(lost)) log(sprintf("Left out (no samples in one of the two groups): %s.", paste(paste0(lost$treated, " vs ", lost$control, if (multi) paste0(" at ", lost$time_point) else ""), collapse = "; ")))
  comps <- comps[comps$g_treated %in% S$group & comps$g_control %in% S$group, , drop = FALSE]
  if (!nrow(comps)) stop("None of the treated / control comparisons has samples in both groups.", call. = FALSE)
  comps$label <- gsub("[^A-Za-z0-9._-]+", "_", paste0(comps$treated, "_vs_", comps$control, if (multi) paste0("_", comps$time_point) else ""))

  # ---- DESeq2 ----
  suppressPackageStartupMessages(library(DESeq2))
  keep <- rowSums(m) >= 10
  if (sum(keep) < 5) stop("Fewer than 5 genes have 10 or more counts in total: nothing to test.", call. = FALSE)
  coldata <- data.frame(group = factor(S$group), row.names = S$sample)
  dds <- DESeqDataSetFromMatrix(countData = m[keep, , drop = FALSE], colData = coldata, design = ~ group)
  dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) {
    sf <- "ratio"
    if (grepl("zero", conditionMessage(e), fixed = TRUE)) {
      sf <- "poscounts"; log("Every gene has at least one zero count (a small gene panel or low depth): using size factors that allow zeros (poscounts).") }
    tryCatch(DESeq(dds, sfType = sf, quiet = TRUE), error = function(e2) {       # e.g. locfit's "out of vertex space" on a very small gene panel
      log("The dispersion curve could not be fitted (very few genes): using one common dispersion value instead.")
      DESeq(dds, sfType = sf, fitType = "mean", quiet = TRUE) }) })
  log(sprintf("DESeq2: %d of %d genes tested (>= 10 counts), %d samples, %d comparison%s.", sum(keep), nrow(m), ncol(m), nrow(comps), if (nrow(comps) == 1) "" else "s"))
  nc <- counts(dds, normalized = TRUE)
  utils::write.csv(data.frame(gene = rownames(nc), round(nc, 3), check.names = FALSE), file.path(out, "normalized_counts.csv"), row.names = FALSE)
  utils::write.csv(data.frame(sample = names(sizeFactors(dds)), size_factor = round(unname(sizeFactors(dds)), 4)), file.path(out, "size_factors.csv"), row.names = FALSE)
  files <- c("normalized_counts.csv", "size_factors.csv")

  alls <- list(); degs <- list()
  for (i in seq_len(nrow(comps))) {
    r <- as.data.frame(results(dds, contrast = c("group", comps$g_treated[i], comps$g_control[i]), alpha = p$q_cutoff))
    r$padj[is.na(r$padj)] <- 1; r$pvalue[is.na(r$pvalue)] <- 1                                 # untestable genes count as not significant (as in easy_deseq2)
    r$mean_treated <- rowMeans(nc[rownames(r), S$sample[S$group == comps$g_treated[i]], drop = FALSE])
    r$mean_control <- rowMeans(nc[rownames(r), S$sample[S$group == comps$g_control[i]], drop = FALSE])
    res <- data.frame(gene = rownames(r), r, check.names = FALSE, row.names = NULL); res <- res[order(res$pvalue), ]
    res$direction <- ifelse(res$padj < p$q_cutoff & res$log2FoldChange >= p$min_log2fc, "up", ifelse(res$padj < p$q_cutoff & res$log2FoldChange <= -p$min_log2fc, "down", "none"))
    f <- sprintf("deseq2_%s.csv", comps$label[i]); utils::write.csv(res, file.path(out, f), row.names = FALSE); files <- c(files, f)
    res$comparison <- paste0(comps$treated[i], "_vs_", comps$control[i]); res$time_point <- comps$time_point[i]; alls[[i]] <- res
    degs[[i]] <- res[res$direction != "none", , drop = FALSE]
    log(sprintf("%s: %d up, %d down (padj < %g, |log2FC| >= %g).", comps$label[i], sum(res$direction == "up"), sum(res$direction == "down"), p$q_cutoff, p$min_log2fc))
  }
  deg <- do.call(rbind, degs); if (is.null(deg)) deg <- data.frame()
  if (nrow(deg)) deg <- deg[, c("gene", "comparison", "time_point", "direction", setdiff(names(deg), c("gene", "comparison", "time_point", "direction")))]
  utils::write.csv(deg, file.path(out, "deseq2_DEGs.csv"), row.names = FALSE); files <- c(files, "deseq2_DEGs.csv")
  utils::write.csv(data.frame(comps[, c("treated", "control", "time_point", "g_treated", "g_control")], up = vapply(alls, function(x) sum(x$direction == "up"), 1L),
                              down = vapply(alls, function(x) sum(x$direction == "down"), 1L)), file.path(out, "deseq2_comparisons_used.csv"), row.names = FALSE)
  files <- c(files, "deseq2_comparisons_used.csv")

  # ---- volcano plot (axes follow the data) and heatmap of the DE genes ----
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    v <- do.call(rbind, alls); v$panel <- paste0(v$comparison, if (multi) paste0(" @ ", v$time_point) else "")
    v$panel <- factor(v$panel, levels = unique(v$panel)); v$y <- -log10(pmax(v$padj, 1e-300)); v$direction <- factor(v$direction, levels = c("up", "down", "none"))
    g <- ggplot2::ggplot(v, ggplot2::aes(log2FoldChange, y, colour = direction)) + ggplot2::geom_point(size = 0.8, alpha = 0.7) +
      ggplot2::geom_vline(xintercept = c(-p$min_log2fc, p$min_log2fc), linetype = "dashed", colour = "grey50") +
      ggplot2::geom_hline(yintercept = -log10(p$q_cutoff), linetype = "dashed", colour = "grey50") +
      ggplot2::scale_colour_manual(values = c(up = "#6F7D1C", down = "#4156A1", none = "#B2B2B2"), drop = FALSE) +
      ggplot2::facet_wrap(~panel) + ggplot2::labs(x = "log2 fold change (treated / control)", y = "-log10 adjusted p-value", colour = NULL) + ggplot2::theme_bw()
    ggplot2::ggsave(file.path(out, "deseq2_volcano.png"), g, width = 3 + 3 * min(4, nlevels(v$panel)), height = 3.2 + 2.6 * ceiling(nlevels(v$panel) / 4), dpi = 120, limitsize = FALSE)
    files <- c(files, "deseq2_volcano.png")
  }
  top <- unique(unlist(lapply(degs, function(x) utils::head(x$gene[order(x$padj)], 50))))
  if (length(top) >= 2 && requireNamespace("pheatmap", quietly = TRUE)) {
    top <- utils::head(top, 100); hm <- log2(nc[top, , drop = FALSE] + 1); hm <- hm[apply(hm, 1, sd) > 0, , drop = FALSE]
    if (nrow(hm) >= 2) {
      ann <- data.frame(group = S$group, row.names = S$sample)
      grDevices::png(file.path(out, "deseq2_DEG_heatmap.png"), width = 900, height = 300 + 9 * nrow(hm), res = 110)
      pheatmap::pheatmap(hm, scale = "row", annotation_col = ann, show_rownames = nrow(hm) <= 60, main = sprintf("Top DE genes (log2 normalized counts, row-scaled)"), fontsize = 8)
      grDevices::dev.off(); files <- c(files, "deseq2_DEG_heatmap.png")
    }
  }
  n_deg <- if (nrow(deg)) length(unique(deg$gene)) else 0L
  list(status = "ok", files = files,
       message = sprintf("%d comparison%s; %d differentially expressed gene%s (padj < %g, |log2 FC| >= %g).", nrow(comps), if (nrow(comps) == 1) "" else "s", n_deg, if (n_deg == 1) "" else "s", p$q_cutoff, p$min_log2fc))
}
