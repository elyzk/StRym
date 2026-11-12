# analysis/derived_diff.R -- q-values for the time-course steps (Shannon entropy, modified Shannon entropy) from a RAW counts table.
#
# The time-course steps compare each treated condition with its control at every time point and need a q-value for every gene in every
# comparison. A Tux (Cuffdiff) run provides them. When there is no Tux run -- the reads were quantified elsewhere (long reads such as
# Oxford Nanopore / Plasmidsaurus, 3'-tag or any other protocol, featureCounts / HTSeq / Salmon tables) -- the same comparisons are made
# here with DESeq2 (one small model per treated-vs-control pair and time point), and the result has the same shape as Cuffdiff's
# gene_exp.diff, so everything downstream is unchanged.
#
# Needs: raw (whole-number) counts and at least 3 samples in each comparison (e.g. 2 control + 1 treated; 2 + 2 or more is better).
# A comparison that cannot be made (too few samples, a model that fails) gets q = 1 for every gene and is reported, never invented.
# The result is remembered in the project folder, so the entropy steps of one run share it.

# The step calls this after loading its data: when there are no Cuffdiff q-values (a counts table), calculate them. Safe to call always.
strym_ensure_diff <- function(data, ctx) {
  if (!is.null(data$diff) || is.null(data$counts)) return(data)
  d <- ctx$design; p <- strym_normalize_params(d$params)
  data$diff <- strym_derive_diff(data$counts, data$samples, p, d$timepoint_names,
                                 cache_file = file.path(ctx$project_dir, "02_differential_expression", ".derived_diff_cache.rds"),
                                 threads = max(1L, min(as.integer(d$threads), 4L)), log = if (is.function(ctx$log)) ctx$log else message)
  data$diff_source <- "DESeq2"
  data
}

strym_derive_possible <- function(d) {                       # "" = can be done for this design, otherwise why not (used by the plan)
  if (!identical(d$data_kind, "counts") && !(identical(d$data_kind, "fastq") && identical(d$read_platform, "long"))) return("Only for an expression table or long reads.")
  if (identical(d$expr_unit, "normalized"))
    return("The table holds normalized values (FPKM / TPM), from which q-values cannot be calculated: add your Cuffdiff results, or give raw counts.")
  if (!requireNamespace("DESeq2", quietly = TRUE)) return("The DESeq2 package is not installed.")
  sh <- if (nzchar(d$sample_sheet) && file.exists(d$sample_sheet)) strym_read_sample_sheet(d$sample_sheet) else NULL
  if (is.data.frame(sh)) {
    cc <- .find_col(sh, "^condition"); tc <- .find_col(sh, "^time")
    if (!is.na(cc)) {
      g <- paste(sh[[cc]], if (!is.na(tc)) sh[[tc]] else "")
      if (!any(table(g) >= 2)) return("No condition has replicates (two or more samples per condition and time point), which the significance tests need.")
    }
  }
  ""
}

# m = genes x samples RAW counts; samples = data.frame(sample, condition, time_point, replicate, group)
# Returns data.frame(test_id, sample_1 (control group), sample_2 (treated group), status, q_value, log2fc), plus attr "skipped" (text).
strym_derive_diff <- function(m, samples, params, tp_names = character(), cache_file = "", threads = 1, log = message) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("The DESeq2 package is needed to calculate q-values from counts.", call. = FALSE)
  S <- samples; S <- S[S$sample %in% colnames(m), , drop = FALSE]
  pairs <- strym_pairs(S, params); tps <- strym_order_timepoints(S$time_point, tp_names)
  if (!length(tps)) stop("The sample sheet has no time points: the entropy steps compare treated with control at several time points (add a time_point column).", call. = FALSE)
  # the key is a checksum of EVERYTHING the result depends on (every count, gene id, sample, group, pair and time point): swapping two columns or
  # renaming genes must never reuse another dataset's q-values
  kf <- tempfile(fileext = ".rds"); on.exit(unlink(kf), add = TRUE)
  saveRDS(list(m, S[, c("sample", "condition", "time_point", "group")], pairs, tps), kf, compress = FALSE)
  key <- unname(tools::md5sum(kf))
  if (nzchar(cache_file) && file.exists(cache_file)) {
    old <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.list(old) && identical(old$key, key)) { log("q-values: reusing the DESeq2 comparisons calculated earlier in this project."); return(old$diff) }
  }
  if (any(abs(m - round(m)) > 1e-6, na.rm = TRUE))
    stop("Calculating q-values needs RAW counts (whole numbers), but this table has decimals: it looks normalized. Give the raw counts, or your Cuffdiff results.", call. = FALSE)
  jobs <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) data.frame(i = i, tp = tps, stringsAsFactors = FALSE)))
  one <- function(k) {
    i <- jobs$i[k]; tp <- jobs$tp[k]; ct <- pairs$control[i]; tr <- pairs$treated[i]
    s_c <- S[S$condition == ct & S$time_point == tp, , drop = FALSE]; s_t <- S[S$condition == tr & S$time_point == tp, , drop = FALSE]
    g1 <- unique(s_c$group)[1]; g2 <- unique(s_t$group)[1]
    if (!nrow(s_c) || !nrow(s_t)) return(list(why = sprintf("%s vs %s at %s: one of them has no samples at this time point", tr, ct, tp)))
    notest <- function(why) list(why = why, diff = data.frame(test_id = rownames(m), sample_1 = g1, sample_2 = g2, status = "NOTEST", q_value = 1, log2fc = 0, stringsAsFactors = FALSE))
    if (nrow(s_c) + nrow(s_t) < 3) return(notest(sprintf("%s vs %s at %s: only %d samples in total (3 are needed)", tr, ct, tp, nrow(s_c) + nrow(s_t))))
    smp <- c(s_c$sample, s_t$sample); cnt <- round(m[, smp, drop = FALSE]); cnt[is.na(cnt)] <- 0
    keep <- rowSums(cnt) > 0; if (sum(keep) < 2) return(notest(sprintf("%s vs %s at %s: no gene has any counts", tr, ct, tp)))
    cond <- factor(c(rep("control", nrow(s_c)), rep("treated", nrow(s_t))), levels = c("control", "treated"))
    dds <- DESeq2::DESeqDataSetFromMatrix(cnt[keep, , drop = FALSE], data.frame(cond = cond, row.names = smp), ~ cond)
    fit <- NULL
    for (ft in c("parametric", "local", "mean")) {                       # the usual fallbacks when the dispersion trend cannot be fitted
      fit <- tryCatch({ d <- DESeq2::estimateSizeFactors(dds, type = "poscounts"); d <- DESeq2::estimateDispersions(d, fitType = ft, quiet = TRUE)
                        DESeq2::nbinomWaldTest(d, quiet = TRUE) }, error = function(e) NULL, warning = function(w) NULL)
      if (!is.null(fit)) break
    }
    if (is.null(fit)) {                                                    # warnings above are not fatal: retry ignoring them
      for (ft in c("parametric", "local", "mean")) {
        fit <- tryCatch(suppressWarnings({ d <- DESeq2::estimateSizeFactors(dds, type = "poscounts"); d <- DESeq2::estimateDispersions(d, fitType = ft, quiet = TRUE)
                                           DESeq2::nbinomWaldTest(d, quiet = TRUE) }), error = function(e) NULL)
        if (!is.null(fit)) break
      }
    }
    if (is.null(fit)) return(notest(sprintf("%s vs %s at %s: the DESeq2 model could not be fitted", tr, ct, tp)))
    r <- suppressWarnings(DESeq2::results(fit, contrast = c("cond", "treated", "control"), independentFiltering = FALSE))
    q <- r$padj; ok <- !is.na(q); q[!ok] <- 1
    list(diff = data.frame(test_id = rownames(r), sample_1 = g1, sample_2 = g2, status = ifelse(ok, "OK", "NOTEST"), q_value = q,
                           log2fc = r$log2FoldChange, stringsAsFactors = FALSE))
  }
  threads <- max(1L, as.integer(threads))
  big <- as.numeric(nrow(m)) * ncol(m) > 5e6                                    # a very large table: every parallel worker holds its own copy, so use fewer workers
  res <- if (threads > 1 && .Platform$OS.type != "windows") parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = min(if (big) 2L else threads, nrow(jobs))) else lapply(seq_len(nrow(jobs)), one)
  bad <- vapply(res, function(x) !is.list(x), TRUE)
  if (any(bad)) {                                                               # a worker that was stopped (usually for lack of memory): redo those comparisons one at a time
    if (is.function(log)) log(sprintf("%d comparison(s) could not run in parallel (probably not enough memory); repeating them one at a time.", sum(bad)))
    res[bad] <- lapply(which(bad), one); bad <- vapply(res, function(x) !is.list(x), TRUE) }
  if (any(bad)) stop("The DESeq2 comparisons failed unexpectedly: ", paste(unique(vapply(res[bad], function(x) as.character(x)[1], "")), collapse = "; "), call. = FALSE)
  why <- unlist(lapply(res, `[[`, "why")); parts <- lapply(res, `[[`, "diff"); parts <- parts[!vapply(parts, is.null, TRUE)]
  if (!length(parts) || !any(vapply(parts, function(d) any(d$status == "OK"), TRUE)))
    stop("No treated-vs-control comparison could be made: ", paste(why, collapse = "; "), ".", call. = FALSE)
  de <- do.call(rbind, parts); attr(de, "skipped") <- why
  log(sprintf("q-values: %d of %d treated-vs-control comparisons calculated with DESeq2 (%s genes tested in total).", length(parts), nrow(jobs), format(nrow(de), big.mark = ",")))
  if (length(why)) log(sprintf("Comparisons that could not be made (their q-values are 1): %s.", paste(why, collapse = "; ")))
  if (nzchar(cache_file)) tryCatch({ dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE); saveRDS(list(key = key, diff = de), cache_file) }, error = function(e) NULL)
  de
}
