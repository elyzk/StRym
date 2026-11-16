# analysis/modified_shannon_entropy.R -- modified Shannon entropy (MSE / ROKU-style outliers) from Tux (Cuffdiff) results
# (from the baseline mse_simple.R, generalized).
#
# For every gene the mean-FPKM time profile is scanned for outlier time points (+1 high, -1 low, 0 normal) with an
# information-criterion model choice (up to s_max outliers). Then:
#   control : genes that change over time WITHOUT treatment (control outliers >= 1.5-fold) are excluded,
#             unless the treatment effect is >= rescue_ratio x that drift;
#   treated : MSE outliers that also differ from the control (|log2 FC| >= min_log2fc and Cuffdiff q < q_cutoff) are called;
#             a sustained response (MSE flags the few unaffected points) is called on the remaining time points.
# Input  : Tux results of this project (needs raw reads -> Quantify) and >= 3 time points; treated/control pairs are found from
#          the condition names (or set in Analysis settings). The first time point is a baseline unless switched off.
# Output : MSE_<treated>.csv per treated condition, MSE_main_genes.csv, MSE_comparisons_used.csv, a summary picture.
strym_step_modified_shannon_entropy <- function(ctx) {
  p <- strym_normalize_params(ctx$design$params)
  data <- strym_drop_samples(strym_load_expression(ctx), p)
  data <- strym_ensure_diff(data, ctx)             # q-values: Tux / Cuffdiff files, or calculated with DESeq2 from raw counts
  if (is.null(data$diff)) stop("Modified Shannon entropy needs q-values for treated vs control: start from raw reads (Tux), give raw counts with replicates, or give your Cuffdiff results (gene_exp.diff files) in 'Differential-expression results'.", call. = FALSE)
  tc <- strym_timecourse(data, p, ctx$design$timepoint_names, log = ctx$log)
  if (length(tc$tps) < 3) stop(sprintf("Modified Shannon entropy needs at least 3 time points to tell outliers from the rest; this experiment has %d.", length(tc$tps)), call. = FALSE)
  run <- strym_run_mse(tc, p)
  files <- strym_write_calls(ctx$out_dir, "MSE", run, tc, "Modified Shannon entropy", params = p)
  for (i in seq_len(nrow(tc$pairs))) {
    r <- run$results[[.pair_keys(tc$pairs)[i]]]
    ctx$log(sprintf("%s vs %s: %d control-variable genes excluded; %d genes called; main set (>= %d time points): %d",
                    tc$pairs$treated[i], tc$pairs$control[i], run$excluded[i], nrow(r), p$min_sig_tp, sum(r$main_set)))
  }
  n_main <- if (is.null(run$main)) 0L else length(unique(run$main$gene))
  list(status = "ok", files = files,
       message = sprintf("%d treated-vs-control comparison%s; %d genes in the main set (changed at >= %d time points).", nrow(tc$pairs), if (nrow(tc$pairs) == 1) "" else "s", n_main, p$min_sig_tp))
}
