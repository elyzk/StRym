# analysis/shannon_entropy.R -- Shannon entropy genes from Tux (Cuffdiff) results (from the baseline shannon_simple.R, generalized).
#
# Standard Shannon entropy of each gene's mean-FPKM time profile (Schug et al. 2005): p_t = FPKM_t / sum(FPKM),
# H = -sum(p_t log2 p_t), H_norm = H / log2(#time points); 1 = flat, lower = changes over time. It does not say WHICH time
# points change, so those come from Cuffdiff (treated vs control fold change + q-value).
#   1. drop genes whose CONTROL profile changes over time (low control entropy) unless the treatment effect is much bigger;
#   2. keep genes whose TREATED profile changes over time (low entropy);
#   3. call the time points where treated differs from control (|log2 FC| and q).
# Input  : Tux results of this project (needs raw reads -> Quantify) and >= 2 time points; treated/control pairs are found from
#          the condition names (or set in Analysis settings). The first time point is a baseline unless switched off.
# Output : shannon_<treated>.csv per treated condition, shannon_main_genes.csv, shannon_comparisons_used.csv, a summary picture.
strym_step_shannon_entropy <- function(ctx) {
  p <- strym_normalize_params(ctx$design$params)
  data <- strym_drop_samples(strym_load_expression(ctx), p)
  data <- strym_ensure_diff(data, ctx)             # q-values: Tux / Cuffdiff files, or calculated with DESeq2 from raw counts
  if (is.null(data$diff)) stop("Shannon entropy needs q-values for treated vs control: start from raw reads (Tux), give raw counts with replicates, or give your Cuffdiff results (gene_exp.diff files) in 'Differential-expression results'.", call. = FALSE)
  tc <- strym_timecourse(data, p, ctx$design$timepoint_names, log = ctx$log)
  run <- strym_run_shannon(tc, p)
  files <- strym_write_calls(ctx$out_dir, "shannon", run, tc, "Shannon entropy")
  for (i in seq_len(nrow(tc$pairs))) {
    r <- run$results[[.pair_keys(tc$pairs)[i]]]
    ctx$log(sprintf("%s vs %s: %d control-variable genes excluded; %d genes called; main set (>= %d time points): %d",
                    tc$pairs$treated[i], tc$pairs$control[i], run$excluded[i], nrow(r), p$min_sig_tp, sum(r$main_set)))
  }
  n_main <- if (is.null(run$main)) 0L else length(unique(run$main$gene))
  list(status = "ok", files = files,
       message = sprintf("%d treated-vs-control comparison%s; %d genes in the main set (changed at >= %d time points).", nrow(tc$pairs), if (nrow(tc$pairs) == 1) "" else "s", n_main, p$min_sig_tp))
}
