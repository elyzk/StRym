# Shannon entropy (on mean FPKM) heat-response genes from Cuffdiff output
#
# Standard Shannon entropy of each gene's mean-FPKM time profile (Schug et al. 2005):
#   p_t = FPKM_t / sum(FPKM),  H = -sum(p_t * log2 p_t),  H_norm = H / log2(#timepoints)
#   H_norm = 1: same level at every timepoint; lower = expression changes over time.
# Unlike modified Shannon entropy (mse_simple.R), it does not say WHICH timepoints change,
# so the timepoints are taken from Cuffdiff (heat vs control fold change + q).
#
# Steps: (1) drop genes whose CONTROL profile changes over time (low control entropy),
#        unless the heat effect is much bigger than that control change;
#        (2) keep genes whose HEAT profile changes over time (low heat entropy);
#        (3) call timepoints where heat differs from control (fold change + q).
#
# Input : same folder layout as mse_simple.R
#           genes.read_group_tracking_<genotype>_T<n>  and  gene_exp.diff_<genotype>_T<n>[.diff]
#         First timepoint = baseline (before heat): reported, never called.
# Output: <out_dir>/shannon_<genotype>.csv, <out_dir>/shannon_main_genes.csv
#
# Run   : Rscript --vanilla shannon_simple.R [data_dir] [out_dir]

#------------------------------------------------------------------------------#
# Settings
#------------------------------------------------------------------------------#
args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "/Users/jlreed4/Downloads/genes_read_group_tracking_files"
out_dir  <- if (length(args) >= 2) args[2] else "shannon_simple_results"

heat_pattern  <- "heat"   # condition names containing this = heat; everything else = control
min_fpkm      <- 2        # gene must reach this mean FPKM in some sample
min_log2fc    <- 1        # heat timepoint needs |log2 (heat+1)/(control+1)| >= this ...
q_cutoff      <- 0.05     # ... and Cuffdiff q < this
min_sig_tp    <- 2        # main set: significant at >= this many timepoints
h_control     <- 0.99     # control H_norm below this = changes over time in control -> excluded, those are genes that fluctuate naturally
h_heat        <- 0.99     # heat H_norm must be below this = changes over time under heat
rescue_ratio  <- 2        # keep a control-variable gene if heat effect >= this x control drift
drop_reps     <- character(0)   # replicates to leave out, e.g. c("Pi_Heated_T1:0"), outliers from pca potentially


# Read data -> mean FPKM matrices (genes x timepoints), log2FC and q-values

dir.create(out_dir, showWarnings = FALSE)
files <- list.files(data_dir, pattern = "^genes\\.read_group_tracking_.+_T[0-9]+", full.names = TRUE)
info  <- data.frame(file = files,
                    genotype = sub(".*tracking_(.+)_T[0-9]+.*", "\\1", basename(files)),
                    t = as.integer(sub(".*_T([0-9]+).*", "\\1", basename(files))))
genotypes  <- unique(info$genotype)
tp_numbers <- sort(unique(info$t))
timepoints <- paste0("T", seq_along(tp_numbers) - 1)        # first file timepoint = T0 (baseline)
cat("Found", nrow(info), "tracking files; genotypes:", paste(genotypes, collapse = ", "),
    "; timepoints:", paste(timepoints, collapse = " "), "\n")

reps <- do.call(rbind, lapply(seq_len(nrow(info)), function(i) {
  d <- read.delim(info$file[i], stringsAsFactors = FALSE)
  d$FPKM[d$status %in% c("FAIL", "HIDATA")] <- NA
  heat <- grepl(heat_pattern, d$condition, ignore.case = TRUE)
  data.frame(gene = d$tracking_id, FPKM = d$FPKM, heat = heat, genotype = info$genotype[i],
             tp = timepoints[match(info$t[i], tp_numbers)],
             rep_id = sprintf("%s_%s_T%d:%s", info$genotype[i], ifelse(heat, "Heated", "Control"),
                              info$t[i], d$replicate))
}))
reps <- reps[!reps$rep_id %in% drop_reps, ]
all_genes <- sort(unique(reps$gene))

mean_matrix <- function(g, is_heat) {
  x <- reps[reps$genotype == g & reps$heat == is_heat, ]
  m <- tapply(x$FPKM, list(x$gene, x$tp), mean, na.rm = TRUE)
  m <- m[all_genes, timepoints]; m[is.na(m)] <- 0
  m
}
heat <- sapply(genotypes, mean_matrix, is_heat = TRUE,  simplify = FALSE)
ctrl <- sapply(genotypes, mean_matrix, is_heat = FALSE, simplify = FALSE)
lfc  <- sapply(genotypes, function(g) log2((heat[[g]] + 1) / (ctrl[[g]] + 1)), simplify = FALSE)

qval <- sapply(genotypes, function(g) {
  m <- matrix(1, length(all_genes), length(timepoints), dimnames = list(all_genes, timepoints))
  for (k in seq_along(tp_numbers)) {
    f <- list.files(data_dir, pattern = sprintf("^gene_exp\\.diff_%s_T%d(\\.|$)", g, tp_numbers[k]),
                    full.names = TRUE)[1]
    d <- read.delim(f, stringsAsFactors = FALSE)
    q <- ifelse(d$status == "OK", d$q_value, 1)
    keep <- d$test_id %in% all_genes
    m[d$test_id[keep], k] <- q[keep]
  }
  m
}, simplify = FALSE)

universe <- all_genes[apply(do.call(cbind, c(heat, ctrl)), 1, max) >= min_fpkm]
cat(length(universe), "expressed genes (mean FPKM >=", min_fpkm, "somewhere)\n")

# Shannon entropy of each gene's profile (1 = flat, lower = changes over time)
entropy <- function(mat) {
  p <- mat / rowSums(mat)
  H <- -rowSums(ifelse(p > 0, p * log2(p), 0))
  ifelse(rowSums(mat) > 0, H / log2(ncol(mat)), 1)
}

# Control: genes with low control entropy are excluded (in both genotypes)

H_ctrl <- sapply(genotypes, function(g) entropy(ctrl[[g]][universe, ]))
drift  <- sapply(genotypes, function(g) {                      # largest control shift, log2
  x <- log2(ctrl[[g]][universe, ] + 1); apply(abs(x - apply(x, 1, median)), 1, max) })
drift  <- apply(ifelse(H_ctrl < h_control, drift, 0), 1, max)  # 0 = stable in both genotypes


# Heat: low heat entropy + timepoints where heat differs from control

results <- list()
for (g in genotypes) {
  sig <- abs(lfc[[g]]) >= min_log2fc & qval[[g]] < q_cutoff
  sig[, 1] <- FALSE                                          # baseline is never a call
  heat_effect <- apply(ifelse(sig, abs(lfc[[g]]), 0), 1, max)
  excluded <- universe[drift[universe] > 0 & heat_effect[universe] < rescue_ratio * drift[universe]]
  genes <- setdiff(universe, excluded)

  H_heat <- entropy(heat[[g]][genes, ])
  genes  <- genes[H_heat < h_heat]
  calls  <- sign(lfc[[g]][genes, , drop = FALSE]) * sig[genes, , drop = FALSE]
  calls  <- calls[rowSums(calls != 0) > 0, , drop = FALSE]

  res <- data.frame(gene = rownames(calls), genotype = g,
                    H_heat = round(H_heat[rownames(calls)], 4),
                    H_control = round(H_ctrl[rownames(calls), g], 4),
                    n_sig_timepoints = rowSums(calls != 0),
                    main_set = rowSums(calls != 0) >= min_sig_tp,
                    direction = ifelse(rowSums(calls == 1) > 0 & rowSums(calls == -1) > 0, "mixed",
                                       ifelse(rowSums(calls == 1) > 0, "up", "down")),
                    setNames(as.data.frame(calls), paste0("call_", timepoints)),
                    setNames(as.data.frame(round(lfc[[g]][rownames(calls), ], 3)), paste0("log2FC_", timepoints)),
                    setNames(as.data.frame(signif(qval[[g]][rownames(calls), ], 3)), paste0("q_", timepoints)))
  res <- res[order(!res$main_set, -res$n_sig_timepoints), ]
  write.csv(res, file.path(out_dir, sprintf("shannon_%s.csv", g)), row.names = FALSE)
  results[[g]] <- res
  cat(sprintf("%s: %d control-variable genes excluded; %d genes called; main set (>= %d timepoints): %d\n",
              g, length(excluded), nrow(res), min_sig_tp, sum(res$main_set)))
}
main <- do.call(rbind, lapply(results, function(r) r[r$main_set, ]))
write.csv(main, file.path(out_dir, "shannon_main_genes.csv"), row.names = FALSE)
cat("Main set, either genotype:", length(unique(main$gene)), "genes ->", out_dir, "\n")
