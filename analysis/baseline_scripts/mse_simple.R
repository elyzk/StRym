# Modified Shannon entropy (MSE / ROKU-style) heat-response genes from Cuffdiff output
#
# Input : a folder with, per genotype x timepoint,
#           genes.read_group_tracking_<genotype>_T<n>   (replicate FPKM)
#           gene_exp.diff_<genotype>_T<n>[.diff]        (heat vs control q-values)
#         Each file holds one heat and one control condition. The first timepoint
#         is the baseline (before heat): reported, never called as a heat response.
# Output: <out_dir>/MSE_<genotype>.csv  per-gene timepoint calls (+1 up, -1 down)
#         <out_dir>/MSE_main_genes.csv  all main-set genes, both genotypes
#
# Run   : Rscript --vanilla mse_simple.R [data_dir] [out_dir]

#------------------------------------------------------------------------------#
# Settings
#------------------------------------------------------------------------------#
args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "/Users/jlreed4/Downloads/genes_read_group_tracking_files"
out_dir  <- if (length(args) >= 2) args[2] else "MSE_simple_results"

heat_pattern  <- "heat"   # condition names containing this = heat; everything else = control
min_fpkm      <- 2        # gene must reach this mean FPKM in some sample
min_log2fc    <- 1        # heat timepoint needs |log2 (heat+1)/(control+1)| >= this ...
q_cutoff      <- 0.05     # ... and Cuffdiff q < this
min_sig_tp    <- 2        # main set: significant at >= this many timepoints
ctrl_log2fc   <- 0.585    # control outlier must differ >= 1.5-fold from the rest
rescue_ratio  <- 2        # keep a control-variable gene if heat effect >= this x control drift
s_max         <- 3        # max outlier timepoints per gene
drop_reps     <- character(0)   # replicates to leave out, e.g. c("Pi_Heated_T1:0")

#------------------------------------------------------------------------------#
# 1. Read data -> mean FPKM matrices (genes x timepoints), log2FC and q-values
#------------------------------------------------------------------------------#
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

# genes expressed in at least one sample (shared by all MSE runs)
universe <- all_genes[apply(do.call(cbind, c(heat, ctrl)), 1, max) >= min_fpkm]
cat(length(universe), "expressed genes (mean FPKM >=", min_fpkm, "somewhere)\n")


# MSE: per gene, which timepoints are outliers (+1 high, -1 low, 0 normal)

tukey_biweight <- function(x, c = 5, eps = 1e-4) {
  u <- (x - median(x)) / (c * median(abs(x - median(x))) + eps)
  w <- ifelse(abs(u) <= 1, (1 - u^2)^2, 0)
  sum(w * x) / sum(w)
}

mse_outliers <- function(mat, s_max = 3) {
  mat <- mat[rowSums(mat != 0) > ncol(mat) / 2 + 0.5 & apply(mat, 1, sd) > 0, , drop = FALSE]
  L   <- ncol(mat)
  pen <- function(n, s) sqrt(2) * s * log10(factorial(n)) / n
  out <- t(apply(mat, 1, function(x) {
    z <- as.numeric(scale(x)); ord <- order(z); zs <- z[ord]
    models <- list(integer(L)); U <- L * log10(sd(zs))            # model: no outliers
    for (s in 1:s_max) {
      lo <- integer(L); lo[ord[1:s]] <- -1L                         # s lowest are outliers
      hi <- integer(L); hi[ord[(L - s + 1):L]] <- 1L                # s highest are outliers
      models <- c(models, list(lo, hi))
      U <- c(U, (L - s) * log10(sd(zs[(s + 1):L])) + pen(L - s, s),
                (L - s) * log10(sd(zs[1:(L - s)])) + pen(L - s, s))
    }
    both <- integer(L); both[ord[1]] <- -1L; both[ord[L]] <- 1L      # one low + one high
    models <- c(models, list(both))
    U <- c(U, (L - 2) * log10(sd(zs[2:(L - 1)])) + pen(L - 2, 2))
    U[!is.finite(U)] <- Inf
    models[[which.min(U)]]
  }))
  dimnames(out) <- dimnames(mat)
  out
}


# Control: genes that change over time without heat are excluded

drift <- sapply(genotypes, function(g) {
  o <- mse_outliers(ctrl[[g]][universe, ], s_max)
  x <- log2(ctrl[[g]][rownames(o), , drop = FALSE] + 1)
  d <- setNames(numeric(length(universe)), universe)
  for (i in which(rowSums(o != 0) > 0)) {
    shift <- abs(x[i, ] - median(x[i, o[i, ] == 0]))[o[i, ] != 0]
    if (any(shift >= ctrl_log2fc)) d[rownames(o)[i]] <- max(shift)
  }
  d
})                                   # genes x genotypes: largest control shift (0 = stable)
drift <- apply(drift, 1, max)        # a gene varying in either genotype counts for both


# Heat: MSE outliers that also differ from control (fold change + q)

results <- list()
for (g in genotypes) {
  sig <- abs(lfc[[g]]) >= min_log2fc & qval[[g]] < q_cutoff
  sig[, 1] <- FALSE                                          # baseline is never a call
  heat_effect <- apply(ifelse(sig, abs(lfc[[g]]), 0), 1, max)
  excluded <- universe[drift[universe] > 0 & heat_effect[universe] < rescue_ratio * drift[universe]]
  genes <- setdiff(universe, excluded)

  o  <- mse_outliers(heat[[g]][genes, ], s_max)
  fc <- lfc[[g]][rownames(o), , drop = FALSE]
  calls <- o * ((o == 1 & fc >= min_log2fc) | (o == -1 & fc <= -min_log2fc))
  # sustained responses: MSE flags the few unstressed timepoints as the outliers;
  # then the response is the remaining timepoints
  flip <- rowSums(calls != 0) == 0 & rowSums(o != 0) > 0
  calls[flip, ] <- (o[flip, ] == 0) * ((fc[flip, ] >= min_log2fc) - (fc[flip, ] <= -min_log2fc))
  calls <- calls * sig[rownames(o), ]
  calls <- calls[rowSums(calls != 0) > 0, , drop = FALSE]

  res <- data.frame(gene = rownames(calls), genotype = g,
                    n_sig_timepoints = rowSums(calls != 0),
                    main_set = rowSums(calls != 0) >= min_sig_tp,
                    direction = ifelse(rowSums(calls == 1) > 0 & rowSums(calls == -1) > 0, "mixed",
                                       ifelse(rowSums(calls == 1) > 0, "up", "down")),
                    setNames(as.data.frame(calls), paste0("call_", timepoints)),
                    setNames(as.data.frame(round(lfc[[g]][rownames(calls), ], 3)), paste0("log2FC_", timepoints)),
                    setNames(as.data.frame(signif(qval[[g]][rownames(calls), ], 3)), paste0("q_", timepoints)))
  res <- res[order(!res$main_set, -res$n_sig_timepoints), ]
  write.csv(res, file.path(out_dir, sprintf("MSE_%s.csv", g)), row.names = FALSE)
  results[[g]] <- res
  cat(sprintf("%s: %d control-variable genes excluded; %d genes called; main set (>= %d timepoints): %d\n",
              g, length(excluded), nrow(res), min_sig_tp, sum(res$main_set)))
}
main <- do.call(rbind, lapply(results, function(r) r[r$main_set, ]))
write.csv(main, file.path(out_dir, "MSE_main_genes.csv"), row.names = FALSE)
cat("Main set, either genotype:", length(unique(main$gene)), "genes ->", out_dir, "\n")
