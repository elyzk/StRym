# DESeq2 step: explicit treated-vs-control contrasts per time point, checked against planted differences and against DESeq2 called directly.
skip_if_not_installed("DESeq2")
source(file.path(here, "analysis", "deseq2.R"))
ds_ctx <- function(path, sheet, out, ..., tps = c("0h", "1h", "4h")) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  d <- strym_default_design(); d$data_path <- path; d$sample_sheet <- sheet; d$project_dir <- dirname(out); d$timepoint_names <- tps; d$n_timepoints <- length(tps)
  a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]
  logs <- new.env(); logs$x <- character()
  list(design = strym_normalize_design(d), data_path = path, sample_sheet = sheet, out_dir = out, project_dir = dirname(out), app_dir = here, prereq = list(),
       log = function(...) logs$x <- c(logs$x, paste0(...)), logs = logs)
}
# Control / Heat x 3 time points x 3 replicates, negative-binomial counts; 40 genes up and 40 down in Heat from 1h on
synth_counts <- function(root, n_genes = 400, reps = 3, seed = 3, conds = c("Control", "Heat")) {
  set.seed(seed); tps <- c("0h", "1h", "4h"); base <- exp(rnorm(n_genes, log(80), 1)); genes <- sprintf("G%04d", seq_len(n_genes))
  up <- 1:40; down <- 41:80
  S <- expand.grid(replicate = seq_len(reps), time_point = tps, condition = conds, stringsAsFactors = FALSE)
  S$sample <- sprintf("%s_%s_r%d", S$condition, S$time_point, S$replicate)
  m <- sapply(seq_len(nrow(S)), function(j) {
    eff <- rep(0, n_genes); if (S$condition[j] == conds[2] && S$time_point[j] != "0h") { eff[up] <- 2; eff[down] <- -2 }
    stats::rnbinom(n_genes, mu = base * 2^eff, size = 20) })
  dimnames(m) <- list(genes, S$sample)
  cf <- file.path(root, "counts.csv"); utils::write.csv(data.frame(gene = genes, m, check.names = FALSE), cf, row.names = FALSE)
  sf <- file.path(root, "sheet.csv"); utils::write.csv(S[, c("sample", "condition", "time_point", "replicate")], sf, row.names = FALSE)
  list(counts = cf, sheet = sf, genes = genes, up = genes[up], down = genes[down], S = S, m = m)
}

test_that("DESeq2 finds the planted genes in the right direction at each time point, and matches DESeq2 run by hand", {
  root <- tempfile("ds"); dir.create(root); sy <- synth_counts(root)
  out <- file.path(root, "out", "deseq2"); ctx <- ds_ctx(sy$counts, sy$sheet, out)
  r <- strym_step_deseq2(ctx); expect_equal(r$status, "ok")
  expect_true(all(file.exists(file.path(out, r$files))))
  expect_true(all(c("deseq2_Heat_vs_Control_0h.csv", "deseq2_Heat_vs_Control_1h.csv", "deseq2_Heat_vs_Control_4h.csv", "deseq2_DEGs.csv", "normalized_counts.csv",
                    "size_factors.csv", "deseq2_comparisons_used.csv", "deseq2_volcano.png", "deseq2_DEG_heatmap.png") %in% r$files))
  r4 <- read.csv(file.path(out, "deseq2_Heat_vs_Control_4h.csv")); r0 <- read.csv(file.path(out, "deseq2_Heat_vs_Control_0h.csv"))
  expect_gt(mean(sy$up %in% r4$gene[r4$direction == "up"]), 0.9); expect_gt(mean(sy$down %in% r4$gene[r4$direction == "down"]), 0.9)
  expect_lt(sum(r4$direction != "none" & !(r4$gene %in% c(sy$up, sy$down))), 5)                               # almost no false positives
  expect_lt(sum(r0$direction != "none"), 5)                                                                   # nothing happens at 0h
  expect_true(all(r4$log2FoldChange[r4$gene %in% sy$up & r4$direction == "up"] >= 1))                         # positive = higher in the treated condition
  expect_gt(median(r4$log2FoldChange[r4$gene %in% sy$up]), 1.5); expect_lt(median(r4$log2FoldChange[r4$gene %in% sy$down]), -1.5)
  deg <- read.csv(file.path(out, "deseq2_DEGs.csv")); expect_equal(names(deg)[1:4], c("gene", "comparison", "time_point", "direction")); expect_setequal(unique(deg$time_point), c("1h", "4h"))
  # DESeq2 called directly with the same filter, design and contrast
  suppressPackageStartupMessages(library(DESeq2))
  S <- sy$S; m <- sy$m; keep <- rowSums(m) >= 10; S$group <- paste(S$condition, S$time_point, sep = "_")
  dds <- DESeq(DESeqDataSetFromMatrix(m[keep, ], data.frame(group = factor(S$group), row.names = S$sample), ~ group), quiet = TRUE)
  ref <- as.data.frame(results(dds, contrast = c("group", "Heat_4h", "Control_4h"), alpha = 0.05)); ref$padj[is.na(ref$padj)] <- 1
  mine <- r4[match(rownames(ref), r4$gene), ]
  expect_equal(mine$log2FoldChange, ref$log2FoldChange); expect_equal(mine$padj, ref$padj)
  nc <- read.csv(file.path(out, "normalized_counts.csv"), check.names = FALSE); expect_equal(dim(nc), c(sum(keep), 19))
  expect_equal(read.csv(file.path(out, "size_factors.csv"))$size_factor, round(unname(sizeFactors(dds)), 4))
  cu <- read.csv(file.path(out, "deseq2_comparisons_used.csv")); expect_equal(cu$g_treated, paste0("Heat_", c("0h", "1h", "4h"))); expect_equal(cu$g_control, paste0("Control_", c("0h", "1h", "4h")))
  expect_match(r$message, "3 comparisons"); expect_true(any(grepl("Heat_vs_Control_4h: [0-9]+ up, [0-9]+ down", ctx$logs$x)))
})

test_that("DESeq2 with one time point, cutoffs from the settings, and samples left out", {
  root <- tempfile("ds"); dir.create(root); sy <- synth_counts(root)
  one <- sy$S[sy$S$time_point == "4h", ]; cf <- file.path(root, "one.csv"); sf <- file.path(root, "one_sheet.csv")
  write.csv(data.frame(gene = sy$genes, sy$m[, one$sample], check.names = FALSE), cf, row.names = FALSE); write.csv(one[, c("sample", "condition", "time_point", "replicate")], sf, row.names = FALSE)
  out <- file.path(root, "o1", "deseq2"); r <- strym_step_deseq2(ds_ctx(cf, sf, out, tps = "4h"))
  expect_true(file.exists(file.path(out, "deseq2_Heat_vs_Control.csv"))); expect_match(r$message, "1 comparison;")        # no time point in the name with one time point
  strict <- strym_step_deseq2(ds_ctx(cf, sf, file.path(root, "o2", "deseq2"), tps = "4h", params = list(min_log2fc = 3, q_cutoff = 0.001)))
  expect_lt(read.csv(file.path(root, "o2", "deseq2", "deseq2_DEGs.csv"))$gene %>% length, nrow(read.csv(file.path(out, "deseq2_DEGs.csv"))) + 1)
  expect_equal(nrow(read.csv(file.path(root, "o2", "deseq2", "deseq2_DEGs.csv"))), 0)                                  # nothing reaches a 8-fold change here (planted: 4-fold)
  dr <- strym_step_deseq2(ds_ctx(cf, sf, file.path(root, "o3", "deseq2"), tps = "4h", params = list(drop_reps = "Heat_4h_r1")))
  expect_equal(ncol(read.csv(file.path(root, "o3", "deseq2", "normalized_counts.csv"), check.names = FALSE)), 1 + 5)
})

test_that("DESeq2 without a control condition compares the 2nd condition with the 1st and says so", {
  root <- tempfile("ds"); dir.create(root); sy <- synth_counts(root, conds = c("B", "D"))
  out <- file.path(root, "out", "deseq2"); ctx <- ds_ctx(sy$counts, sy$sheet, out)
  r <- strym_step_deseq2(ctx); expect_equal(r$status, "ok")
  expect_true(any(grepl("No control condition found: comparing D with B", ctx$logs$x))); expect_true(file.exists(file.path(out, "deseq2_D_vs_B_4h.csv")))
  r4 <- read.csv(file.path(out, "deseq2_D_vs_B_4h.csv")); expect_gt(mean(sy$up %in% r4$gene[r4$direction == "up"]), 0.9)
})

test_that("DESeq2 explains what it cannot do", {
  root <- tempfile("ds"); dir.create(root); sy <- synth_counts(root)
  frac <- file.path(root, "frac.csv"); write.csv(data.frame(gene = sy$genes, sy$m / 3, check.names = FALSE), frac, row.names = FALSE)
  expect_error(strym_step_deseq2(ds_ctx(frac, sy$sheet, file.path(root, "e1"))), "RAW counts.*normalized")
  neg <- sy$m; neg[1, 1] <- -1; ng <- file.path(root, "neg.csv"); write.csv(data.frame(gene = sy$genes, neg, check.names = FALSE), ng, row.names = FALSE)
  expect_error(strym_step_deseq2(ds_ctx(ng, sy$sheet, file.path(root, "e2"))), "negative")
  s1 <- sy$S[!duplicated(paste(sy$S$condition, sy$S$time_point)), ]; c1 <- file.path(root, "c1.csv"); h1 <- file.path(root, "h1.csv")
  write.csv(data.frame(gene = sy$genes, sy$m[, s1$sample], check.names = FALSE), c1, row.names = FALSE); write.csv(s1[, c("sample", "condition", "time_point", "replicate")], h1, row.names = FALSE)
  expect_error(strym_step_deseq2(ds_ctx(c1, h1, file.path(root, "e3"))), "needs replicates")
  only <- sy$S[sy$S$condition == "Heat", ]; c2 <- file.path(root, "c2.csv"); h2 <- file.path(root, "h2.csv")
  write.csv(data.frame(gene = sy$genes, sy$m[, only$sample], check.names = FALSE), c2, row.names = FALSE); write.csv(only[, c("sample", "condition", "time_point", "replicate")], h2, row.names = FALSE)
  expect_error(strym_step_deseq2(ds_ctx(c2, h2, file.path(root, "e4"))), "only one")
})

test_that("DESeq2 with two treated conditions compares each with the control", {
  root <- tempfile("ds"); dir.create(root); sy <- synth_counts(root, conds = c("Control", "Heat", "Cold"))     # only Heat has planted differences
  out <- file.path(root, "out", "deseq2"); r <- strym_step_deseq2(ds_ctx(sy$counts, sy$sheet, out)); expect_equal(r$status, "ok")
  cu <- read.csv(file.path(out, "deseq2_comparisons_used.csv")); expect_equal(unique(cu$control), "Control"); expect_setequal(unique(cu$treated), c("Heat", "Cold")); expect_equal(nrow(cu), 6)
  expect_gt(cu$up[cu$treated == "Heat" & cu$time_point == "4h"], 25); expect_lt(cu$up[cu$treated == "Cold" & cu$time_point == "4h"] + cu$down[cu$treated == "Cold" & cu$time_point == "4h"], 5)
})
