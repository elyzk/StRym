# The generalized analysis steps (PCA, Shannon entropy, MSE, trajectories) on Tux output and on a counts table.
for (f in c("pca", "shannon_entropy", "modified_shannon_entropy", "trajectory")) source(file.path(here, "analysis", paste0(f, ".R")))
logs <- new.env(); logs$x <- character()
step_ctx <- function(d, out, prereq = list()) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE); d <- strym_normalize_design(d); logs$x <- character()
  list(design = d, data_path = d$data_path, sample_sheet = d$sample_sheet, out_dir = out, project_dir = d$project_dir, app_dir = here,
       prereq = prereq, log = function(...) logs$x <- c(logs$x, paste0(...)))
}

test_that("PCA reproduces the calculation of the original pca_plot_script.R on the same Cuffdiff files", {
  root <- tempfile("st"); dir.create(root); syn <- synth_tux(root, orig = FALSE)
  # the original script's data preparation, line by line (data.table + prcomp)
  f <- file.path(syn$project, "01_quantification", "tux", "diff_out", "genes.read_group_tracking")
  d <- data.table::fread(f, select = c("tracking_id", "condition", "replicate", "FPKM")); d[, sample := paste0(condition, "_R", replicate + 1)]
  m <- as.matrix(data.table::dcast(d, tracking_id ~ sample, value.var = "FPKM", fun.aggregate = sum, fill = 0), rownames = "tracking_id")
  m <- log2(m[rowMeans(m) >= 1, ] + 1); m <- m[order(apply(m, 1, var), decreasing = TRUE)[seq_len(min(500, nrow(m)))], ]   # (the original's [1:500] fails with fewer genes)
  ref <- prcomp(t(m)); pve <- 100 * ref$sdev^2 / sum(ref$sdev^2)
  out <- file.path(root, "pca"); res <- strym_step_pca(step_ctx(synth_design(syn), out))
  expect_equal(res$status, "ok"); expect_true(all(file.exists(file.path(out, res$files))))
  v <- read.csv(file.path(out, "pca_variance.csv")); expect_equal(v$percent_variance[1:5], round(pve[1:5], 2)); expect_equal(sum(v$percent_variance), 100, tolerance = 0.1)
  sc <- read.csv(file.path(out, "pca_scores.csv")); sc$key <- sub("_r([0-9]+)$", "_R\\1", sc$sample)
  expect_equal(nrow(sc), 60); expect_setequal(sc$key, rownames(ref$x))
  for (pc in c("PC1", "PC2", "PC3")) expect_equal(abs(sc[[pc]][match(rownames(ref$x), sc$key)]), unname(abs(ref$x[, pc])), tolerance = 1e-6)   # signs of a PC are arbitrary
  expect_setequal(names(sc)[1:4], c("sample", "condition", "time_point", "replicate"))
})

test_that("PCA can leave samples out, and needs at least 3", {
  root <- tempfile("st"); dir.create(root); syn <- synth_tux(root, orig = FALSE)
  drop <- "A_Heated_T2_r1"
  d <- synth_design(syn, params = list(drop_reps = drop))
  out <- file.path(root, "pca"); strym_step_pca(step_ctx(d, out))
  sc <- read.csv(file.path(out, "pca_scores.csv")); expect_equal(nrow(sc), 59); expect_false(drop %in% sc$sample)
  keep <- unique(read.csv(file.path(syn$project, "01_quantification", "sample_sheet.csv"))$sample)
  d2 <- synth_design(syn, params = list(drop_reps = paste(keep[-(1:2)], collapse = ", ")))
  expect_error(strym_step_pca(step_ctx(d2, file.path(root, "pca2"))), "at least 3 samples")
})

test_that("the entropy steps write their tables on StRym's Tux layout, with gene names and the comparisons used", {
  root <- tempfile("st"); dir.create(root); syn <- synth_tux(root, orig = FALSE); d <- synth_design(syn)
  out <- file.path(root, "mse"); r <- strym_step_modified_shannon_entropy(step_ctx(d, out))
  expect_equal(r$status, "ok"); expect_match(r$message, "2 treated-vs-control comparisons")
  expect_true(all(file.exists(file.path(out, c("MSE_A_Heated.csv", "MSE_B_Heated.csv", "MSE_main_genes.csv", "MSE_comparisons_used.csv", "MSE_calls_per_timepoint.png")))))
  a <- read.csv(file.path(out, "MSE_A_Heated.csv")); expect_true(all(c("gene", "gene_name", "comparison", "n_sig_timepoints", "main_set", "direction", "call_T0", "log2FC_T4", "q_T4") %in% names(a)))
  expect_true(all(grepl("^GENE[0-9]{4}$", a$gene_name))); expect_true(all(a$comparison == "A_Heated_vs_A_Control")); expect_true(all(a$call_T0 == 0))   # baseline
  cu <- read.csv(file.path(out, "MSE_comparisons_used.csv")); expect_equal(cu$treated, c("A_Heated", "B_Heated")); expect_equal(cu$control, c("A_Control", "B_Control"))
  # the planted truth: genes that go up after the baseline in both genotypes are found, drifting genes are not
  main <- read.csv(file.path(out, "MSE_main_genes.csv")); up <- names(syn$kind)[syn$kind == "up"]; dr <- names(syn$kind)[syn$kind == "drift"]
  expect_gt(mean(up %in% main$gene), 0.6); expect_lt(mean(dr %in% main$gene), 0.1)
  s <- strym_step_shannon_entropy(step_ctx(d, file.path(root, "sh")))
  expect_equal(s$status, "ok"); expect_true(file.exists(file.path(root, "sh", "shannon_main_genes.csv")))
  expect_gt(mean(up %in% read.csv(file.path(root, "sh", "shannon_main_genes.csv"))$gene), 0.6)
})

test_that("the entropy steps explain themselves when the data cannot support them", {
  root <- tempfile("st"); dir.create(root); demo <- strym_make_demo_data(root)
  cd <- strym_default_design(); cd$project_dir <- file.path(root, "p"); cd$data_path <- demo$counts; cd$sample_sheet <- demo$sample_sheet; cd$n_conditions <- 2L; cd$n_timepoints <- 3L
  expect_equal(strym_step_shannon_entropy(step_ctx(cd, file.path(root, "o1")))$status, "ok")               # a counts table with replicates: q-values are calculated with DESeq2
  expect_equal(strym_step_modified_shannon_entropy(step_ctx(cd, file.path(root, "o2")))$status, "ok")
  cd$expr_unit <- "normalized"                                                                              # normalized values: no q-values can be made
  expect_error(strym_step_shannon_entropy(step_ctx(cd, file.path(root, "o1b"))), "q-values")
  expect_error(strym_step_modified_shannon_entropy(step_ctx(cd, file.path(root, "o2b"))), "q-values")
  syn <- synth_tux(file.path(root, "s"), n_tp = 2, orig = FALSE, n_genes = 150)
  expect_error(strym_step_modified_shannon_entropy(step_ctx(synth_design(syn), file.path(root, "o3"))), "at least 3 time points")
  r <- strym_step_shannon_entropy(step_ctx(synth_design(syn), file.path(root, "o4"))); expect_equal(r$status, "ok")     # Shannon works with 2
  d1 <- synth_design(syn, params = list(control = "Nope")); expect_error(strym_step_shannon_entropy(step_ctx(d1, file.path(root, "o5"))), "not one of the conditions")
  d2 <- synth_design(syn, params = list(pairs = "A_Heated:A_Control")); r2 <- strym_step_shannon_entropy(step_ctx(d2, file.path(root, "o6")))   # one pair only
  expect_match(r2$message, "1 treated-vs-control comparison;")
})

test_that("trajectories plot the entropy main set, or your own gene list (ids and names, with a note about the missing)", {
  root <- tempfile("st"); dir.create(root); syn <- synth_tux(root, orig = FALSE); d <- synth_design(syn)
  mse <- file.path(root, "02", "mse"); strym_step_modified_shannon_entropy(step_ctx(d, mse))
  out <- file.path(root, "traj"); r <- strym_step_trajectory(step_ctx(d, out, prereq = list(modified_shannon_entropy = mse)))
  expect_equal(r$status, "ok"); expect_true(all(file.exists(file.path(out, c("trajectories.png", "trajectories.pdf", "trajectory_values.csv")))))
  tv <- read.csv(file.path(out, "trajectory_values.csv")); main <- read.csv(file.path(mse, "MSE_main_genes.csv"))
  expect_true(all(tv$gene %in% main$gene)); expect_lte(length(unique(tv$gene)), 200)
  expect_setequal(unique(tv$panel), c("A_Heated vs A_Control", "B_Heated vs B_Control")); expect_setequal(unique(tv$role), c("Control", "A_Heated", "B_Heated"))
  expect_setequal(unique(tv$time_point), syn$tps)
  # own list: an id, a name (any case), a gene that is not there, a header line
  gl <- file.path(root, "genes.csv"); writeLines(c("gene", syn$genes[1], tolower(syn$gene_names[2]), "NOT_A_GENE", syn$genes[1]), gl)
  out2 <- file.path(root, "traj2"); r2 <- strym_step_trajectory(step_ctx(synth_design(syn, params = list(trajectory_genes = gl)), out2))
  tv2 <- read.csv(file.path(out2, "trajectory_values.csv")); expect_setequal(unique(tv2$gene), syn$genes[1:2])
  expect_match(r2$message, "Plotted 2 of 3 genes from your gene list"); expect_true(any(grepl("1 of 3 genes not found.*NOT_A_GENE", logs$x)))
  expect_true(all(grepl("^GENE[0-9]{4} \\(XLOC_[0-9]{6}\\)$", unique(tv2$label))))                              # readable names with the id
  # nothing to plot
  expect_error(strym_step_trajectory(step_ctx(d, file.path(root, "t3"))), "No genes to plot")
  writeLines("zzz", gl); expect_error(strym_step_trajectory(step_ctx(synth_design(syn, params = list(trajectory_genes = gl)), file.path(root, "t4"))), "None of the 1 genes")
  expect_error(strym_step_trajectory(step_ctx(synth_design(syn, params = list(trajectory_genes = file.path(root, "nope.txt"))), file.path(root, "t5"))), "not found")
})

test_that("PCA and trajectories also work from a counts table (counts per million)", {
  root <- tempfile("st"); dir.create(root); demo <- strym_make_demo_data(root)
  d <- strym_default_design(); d$project_dir <- file.path(root, "p"); d$data_path <- demo$counts; d$sample_sheet <- demo$sample_sheet
  d$n_conditions <- 2L; d$condition_names <- c("control", "heat"); d$n_timepoints <- 3L; d$timepoint_names <- c("0h", "1h", "4h")
  data <- strym_load_counts(demo$counts, demo$sample_sheet)
  expect_equal(data$unit, "CPM"); expect_equal(nrow(data$samples), 12); expect_equal(sort(unique(data$samples$condition)), c("control", "heat"))
  cpm <- strym_wide(data); expect_equal(unname(colSums(cpm)), rep(1e6, 12), tolerance = 1e-9)                    # library-size normalized
  out <- file.path(root, "pca"); r <- strym_step_pca(step_ctx(d, out)); expect_equal(r$status, "ok"); expect_equal(nrow(read.csv(file.path(out, "pca_scores.csv"))), 12)
  expect_match(r$message, "CPM")
  gl <- file.path(root, "g.txt"); writeLines(c("GENE001", "GENE002"), gl); d$params <- list(trajectory_genes = gl)
  out2 <- file.path(root, "tr"); r2 <- strym_step_trajectory(step_ctx(d, out2)); expect_equal(r2$status, "ok")
  tv <- read.csv(file.path(out2, "trajectory_values.csv")); expect_setequal(unique(tv$time_point), c("0h", "1h", "4h")); expect_setequal(unique(tv$role), c("Control", "heat"))
  # a counts table without a sample sheet: every sample is its own condition, PCA still works
  d$sample_sheet <- ""; r3 <- strym_step_pca(step_ctx(d, file.path(root, "pca3"))); expect_equal(r3$status, "ok")
})
