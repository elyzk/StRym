# Single cell / single nucleus: the Seurat step (generalized from the baseline notebooks), checked against synthetic data with known truth.
skip_if_not_installed("Seurat"); skip_if_not_installed("harmony")
source(file.path(here, "analysis", "seurat_clustering.R"))
sc_ctx <- function(path, out, ..., sample_sheet = "") {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  d <- strym_default_design(); d$data_type <- "single_cell"; d$data_kind <- "10x"; d$data_path <- path; d$project_dir <- dirname(out)
  a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]
  logs <- new.env(); logs$x <- character()
  list(design = strym_normalize_design(d), data_path = path, sample_sheet = sample_sheet, out_dir = out, project_dir = dirname(out), app_dir = here, prereq = list(),
       log = function(...) logs$x <- c(logs$x, paste0(...)), logs = logs)
}

test_that("10x folders are found (single matrix, or one per sample) and conditions come from the sample names or the sheet", {
  root <- tempfile("sc"); dir.create(root); syn <- synth_10x(root, cells = 20, n_genes = 100)
  f <- .find_10x(root); expect_equal(sort(f$sample), sort(syn$samples)); expect_equal(basename(f$dir), rep("raw_matrix", 4))     # "raw_matrix" is not a sample name
  one <- .find_10x(file.path(root, "Control-1", "raw_matrix")); expect_equal(one$sample, "Control-1")
  nothing <- file.path(root, "nothing"); dir.create(nothing); expect_equal(nrow(.find_10x(nothing)), 0)
  expect_equal(.sc_conditions(c("Control-1", "FCP1.2", "WT_3", "Solo"), ""), c("Control", "FCP1", "WT", "Solo"))
  ss <- file.path(root, "ss.csv"); write.csv(data.frame(sample = c("control-1", "Treated-1"), condition = c("ctrl", "heat")), ss, row.names = FALSE)
  expect_equal(.sc_conditions(c("Control-1", "Control-2", "Treated-1"), ss), c("ctrl", "Control", "heat"))            # the sheet wins where it knows the sample
})

test_that("the elbow rule of the notebooks picks a sensible number of PCs", {
  expect_equal(.elbow_pcs(c(10, 9, 8, 1, 0.9, 0.8, 0.7, 0.6)), 4)                          # the rule: first PC past 90% cumulative variance (4th) vs last drop > 0.1% (8th) -> the smaller
  expect_gte(.elbow_pcs(rep(1, 30)), 2); expect_lte(.elbow_pcs(rep(1, 30)), 30)
})

test_that("Seurat step: clusters match the planted cell types, markers and condition effects are found (CCA integration, UMAP + t-SNE)", {
  root <- tempfile("sc"); dir.create(root); syn <- synth_10x(root)                                        # one sample gz, the next plain: both read
  out <- file.path(root, "out", "seurat"); ctx <- sc_ctx(syn$root, out, plots = c("umap", "tsne"), n_conditions = 2L)
  r <- strym_step_seurat_clustering(ctx)
  expect_equal(r$status, "ok"); expect_true(all(file.exists(file.path(out, r$files))))
  expect_true(all(c("umap.png", "tsne.png", "umap_by_condition.png", "qc_violin.png", "elbow.png", "cells_per_sample.csv", "cluster_markers.csv",
                    "condition_DEGs_by_cluster.csv", "cluster_composition.csv", "seurat_object.rds") %in% r$files))
  expect_match(r$message, "4 samples \\(2 conditions\\)"); expect_match(r$message, "UMAP \\+ TSNE")
  qc <- read.csv(file.path(out, "cells_per_sample.csv")); expect_equal(qc$sample, syn$samples); expect_equal(qc$condition, c("Control", "Control", "Treated", "Treated"))
  expect_true(all(qc$cells_kept > 150)); expect_true(any(grepl("Integrated 4 samples with CCA", ctx$logs$x)))
  obj <- readRDS(file.path(out, "seurat_object.rds")); truth <- sub("-.*", "", colnames(obj)); truth <- sub("^[^_]*_", "", truth)
  type <- sub("^.*_(ct[0-9])-.*$", "\\1", colnames(obj)); expect_true(all(grepl("^ct[1-3]$", type)))
  tab <- table(cluster = obj$clusters, type = type)
  purity <- sum(apply(tab, 1, max)) / sum(tab); expect_gt(purity, 0.95)                                    # clusters are (nearly) pure cell types
  expect_gte(length(unique(obj$clusters)), 3); expect_lte(length(unique(obj$clusters)), 8)
  expect_true(all(c("umap", "tsne") %in% names(obj@reductions)))
  # markers: the top gene of each cluster comes from that cell type's marker block
  mk <- read.csv(file.path(out, "cluster_markers.csv")); expect_equal(names(mk)[1:2], c("gene", "cluster"))
  dom <- apply(tab, 1, function(x) colnames(tab)[which.max(x)])
  for (cl in rownames(tab)) {
    top <- head(mk[mk$cluster == cl, "gene"], 5); k <- as.integer(sub("ct", "", dom[[cl]])); idx <- as.integer(sub("Gene", "", top))
    expect_true(mean(idx >= (k - 1) * 30 + 1 & idx <= k * 30) >= 0.6, info = paste("cluster", cl))
  }
  # condition effect: Gene700-720 higher in Treated, inside the cell-type-2 cluster(s)
  dg <- read.csv(file.path(out, "condition_DEGs_by_cluster.csv")); expect_true(all(dg$comparison == "Treated_vs_Control"))
  expect_true(all(dg$p_val_adj < 0.05))                                                                     # only significant genes are kept
  hit <- dg[dg$gene %in% paste0("Gene", 700:720), ]; expect_gt(length(unique(hit$gene)), 10); expect_true(all(hit$avg_log2FC > 0))
  c2 <- rownames(tab)[dom == "ct2"]; expect_true(all(hit$cluster %in% c2))
  expect_gt(mean(dg$gene %in% paste0("Gene", 700:720)), 0.5)                                                # ... and they are the planted ones
  comp <- read.csv(file.path(out, "cluster_composition.csv")); expect_equal(sum(comp[, -1]), ncol(obj))
})

test_that("Seurat step: other integration choices, one sample, an .rds object, and clear errors", {
  root <- tempfile("sc"); dir.create(root); syn <- synth_10x(root, cells = 100)
  # Harmony, UMAP only, a typed number of PCs, a stricter UMI filter
  d1 <- file.path(root, "o1"); c1 <- sc_ctx(syn$root, d1, plots = "umap", params = list(sc_integration = "harmony", sc_npcs = 8, sc_min_counts = 400, sc_resolution = 0.5))
  r1 <- strym_step_seurat_clustering(c1); expect_equal(r1$status, "ok"); expect_false(file.exists(file.path(d1, "tsne.png"))); expect_true(file.exists(file.path(d1, "umap.png")))
  expect_true(any(grepl("Using 8 principal components \\(set in Analysis", c1$logs$x))); expect_true(any(grepl("Integrated 4 samples with HARMONY", c1$logs$x)))
  # no integration: samples only merged
  c2 <- sc_ctx(syn$root, file.path(root, "o2"), plots = "umap", params = list(sc_integration = "none")); r2 <- strym_step_seurat_clustering(c2); expect_equal(r2$status, "ok")
  expect_true(any(grepl("Integration switched off", c2$logs$x)))
  # one sample (a single matrix folder): no integration, no condition comparison
  c3 <- sc_ctx(file.path(root, "Control-1", "raw_matrix"), file.path(root, "o3")); r3 <- strym_step_seurat_clustering(c3); expect_equal(r3$status, "ok")
  expect_false(file.exists(file.path(root, "o3", "condition_DEGs_by_cluster.csv"))); expect_false(file.exists(file.path(root, "o3", "umap_by_condition.png")))
  # a Seurat object saved earlier (several samples in it) is processed too
  rds <- file.path(root, "obj.rds"); file.copy(file.path(root, "o2", "seurat_object.rds"), rds)
  c4 <- sc_ctx(rds, file.path(root, "o4"), plots = "umap"); r4 <- strym_step_seurat_clustering(c4); expect_equal(r4$status, "ok"); expect_match(r4$message, "4 samples")
  expect_true(any(grepl("Integrated 4 samples with CCA", c4$logs$x)))                                      # the object's samples were split for integration
  # clear errors
  expect_error(strym_step_seurat_clustering(sc_ctx(file.path(root, "nope"), file.path(root, "e1"))), "Choose your single cell")
  empty <- file.path(root, "empty"); dir.create(empty); expect_error(strym_step_seurat_clustering(sc_ctx(empty, file.path(root, "e2"))), "No 10x matrix")
  expect_error(strym_step_seurat_clustering(sc_ctx(syn$root, file.path(root, "e3"), params = list(sc_min_counts = 1e6))), "Only 0 cells pass the UMI filter")
  bad <- file.path(root, "bad.rds"); saveRDS(1:3, bad); expect_error(strym_step_seurat_clustering(sc_ctx(bad, file.path(root, "e4"))), "not a Seurat object")
})

test_that("Seurat step runs through strym_run for a single cell design and is blocked for whole organ data", {
  root <- tempfile("sc"); dir.create(root); syn <- synth_10x(root, samples = c("Control-1", "Treated-1"), cells = 100)
  d <- strym_default_design(); d$project_name <- "sc"; d$project_dir <- file.path(root, "proj"); d$data_type <- "single_cell"; d$data_kind <- "10x"; d$data_path <- syn$root
  d$n_conditions <- 2L; d$plots <- "umap"
  res <- quiet(strym_run(d, steps = "seurat_clustering", app_dir = here)); expect_equal(res$summary$status, "ok", info = res$summary$message)
  expect_true(file.exists(file.path(d$project_dir, "03_clustering", "seurat_umap_tsne", "umap.png")))
  expect_true(strym_plan(d)$applicable[strym_plan(d)$id == "seurat_clustering"])
  bulk <- d; bulk$data_type <- "bulk"; bulk$data_kind <- "counts"; expect_false(strym_plan(bulk)$applicable[strym_plan(bulk)$id == "seurat_clustering"])
})
