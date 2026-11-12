# Single cell: ONE 10x folder that holds several samples (told apart by the label in front of each barcode), damaged counts, empty samples,
# a low-memory matrix reader, and the plan / draft sample sheet for such data.
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
sheet <- function(root, labels, cond, name = "ss.csv") { f <- file.path(root, name); write.csv(data.frame(sample = labels, condition = cond), f, row.names = FALSE); f }

# Four synthetic samples (cell types and a condition effect planted, see helper-10x.R) put into ONE matrix with barcodes like Ctrl1-ACGT...(16 letters).
# damage = TRUE puts a missing, an infinite, a negative, an absurd (1e308) and a fractional value into the file; `empty` adds a sample with no counts at all.
combined_10x <- function(root, labels = c("Ctrl1", "Ctrl2", "Trt1", "Trt2"), cells = 100, damage = FALSE, empty = NULL, name = "combined") {
  syn <- synth_10x(file.path(root, "parts"), cells = cells); mats <- lapply(syn$dirs, .read_10x); set.seed(5)
  code <- function(n) vapply(seq_len(n), function(i) paste(sample(c("A", "C", "G", "T"), 16, TRUE), collapse = ""), "")
  cols <- lapply(seq_along(mats), function(i) { m <- mats[[i]]; colnames(m) <- paste0(labels[i], "-", code(ncol(m))); m })
  if (!is.null(empty)) { z <- Matrix::Matrix(0, nrow(mats[[1]]), 30, sparse = TRUE); dimnames(z) <- list(rownames(mats[[1]]), paste0(empty, "-", code(30))); cols[[length(cols) + 1]] <- z }
  m <- do.call(cbind, cols); d <- file.path(root, name); dir.create(d)
  tmp <- tempfile(fileext = ".mtx"); Matrix::writeMM(m, tmp); L <- readLines(tmp)
  if (damage) for (k in 1:5) L[2 + k] <- sub("[^ ]+$", c("NaN", "inf", "-4", "1e308", "2.5")[k], L[2 + k])
  con <- gzfile(file.path(d, "matrix.mtx.gz"), "wt"); writeLines(L, con); close(con)
  con <- gzfile(file.path(d, "features.tsv.gz"), "wt"); writeLines(sprintf("%s\t%s", rownames(m), rownames(m)), con); close(con)    # two columns only, as in many archives
  con <- gzfile(file.path(d, "barcodes.tsv.gz"), "wt"); writeLines(colnames(m), con); close(con)
  list(dir = d, labels = c(labels, empty), n = ncol(m), root = root, syn = syn)
}

test_that("sample labels are read from barcodes: prefixes, aggr library numbers, and no label at all", {
  expect_equal(strym_barcode_labels(c("SAMN15315800-AAACACGACGCTCTTC", "Ctrl_AAACCTGAGAAACCAT", "AAACCTGAGAAACCAT-2", "AAACCTGAGAAACCAT-1", "AAACCTGAGAAACCAT", "ct1-AAAC0001", NA)),
               c("SAMN15315800", "Ctrl", "2", "1", NA, NA, NA))
  bc <- function(l, n) paste0(l, "-", rep("AAACCTGAGAAACCAT", n))
  expect_equal(strym_split_labels(c(bc("A", 20), bc("B", 20))), c("A", "B"))                                      # two samples in the order of first appearance
  expect_equal(strym_split_labels(bc("A", 40)), character()); expect_equal(strym_split_labels(c(bc("A", 20), "AAACCTGAGAAACCAT")), character())
  expect_equal(strym_split_labels(sprintf("S%d-AAACCTGAGAAACCAT", 1:40)), character())                           # a "label" per cell is not a sample name
  expect_equal(strym_split_labels(paste0("AAACCTGAGAAACCAT-", rep(1:2, each = 20))), c("1", "2"))                 # cellranger aggr
  expect_equal(strym_split_labels(paste0("AAACCTGAGAAACCAT-", rep(1, 20))), character())                         # an ordinary run: one library
})

test_that("impossible counts are zeroed and counted; real counts are left alone", {
  m <- Matrix::sparseMatrix(i = 1:6, j = c(1, 1, 2, 2, 3, 3), x = c(NA, Inf, -3, 1e308, 2.5, 4), dims = c(6, 3), dimnames = list(letters[1:6], LETTERS[1:3]))
  cl <- .clean_counts(m); expect_equal(cl$bad, 4); expect_equal(cl$fractional, 1); expect_equal(cl$entries, 6)
  expect_equal(unname(Matrix::colSums(cl$m)), c(0, 0, 6.5)); expect_equal(length(cl$m@x), 2)                      # the zeros are not stored
  ok <- .clean_counts(Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(1, 2, 3))); expect_equal(c(ok$bad, ok$fractional), c(0, 0))
})

test_that(".read_mtx reads like Matrix::readMM (plain, gz, unsorted), drops entries outside the matrix, and falls back for other layouts", {
  set.seed(2); m <- Matrix::rsparsematrix(40, 25, 0.2, rand.x = function(n) rpois(n, 3) + 1); d <- tempfile(); dir.create(d)
  f <- file.path(d, "m.mtx"); Matrix::writeMM(m, f); L <- readLines(f); hdr <- L[1:2]; body <- L[-(1:2)]; ref <- as.matrix(m)
  expect_equal(as.matrix(.read_mtx(f)), ref)
  con <- gzfile(file.path(d, "m.mtx.gz"), "wt"); writeLines(L, con); close(con); expect_equal(as.matrix(.read_mtx(file.path(d, "m.mtx.gz"))), ref)
  writeLines(c(hdr, sample(body)), file.path(d, "u.mtx")); expect_equal(as.matrix(.read_mtx(file.path(d, "u.mtx"))), ref)
  writeLines(c(hdr, body, "999 1 5", "1 999 5"), file.path(d, "o.mtx")); expect_equal(as.matrix(.read_mtx(file.path(d, "o.mtx"))), ref)
  pm <- m != 0; Matrix::writeMM(pm, file.path(d, "p.mtx")); expect_equal(dim(.read_mtx(file.path(d, "p.mtx"))), c(40, 25))      # no values: the standard reader
  writeLines(c(hdr, gsub("e", "E", body)), file.path(d, "w.mtx")); expect_equal(dim(.read_mtx(file.path(d, "w.mtx"))), c(40, 25))
})

test_that("one matrix becomes its samples: split by the sheet, labels left out noted, empty samples skipped with the reason", {
  root <- tempfile("mc"); dir.create(root); cb <- combined_10x(root, empty = "Empty1"); m <- .read_10x(cb$dir); expect_equal(ncol(m), cb$n)
  p <- strym_normalize_params(list()); logs <- character(); lg <- function(...) logs <<- c(logs, paste0(...))
  pr <- .sc_prepare(m, "combined", sheet(root, c("ctrl1", "TRT1", "Empty1"), c("a", "b", "c")), p, lg)                     # labels are matched without regard to case
  expect_equal(names(pr$mats), c("Ctrl1", "Trt1")); expect_equal(pr$read, c(Ctrl1 = 100, Trt1 = 100, Empty1 = 30)[c("Ctrl1", "Trt1", "Empty1")])
  expect_equal(pr$skipped$sample, "Empty1"); expect_match(pr$skipped$reason, "none of its 30 cells has more than 500 UMIs \\(the most any cell has is 0\\)")
  expect_true(any(grepl("Left out \\(not in the sample sheet\\): Ctrl2, Trt2", logs))); expect_true(any(grepl("combined holds 5 samples; 3 of them are analyzed", logs)))
  expect_true(any(grepl("Sample Empty1 skipped: none of its 30 cells", logs)))
  expect_false(any(grepl("^(Ctrl1|Trt1)", colnames(pr$mats$Ctrl1)))); expect_true(all(nchar(colnames(pr$mats$Ctrl1)) == 16))   # the label is added back when merging
  expect_true(all(Matrix::colSums(pr$mats$Ctrl1) > 500))                                                            # cells below the UMI filter are gone before Seurat
  logs <- character(); pr0 <- .sc_prepare(m, "combined", "", p, lg)                                                  # no sheet: ONE sample, and the log says how to change that
  expect_equal(names(pr0$mats), "combined"); expect_true(any(grepl("5 sample labels \\(Ctrl1, Ctrl2, Trt1, Trt2, \\.\\.\\.\\).*there is no sample sheet.*ONE sample.*“sample” column", logs)))
  pr1 <- .sc_prepare(m, "combined", sheet(root, c("other"), "x", "o.csv"), p, lg); expect_equal(names(pr1$mats), "combined")   # a sheet that lists none of them
  expect_true(any(grepl("does not list any of them", logs)))
})

test_that("Seurat step: one matrix with four samples is split by the sample sheet, analyzed, and combined (CCA, UMAP)", {
  root <- tempfile("mc"); dir.create(root); cb <- combined_10x(root); ss <- sheet(root, cb$labels, c("ctrl", "ctrl", "heat", "heat"))
  out <- file.path(root, "out", "seurat"); ctx <- sc_ctx(cb$dir, out, plots = "umap", n_conditions = 2L, sample_sheet = ss)
  r <- strym_step_seurat_clustering(ctx); expect_equal(r$status, "ok"); expect_match(r$message, "4 samples \\(2 conditions\\)")
  qc <- read.csv(file.path(out, "cells_per_sample.csv")); expect_equal(qc$sample, cb$labels); expect_equal(qc$condition, c("ctrl", "ctrl", "heat", "heat")); expect_true(all(qc$cells_kept > 60))
  expect_true(any(grepl("combined holds 4 samples; all are analyzed separately", ctx$logs$x))); expect_true(any(grepl("Integrated 4 samples with CCA", ctx$logs$x)))
  obj <- readRDS(file.path(out, "seurat_object.rds")); expect_setequal(unique(obj$sample), cb$labels); expect_true(all(grepl("^(Ctrl1|Ctrl2|Trt1|Trt2)_", colnames(obj))))
  expect_true(file.exists(file.path(out, "condition_DEGs_by_cluster.csv")))
})

test_that("Seurat step: damaged counts and an empty sample no longer stop the analysis; without a sheet the cells stay one sample", {
  root <- tempfile("mc"); dir.create(root); cb <- combined_10x(root, damage = TRUE, empty = "Empty1")
  ss <- sheet(root, cb$labels, c("ctrl", "ctrl", "heat", "heat", "heat")); out <- file.path(root, "o1")
  ctx <- sc_ctx(cb$dir, out, plots = "umap", n_conditions = 2L, sample_sheet = ss, params = list(sc_integration = "none"))
  r <- strym_step_seurat_clustering(ctx); expect_equal(r$status, "ok"); expect_match(r$message, "4 samples")
  expect_true(any(grepl("entries that were missing, infinite, negative or above 10 million were set to 0", ctx$logs$x))); expect_true(any(grepl("not whole numbers", ctx$logs$x)))
  expect_true(any(grepl("Sample Empty1 skipped: none of its 30 cells has more than 500 UMIs", ctx$logs$x))); expect_true(any(grepl("1 sample skipped", ctx$logs$x)))
  qc <- read.csv(file.path(out, "cells_per_sample.csv")); expect_equal(qc$cells_kept[qc$sample == "Empty1"], 0); expect_equal(sum(qc$cells_kept > 0), 4)
  obj <- readRDS(file.path(out, "seurat_object.rds")); cnt <- SeuratObject::LayerData(obj, layer = "counts"); expect_true(all(is.finite(cnt@x))); expect_lt(max(cnt@x), 1e7); expect_gte(min(cnt@x), 0)
  # no sheet: the matrix is analyzed as one sample, with a log line that says how to split it
  c2 <- sc_ctx(cb$dir, file.path(root, "o2"), plots = "umap"); r2 <- strym_step_seurat_clustering(c2); expect_equal(r2$status, "ok"); expect_match(r2$message, "in 1 samples")
  expect_true(any(grepl("ONE sample", c2$logs$x))); expect_equal(nrow(read.csv(file.path(root, "o2", "cells_per_sample.csv"))), 1)
  # a sheet that lists only two of the labels: those two are analyzed, the others left out
  s3 <- sheet(root, c("Ctrl1", "Trt1"), c("ctrl", "heat"), "two.csv"); c3 <- sc_ctx(cb$dir, file.path(root, "o3"), plots = "umap", n_conditions = 2L, sample_sheet = s3, params = list(sc_integration = "none"))
  r3 <- strym_step_seurat_clustering(c3); expect_match(r3$message, "2 samples \\(2 conditions\\)"); expect_true(any(grepl("Left out \\(not in the sample sheet\\): Ctrl2, Trt2, Empty1", c3$logs$x)))
  # nothing passes the filter at all: the same plain error as before
  expect_error(strym_step_seurat_clustering(sc_ctx(cb$dir, file.path(root, "o4"), sample_sheet = ss, params = list(sc_min_counts = 1e6))), "Only 0 cells pass the UMI filter")
})

test_that("the plan and the draft sample sheet for single cell data", {
  root <- tempfile("mc"); dir.create(root); cb <- combined_10x(root)
  d <- strym_default_design(); d$data_type <- "single_cell"; d$data_kind <- "10x"; d$data_path <- cb$dir; d$n_conditions <- 2L
  a <- strym_auto_sample_sheet(d)                                                                                 # the samples inside the barcodes are listed for the researcher to check
  expect_equal(a$sheet$sample, cb$labels); expect_equal(names(a$sheet), c("sample", "condition", "time_point", "replicate")); expect_match(a$note, "4 samples were found inside the cell barcodes")
  d$condition_names <- c("ctrl", "heat"); a2 <- strym_auto_sample_sheet(d); expect_equal(nrow(a2$sheet), 4)
  syn <- cb$syn; d2 <- d; d2$data_path <- syn$root; d2$condition_names <- character()                              # one folder per sample: named after the folders
  a3 <- strym_auto_sample_sheet(d2); expect_equal(a3$sheet$sample, syn$samples); expect_equal(a3$sheet$condition, c("Control", "Control", "Treated", "Treated")); expect_match(a3$note, "4 samples found")
  expect_match(strym_auto_sample_sheet(modifyList(d, list(data_path = file.path(root, "parts"))))$note, "samples found")
  h5 <- file.path(root, "x.h5"); writeLines("x", h5); expect_match(strym_auto_sample_sheet(modifyList(d, list(data_path = h5)))$note, "cannot list the samples inside this file")
  # DESeq2 reads a counts table: it is grayed out, with the reason, for 10x folders, .h5 files and Seurat objects (and still offered for whole organ data)
  pl <- strym_plan(d); expect_false(pl$applicable[pl$id == "deseq2"]); expect_match(pl$note[pl$id == "deseq2"], "reads a counts table.*not a 10x folder, .h5 file or Seurat object")
  expect_false(pl$default[pl$id == "deseq2"]); expect_false(strym_plan(modifyList(d, list(data_kind = "seurat", de_methods = "deseq2")))$default[strym_plan(d)$id == "deseq2"])
  bulk <- d; bulk$data_type <- "bulk"; bulk$data_kind <- "counts"; expect_true(strym_plan(bulk)$applicable[strym_plan(bulk)$id == "deseq2"])
})
