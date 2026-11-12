# Single cell input handling that must work for any study, not just one dataset: many barcode naming styles, matrices stored cells x genes,
# older / unusual Seurat objects, and file types the step cannot read. Quick checks only (no clustering run).
skip_if_not_installed("Seurat")
source(file.path(here, "analysis", "seurat_clustering.R"))
sc_ctx <- function(path, out, ..., sample_sheet = "") {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  d <- strym_default_design(); d$data_type <- "single_cell"; d$data_kind <- "10x"; d$data_path <- path; d$project_dir <- dirname(out)
  a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]
  logs <- new.env(); logs$x <- character()
  list(design = strym_normalize_design(d), data_path = path, sample_sheet = sample_sheet, out_dir = out, project_dir = dirname(out), app_dir = here, prereq = list(),
       log = function(...) logs$x <- c(logs$x, paste0(...)), logs = logs)
}
B <- "AAACCTGAGAAACCAT"

test_that("sample labels are found in many barcode styles, in front of the barcode or behind it, in any case", {
  expect_equal(strym_barcode_labels(c(paste0("SAMN15315800-", B), paste0("Ctrl_", B), paste0("Sample1#", B, "-1"), paste0("S.2:", B), paste0("lib 3|", B))),
               c("SAMN15315800", "Ctrl", "Sample1", "S.2", "lib 3"))
  expect_equal(strym_barcode_labels(c(paste0(B, "-1_SampleA"), paste0(B, "_B2"), paste0(tolower(B), ".B2"), paste0(B, "-3#Rep 4"), "ACGTACGTAC-S1")),
               c("SampleA", "B2", "B2", "Rep 4", "S1"))
  expect_equal(strym_barcode_labels(c(paste0(B, "-2"), paste0(B, "-1"), B, "ct1-AAAC0001", "", NA, "CELL_12", "gene")), c("2", "1", NA, NA, NA, NA, NA, NA))   # no label, or not a barcode
  expect_equal(strym_split_labels(c(paste0(B, "-1_A"), paste0(B, "-1_A"), paste0(B, "-1_B"), paste0(B, "-1_B"), rep(paste0(B, "-1_A"), 20), rep(paste0(B, "-1_B"), 20))), c("A", "B"))
  expect_equal(strym_split_labels(rep(paste0("Ctrl#", B), 30)), character())                                   # one sample only: nothing to split
  x <- c(paste0("Ctrl_", B), paste0(B, "-1_Ctrl"), paste0(B, "-2"), paste0("Trt_", B))                           # the label comes off whichever side it was on
  expect_equal(strym_strip_label(x, "Ctrl"), c(B, paste0(B, "-1"), paste0(B, "-2"), paste0("Trt_", B))); expect_equal(strym_strip_label(paste0(B, "-2"), "2"), B)
})

test_that("labels in any style split a matrix by the sample sheet", {
  set.seed(1); genes <- paste0("G", 1:30)
  code <- function(n) vapply(seq_len(n), function(i) paste(sample(c("A", "C", "G", "T"), 16, TRUE), collapse = ""), "")                 # 16 letters, as in a real barcode
  cells <- c(paste0(code(30), "-1_Wt"), paste0(code(30), "-1_Mut"))                                                                    # the label comes AFTER the barcode here
  m <- Matrix::Matrix(matrix(rpois(30 * 60, 3) + 1, 30, 60), sparse = TRUE, dimnames = list(genes, cells))
  root <- tempfile("sl"); dir.create(root); ss <- file.path(root, "s.csv"); write.csv(data.frame(sample = c("Wt", "Mut"), condition = c("a", "b")), ss, row.names = FALSE)
  logs <- character(); pr <- .sc_prepare(m, "x", ss, strym_normalize_params(list(sc_min_counts = 10)), function(...) logs <<- c(logs, paste0(...)))
  expect_equal(names(pr$mats), c("Wt", "Mut")); expect_equal(vapply(pr$mats, ncol, 1L), c(Wt = 30, Mut = 30)); expect_true(all(nchar(colnames(pr$mats$Wt)) == 18))   # 16 letters + "-1"
  expect_true(any(grepl("x holds 2 samples; all are analyzed separately", logs)))
})

test_that("a matrix stored cells x genes is turned around; one that fits neither way is refused", {
  root <- tempfile("ct"); dir.create(root); set.seed(3); genes <- paste0("G", 1:12); cells <- sprintf("%s%02d", substr(B, 1, 14), 1:7)
  m <- Matrix::Matrix(matrix(rpois(12 * 7, 2) + 1, 12, 7), sparse = TRUE)
  wr <- function(d, mat, ft, bc) { dir.create(d); Matrix::writeMM(mat, file.path(d, "matrix.mtx")); writeLines(sprintf("%s\t%s\tGene Expression", ft, ft), file.path(d, "features.tsv")); writeLines(bc, file.path(d, "barcodes.tsv")) }
  wr(file.path(root, "ok"), m, genes, cells); wr(file.path(root, "flipped"), Matrix::t(m), genes, cells); wr(file.path(root, "bad"), m, genes[1:10], cells)
  logs <- character(); lg <- function(...) logs <<- c(logs, paste0(...))
  a <- .read_10x(file.path(root, "ok"), lg); b <- .read_10x(file.path(root, "flipped"), lg)
  expect_equal(dim(a), c(12, 7)); expect_equal(as.matrix(a), as.matrix(b)); expect_equal(rownames(b), genes); expect_equal(colnames(b), cells)
  expect_length(logs, 1); expect_match(logs, "stored as cells x genes \\(12 x 7\\): it was turned around")
  expect_error(.read_10x(file.path(root, "bad")), "do not agree: matrix is 12 x 7, features 10, barcodes 7")
})

test_that("unusual input files get a plain message that says what to do", {
  root <- tempfile("fx"); dir.create(root); out <- file.path(root, "o")
  for (f in c("x.h5ad", "x.loom", "x.csv.gz", "x.xyz")) writeLines("x", file.path(root, f))
  e <- function(f) tryCatch(strym_step_seurat_clustering(sc_ctx(file.path(root, f), out)), error = function(e) conditionMessage(e))
  expect_match(e("x.h5ad"), "cannot read “x.h5ad” as single cell data.*10x folder.*10x .h5 file.*Seurat .rds.*write_10x_mtx")
  expect_match(e("x.loom"), "SeuratDisk or zellkonverter"); expect_match(e("x.csv.gz"), "write10xCounts"); expect_match(e("x.xyz"), "cannot read “x.xyz”")
  m <- Matrix::Matrix(matrix(rpois(20 * 30, 3), 20, 30, dimnames = list(paste0("g", 1:20), paste0("c", 1:30))), sparse = TRUE)
  atac <- Seurat::CreateSeuratObject(counts = m, assay = "ATAC"); saveRDS(atac, file.path(root, "atac.rds"))
  expect_match(e("atac.rds"), "no “RNA” assay \\(it has: ATAC\\)")
  saveRDS(1:3, file.path(root, "no.rds")); expect_match(e("no.rds"), "not a Seurat object")
})
