# An expression table (FPKM) plus Cuffdiff comparison files made elsewhere: the entropy steps must work whatever the naming / layout.
mk <- function(genotypes = c("A", "B"), tps = c("T1", "T2", "T3"), ngenes = 60, seed = 1) {
  set.seed(seed); dir <- tempfile("diffin_"); dir.create(dir)
  ids <- sprintf("XLOC_%06d", seq_len(ngenes)); nm <- sprintf("Gene%03d", seq_len(ngenes))
  grid <- expand.grid(cond = c("Control", "Heated"), tp = tps, g = genotypes, stringsAsFactors = FALSE)
  smp <- paste(grid$g, grid$cond, grid$tp, sep = "_")
  M <- matrix(round(rexp(ngenes * length(smp), 1/5), 4), ngenes, dimnames = list(ids, smp))
  tab <- data.frame(tracking_id = ids, gene_short_name = nm, M, check.names = FALSE); write.csv(tab, file.path(dir, "fpkm.csv"), row.names = FALSE, quote = FALSE)
  sheet <- data.frame(sample = smp, condition = paste(grid$g, grid$cond, sep = "_"), time_point = grid$tp, replicate = 1)
  write.csv(sheet, file.path(dir, "sheet.csv"), row.names = FALSE)
  one <- function(g, tp, label_g = g) {
    c1 <- paste(g, "Control", tp, sep = "_"); c2 <- paste(g, "Heated", tp, sep = "_")
    data.frame(test_id = ids, gene_id = ids, gene = nm, locus = "x", sample_1 = paste0(tp, "_", label_g, "_Control"), sample_2 = paste0(tp, "_", label_g, "_Heated"),
               status = "OK", value_1 = M[, c1], value_2 = M[, c2], `log2(fold_change)` = log2((M[, c2] + 1e-3) / (M[, c1] + 1e-3)),
               test_stat = 0, p_value = 0.01, q_value = ifelse(seq_len(ngenes) %% 3 == 0, 0.001, 0.9), significant = "no", check.names = FALSE)
  }
  list(dir = dir, sheet = file.path(dir, "sheet.csv"), table = file.path(dir, "fpkm.csv"), one = one, ids = ids)
}
wr <- function(d, f) { dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE); write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE) }

test_that("per-timepoint files in folders, a zip, and labels that differ from the sheet are all resolved", {
  x <- mk(); root <- file.path(x$dir, "diffs")
  for (tp in c("T1", "T2", "T3")) wr(x$one("A", tp), file.path(root, paste0(tp, "_diff_out"), "gene_exp.diff"))        # folder per time point
  zdir <- file.path(x$dir, "zsrc"); for (tp in c("T1", "T2", "T3")) wr(x$one("B", tp, label_g = "Renamed"), file.path(zdir, paste0(tp, "_diff_out"), "gene_exp.diff"))  # genotype called something else
  zf <- zdir                                                                                                   # as a folder; as a .zip too when a zip tool exists
  if (nzchar(Sys.which("zip"))) { zf <- file.path(x$dir, "B.zip"); old <- setwd(zdir); utils::zip(zf, list.files(".", recursive = TRUE), flags = "-q"); setwd(old) }
  smp <- read.csv(x$sheet, stringsAsFactors = FALSE); smp$group <- paste(smp$condition, smp$time_point, sep = "_")
  ld <- quiet(strym_load_counts(x$table, x$sheet, diff_path = c(root, zf)))
  expect_equal(ld$unit, "FPKM/TPM"); expect_false(is.null(ld$diff))
  expect_setequal(unique(c(ld$diff$sample_1, ld$diff$sample_2)), smp$group)           # every label became a sheet group (B by its values)
  expect_equal(unname(ld$gene_names[x$ids[1]]), "Gene001")                              # readable names kept from the table
  tc <- quiet(strym_timecourse(ld, strym_normalize_params(list()), c("T1", "T2", "T3")))
  expect_equal(nrow(tc$pairs), 2); expect_true(any(tc$qval[[1]] < 0.05))
})

test_that("zips inside folders inside zips, and a zip next to its own unpacked copy, are all found once", {
  skip_if(!nzchar(Sys.which("zip")), "no zip tool")
  x <- mk(genotypes = "A"); base <- file.path(x$dir, "nest"); src <- file.path(x$dir, "src")
  for (tp in c("T1", "T2", "T3")) {
    wr(x$one("A", tp), file.path(src, tp, paste0(tp, "_diff_out"), "gene_exp.diff"))
    old <- setwd(file.path(src, tp)); dir.create(file.path(base, tp), recursive = TRUE)
    utils::zip(file.path(base, tp, paste0(tp, ".zip")), paste0(tp, "_diff_out"), flags = "-qr"); setwd(old)
  }
  expect_equal(nrow(strym_diff_files(base)), 3)                                      # a zip per time-point folder
  old <- setwd(base); utils::zip(file.path(x$dir, "outer.zip"), ".", flags = "-qr"); setwd(old)
  expect_equal(nrow(strym_diff_files(file.path(x$dir, "outer.zip"))), 3)             # ... inside one more zip
  ld <- quiet(strym_load_counts(x$table, x$sheet, diff_path = c(base, file.path(x$dir, "outer.zip"))))
  expect_equal(nrow(ld$diff), 3 * 60)                                                # the same results twice count once
})

test_that("a design with Cuffdiff results enables the entropy steps, and without them they stay blocked", {
  d <- list(data_type = "bulk", data_kind = "counts", n_conditions = 4, n_timepoints = 3, expr_unit = "normalized")    # FPKM/TPM: no q-values can be calculated
  expect_false(strym_plan(d)$applicable[strym_plan(d)$id == "modified_shannon_entropy"])
  d$diff_path <- "somewhere"
  expect_true(strym_plan(d)$applicable[strym_plan(d)$id == "modified_shannon_entropy"])
})

test_that("time points are read from the end of sample names when none were typed", {
  it <- .infer_time_suffix(c("A_Control_T1", "A_Heated_T1", "A_Control_T2", "A_Heated_T2"))
  expect_equal(it$tp, c("T1", "T1", "T2", "T2")); expect_equal(it$cond, c("A_Control", "A_Heated", "A_Control", "A_Heated"))
  it2 <- .infer_time_suffix(c("ctrl_0h_rep1", "ctrl_0h_rep2", "heat_4h_rep1", "heat_4h_rep2"))
  expect_equal(it2$tp, c("0h", "0h", "4h", "4h")); expect_equal(it2$cond, c("ctrl", "ctrl", "heat", "heat"))
  expect_null(.infer_time_suffix(c("ctrl_rep1", "ctrl_rep2", "heat_rep1", "heat_rep2")))
})
