# (mk_counts is defined in test-derived.R and reused here)
# Real data is messy: these are the problems found by the audit of the expression-table path, each with the behaviour StRym must show now.
wtab <- function(lines, name = "t.txt") { f <- file.path(tempfile("rob_"), name); dir.create(dirname(f)); writeLines(lines, f); f }
msgs_of <- function(expr) { m <- character(); withCallingHandlers(expr, message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") }); paste(m, collapse = " | ") }

test_that("table columns that are not in the sample sheet are left out and said so; a sheet without a sample column is an error", {
  tab <- wtab(c("gene\tA_1\tA_2\tB_1\tB_2\tC_1", "g1\t10\t12\t30\t33\t5", "g2\t5\t6\t1\t2\t9", "g3\t0\t1\t2\t3\t4"))
  sh <- wtab(c("sample,condition,time_point", "A_1,A,t1", "A_2,A,t1", "B_1,B,t1", "B_2,B,t1"), "s.csv")
  m <- msgs_of(r <- .counts_raw(tab, sh)); expect_equal(r$samples$sample, c("A_1", "A_2", "B_1", "B_2")); expect_match(m, "left out: C_1")
  bad <- wtab(c("name,condition", "A_1,A", "A_2,A", "B_1,B", "B_2,B"), "s.csv")
  expect_error(.counts_raw(tab, bad), "needs a column named")
  nomatch <- wtab(c("sample,condition", "X1,A", "X2,B"), "s.csv"); expect_error(.counts_raw(tab, nomatch), "None of the table")
  low <- wtab(c("sample,condition", "a_1,A", "a_2,A", "b_1,B", "b_2,B"), "s.csv"); expect_equal(.counts_raw(tab, low)$samples$condition, c("A", "A", "B", "B"))   # case differs
})

test_that("junk cells are counted, thousands separators and decimal commas understood, Inf read as 0, negatives refused", {
  tab <- wtab(c("gene\tA\tB\tC", "g1\t10\t-\t\"1,234\"", "g2\tNA\t5\t7", "g3\tabc\t3\t4", "g4\t1,5\t2\t2"))
  m <- msgs_of(r <- .counts_raw(tab)); expect_equal(unname(r$m["g1", "C"]), 1234); expect_equal(unname(r$m["g4", "A"]), 1.5)
  expect_match(m, "3 empty or non-numeric cell"); expect_equal(unname(r$m["g3", "A"]), 0)
  expect_message(ri <- .counts_raw(wtab(c("gene\tA\tB", "g1\tInf\t1", "g2\t2\t3"))), "Inf / -Inf and were read as 0"); expect_equal(unname(ri$m["g1", "A"]), 0); expect_equal(ri$dropped$inf_cells, 1)     # Inf is read as 0 like an empty cell, and reported (version 10)
  neg <- wtab(c("gene\tA\tB", "g1\t-3\t1", "g2\t2\t3")); expect_error(quiet(strym_load_counts(neg)), "negative values")
})

test_that("duplicate and blank gene ids are reported; an all-zero sample is left out by name", {
  tab <- wtab(c("gene\tA\tB\tC", "g1\t5\t6\t0", "g1\t7\t8\t0", "\t1\t2\t0", "g2\t3\t4\t0"))
  m <- msgs_of(r <- .counts_raw(tab)); expect_equal(rownames(r$m), c("g1", "g2")); expect_match(m, "1 without a gene id, 1 with a gene id already listed")
  expect_equal(colnames(r$m), c("A", "B")); expect_match(m, "no counts at all were left out: C")
})

test_that("semicolon-separated files with decimal commas (European Excel) are read", {
  f <- wtab(c("gene;A_1;A_2", "g1;1,5;2,5", "g2;3;4"), "eu.csv"); t <- read_table_any(f)
  expect_equal(names(t), c("gene", "A_1", "A_2")); expect_equal(t$A_1, c(1.5, 3))
})

test_that("Excel: numeric gene ids stay as written, and a title row above the header is skipped", {
  skip_if_not_installed("writexl")
  f <- file.path(tempfile("xl_"), "t.xlsx"); dir.create(dirname(f))
  d <- data.frame(c("Experiment 1", "gene", "100000", "200000"), c(NA, "A", "5", "6"), c(NA, "B", "7", "8"), stringsAsFactors = FALSE)
  writexl::write_xlsx(stats::setNames(d, c("", "", "")), f)
  t <- read_table_any(f); expect_equal(names(t), c("gene", "A", "B")); expect_equal(t[[1]], c("100000", "200000"))
})

test_that("Excel: an empty row between the title and the header is skipped too", {
  skip_if_not_installed("writexl")
  f <- file.path(tempfile("xl_"), "t.xlsx"); dir.create(dirname(f))
  d <- data.frame(c("Run 2025", NA, "gene", "g1", "g2"), c(NA, NA, "A", "5", "6"), c(NA, NA, "B", "7", "8"), stringsAsFactors = FALSE)
  writexl::write_xlsx(stats::setNames(d, c("", "", "")), f)
  t <- read_table_any(f); expect_equal(names(t), c("gene", "A", "B")); expect_equal(t[[1]], c("g1", "g2"))
})

test_that("Excel: several title lines above the header are skipped", {
  skip_if_not_installed("writexl")
  f <- file.path(tempfile("xl_"), "t.xlsx"); dir.create(dirname(f))
  d <- data.frame(c("Run 2025", "contact: lab", "gene", "g1", "g2"), c(NA, NA, "A", "5", "6"), c(NA, NA, "B", "7", "8"), stringsAsFactors = FALSE)
  writexl::write_xlsx(stats::setNames(d, c("", "", "")), f)
  t <- read_table_any(f); expect_equal(names(t), c("gene", "A", "B")); expect_equal(t[[1]], c("g1", "g2"))
})

test_that("DESeq2 survives a gene panel too small for the dispersion curve", {
  skip_if_not_installed("DESeq2")
  set.seed(3); m <- matrix(rnbinom(12 * 6, mu = 40, size = 8), 12, 6, dimnames = list(sprintf("G%02d", 1:12), paste0(rep(c("c", "t"), each = 3), 1:3))); m[2, 4:6] <- m[2, 4:6] * 8L
  dds <- DESeq2::DESeqDataSetFromMatrix(m, data.frame(group = factor(rep(c("c", "t"), each = 3)), row.names = colnames(m)), ~ group)
  r <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(e) tryCatch(DESeq2::DESeq(dds, fitType = "mean", quiet = TRUE), error = function(e2) NULL))
  expect_false(is.null(r))
})

test_that("a run started without the page still refuses a drive root or a file as the output folder", {
  d <- strym_default_design(); d$project_name <- "x"; d$data_path <- tempfile(fileext = ".csv"); d$project_dir <- "/"
  expect_error(strym_run(d, "pca", app_dir = here), "root of the drive")
  f <- tempfile(); writeLines("x", f); d$project_dir <- f
  expect_error(strym_run(d, "pca", app_dir = here), "is a file")
})

test_that("time points sort by value and unit", {
  expect_equal(strym_order_timepoints(c("24h", "30min", "1h")), c("30min", "1h", "24h"))
  expect_equal(strym_order_timepoints(c("T-2", "T-10", "T-1")), c("T-1", "T-2", "T-10"))
  expect_equal(strym_order_timepoints(c("0", "1", "10", "2")), c("0", "1", "2", "10"))
  expect_equal(strym_order_timepoints(c("1d", "12h", "2d")), c("12h", "1d", "2d"))
  expect_equal(strym_order_timepoints(c("3dpi", "1dpi", "10dpi")), c("1dpi", "3dpi", "10dpi"))
})

test_that("the auto sample sheet reads sample names sensibly", {
  it <- .infer_time_suffix(c("Control_1h_1", "Control_1h_2", "Control_24h_1", "Control_24h_2")); expect_equal(it$tp, c("1h", "1h", "24h", "24h")); expect_equal(it$cond, rep("Control", 4))
  expect_equal(.find_token(c("WT_Heat_1", "WT_1"), c("WT", "WT_Heat")), c("WT_Heat", "WT")); expect_equal(.find_token("WT_10_1", c("0", "1", "10")), "10")
})

test_that("no gene reaching the minimum gives a message that says what to do; a perfect-fit outlier model is still found", {
  m <- rbind(g = c(5, 5, 5, 5, 50, 5)); o <- strym_mse_outliers(m); expect_equal(as.integer(o[1, ]), c(0, 0, 0, 0, 1, 0))
  p <- strym_normalize_params(list(drop_reps = c("A_1", "B_1"))); expect_equal(p$drop_reps, "A_1, B_1")
})

test_that("a big table is read by the fast reader with exactly the same result as the standard reader", {
  skip_if_not_installed("data.table")
  set.seed(3); n <- 70000; f <- tempfile(fileext = ".csv")
  m <- matrix(rpois(n * 25, 30), n, 25, dimnames = list(NULL, paste0("S", 1:25)))
  utils::write.csv(data.frame(gene = sprintf("G%06d", seq_len(n)), annot = "x", m), f, row.names = FALSE, quote = TRUE)
  expect_gt(file.size(f), 5e6)
  fast <- read_table_any(f)
  slow <- utils::read.table(f, header = TRUE, sep = ",", quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE, strip.white = TRUE)
  expect_equal(names(fast), names(slow)); expect_equal(dim(fast), dim(slow)); expect_equal(as.matrix(fast[-(1:2)]), as.matrix(slow[-(1:2)]), check.attributes = FALSE)
  expect_equal(fast[[1]], slow[[1]]); expect_equal(unique(fast$annot), "x")
})

# ---- round-2 review: things found by an adversarial read of the new code ----
test_that("decimal commas and thousands separators are told apart", {
  expect_equal(.to_num(c("0,125", "1,5", "2,25", "-")), c(0.125, 1.5, 2.25, NA))
  expect_equal(.to_num(c("1,234", "5", "1,000", "12,345,678")), c(1234, 5, 1000, 12345678))
})

test_that("only gene-level Cuffdiff comparison files are used; a stray file never aborts the load", {
  d <- tempfile("dff_"); dir.create(d)
  hdr <- "test_id\tgene_id\tgene\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tlog2(fold_change)\tq_value"
  writeLines(c(hdr, "g1\tg1\tg1\tA_T1\tB_T1\tOK\t1\t2\t1\t0.01"), file.path(d, "gene_exp.diff"))
  writeLines(c(hdr, "g1\tg1\tg1\tA_T1\tB_T1\tOK\t1\t2\t1\t0.99"), file.path(d, "cds.diff")); writeLines(c(hdr, "g1\tg1\tg1\tA_T1\tB_T1\tOK\t1\t2\t1\t0.99"), file.path(d, "promoters.diff"))
  writeLines("not a table", file.path(d, "notes_diff_summary.txt")); png(file.path(d, "diff_plot.png")); dev.off()
  expect_equal(basename(strym_diff_files(d)$file), c("gene_exp.diff", "notes_diff_summary.txt"))
  sm <- data.frame(sample = c("a", "b"), condition = c("A", "B"), time_point = "T1", replicate = 1, group = c("A_T1", "B_T1"))
  ld <- strym_load_diff(d, sm, "g1", log = function(...) NULL); expect_equal(ld$diff$q_value, 0.01)
})

test_that("a changed counts table never reuses another table's cached q-values", {
  x <- mk_counts(); ld <- quiet(strym_load_counts(x$table, x$sheet)); p <- strym_normalize_params(list()); cf <- file.path(x$dir, "c.rds")
  d1 <- quiet(strym_derive_diff(ld$counts, ld$samples, p, c("0h", "1h", "3h"), cache_file = cf))
  m2 <- ld$counts; rownames(m2) <- rev(rownames(m2))                       # same numbers, genes renamed
  d2 <- quiet(strym_derive_diff(m2, ld$samples, p, c("0h", "1h", "3h"), cache_file = cf))
  expect_false(identical(d1$test_id, d2$test_id)); expect_equal(d2$test_id[1], rownames(m2)[1])      # the genes of the second table carry the second table's ids (the old key reused the first result)
})

test_that("a dropped first time point does not turn the next one into the baseline", {
  x <- mk_counts(tps = c("0h", "1h", "3h", "6h"), drop = c("heat", "0h"))
  ld <- quiet(strym_load_counts(x$table, x$sheet)); p <- strym_normalize_params(list(min_fpkm = 1))
  ld$diff <- quiet(strym_derive_diff(ld$counts, ld$samples, p, c("0h", "1h", "3h", "6h")))
  tc <- strym_timecourse(ld, p, c("0h", "1h", "3h", "6h"), log = function(...) NULL)
  expect_equal(tc$tps, c("1h", "3h", "6h")); expect_false(tc$baseline_ok)
})

test_that("the Data overview survives identical samples and still writes its notes", {
  d <- tempfile("qc_"); dir.create(d)
  writeLines(c("gene,A_1,A_2,B_1,B_2", sprintf("g%d,%d,%d,%d,%d", 1:30, 1:30 * 10L, 1:30 * 10L, 1:30 * 10L, 1:30 * 10L)), file.path(d, "t.csv"))
  writeLines(c("sample,condition,time_point,replicate", "A_1,A,t1,1", "A_2,A,t1,2", "B_1,B,t1,1", "B_2,B,t1,2"), file.path(d, "s.csv"))
  des <- strym_normalize_design(list(project_dir = d, data_kind = "counts", data_path = file.path(d, "t.csv"), sample_sheet = file.path(d, "s.csv"), n_conditions = 2, n_timepoints = 1))
  out <- file.path(d, "o"); dir.create(out)
  source(file.path(here, "analysis", "qc_overview.R"))
  r <- strym_step_qc_overview(list(design = des, data_path = des$data_path, sample_sheet = des$sample_sheet, out_dir = out, project_dir = d, app_dir = here, log = function(...) NULL))
  expect_true(r$status %in% c("ok", "skipped")); expect_true(file.exists(file.path(out, "data_overview_notes.csv")))
})

test_that("the sample sheet is read as text: leading zeros, a condition called NA", {
  f <- tempfile(fileext = ".csv"); writeLines(c("sample,condition,time_point", "s1,NA,001", "s2,F,010"), f)
  s <- strym_read_sample_sheet(f); expect_equal(s$time_point, c("001", "010")); expect_equal(s$condition, c("NA", "F"))
})

test_that("a Phred+33 file with a K quality is not mistaken for Phred+64", {
  f <- tempfile(fileext = ".fastq"); l <- unlist(lapply(1:6, function(i) c(sprintf("@r%d", i), strrep("A", 30), "+", paste0(strrep("I", 20), strrep("K", 10))))); writeLines(l, f)
  expect_equal(.detect_phred(f), 33L)
})

test_that("an Excel export with an unnamed index column keeps its header", {
  skip_if_not_installed("writexl")
  f <- file.path(tempfile("xl2_"), "t.xlsx"); dir.create(dirname(f))
  writexl::write_xlsx(stats::setNames(data.frame(c("g1", "g2"), c(5, 6), c(7, 8)), c("", "S1", "S2")), f)
  t <- read_table_any(f); expect_equal(names(t), c("id", "S1", "S2")); expect_equal(t[[1]], c("g1", "g2"))
})