# q-values calculated from raw counts (DESeq2) so the time-course steps work without a Cuffdiff run; time points a condition lacks are left out.
test_that("DESeq2-derived comparisons have the shape of Cuffdiff's and find the true genes", {
  x <- mk_counts(); ld <- quiet(strym_load_counts(x$table, x$sheet)); expect_null(ld$diff); expect_false(is.null(ld$counts))
  p <- strym_normalize_params(list()); cf <- file.path(x$dir, "cache", "d.rds")
  de <- quiet(strym_derive_diff(ld$counts, ld$samples, p, c("0h", "1h", "3h"), cache_file = cf))
  expect_true(all(c("test_id", "sample_1", "sample_2", "status", "q_value", "log2fc") %in% names(de)))
  expect_setequal(unique(de$sample_2), c("heat_0h", "heat_1h", "heat_3h")); expect_setequal(unique(de$sample_1), c("control_0h", "control_1h", "control_3h"))
  late <- de[de$sample_2 == "heat_3h", ]; expect_true(all(late$q_value[late$test_id %in% sprintf("G%03d", 1:20)] < 0.05))
  expect_true(mean(late$q_value[!late$test_id %in% sprintf("G%03d", 1:20)] < 0.05) < 0.1)             # almost no false positives
  expect_true(file.exists(cf))
  msgs <- character(); de2 <- strym_derive_diff(ld$counts, ld$samples, p, c("0h", "1h", "3h"), cache_file = cf, log = function(m) msgs <<- c(msgs, m))
  expect_match(paste(msgs, collapse = " "), "reusing"); expect_equal(nrow(de2), nrow(de))
})

test_that("the entropy time course runs from a counts table, and leaves out a time point a condition lacks", {
  x <- mk_counts(tps = c("0h", "1h", "3h", "6h"), drop = c("heat", "3h"))
  ld <- quiet(strym_load_counts(x$table, x$sheet)); p <- strym_normalize_params(list(min_fpkm = 1))
  ld$diff <- quiet(strym_derive_diff(ld$counts, ld$samples, p, c("0h", "1h", "3h", "6h")))
  msgs <- character(); tc <- strym_timecourse(ld, p, c("0h", "1h", "3h", "6h"), log = function(m) msgs <<- c(msgs, m))
  expect_equal(tc$tps, c("0h", "1h", "6h")); expect_match(paste(msgs, collapse = " "), "Left out time point 3h")
  expect_error(strym_timecourse(ld, p, c("0h", "1h", "3h", "6h")[c(1, 3)]) , NA)                       # fewer names typed than exist: still fine
})

test_that("a comparison that cannot be made is reported, not invented", {
  x <- mk_counts(reps = 1); ld <- quiet(strym_load_counts(x$table, x$sheet)); p <- strym_normalize_params(list())
  expect_error(quiet(strym_derive_diff(ld$counts, ld$samples, p, c("0h", "1h", "3h"))), "3 are needed|No comparison")
  expect_match(strym_derive_possible(list(data_kind = "counts", expr_unit = "normalized", sample_sheet = "")), "normalized")
  expect_equal(strym_derive_possible(list(data_kind = "counts", expr_unit = "auto", sample_sheet = x$sheet)), "No condition has replicates (two or more samples per condition and time point), which the significance tests need.")
})

test_that("the plan enables the entropy steps for a counts table with replicates and explains when it cannot", {
  x <- mk_counts()
  d <- list(data_type = "bulk", data_kind = "counts", n_conditions = 2, n_timepoints = 3, sample_sheet = x$sheet, data_path = x$table)
  pl <- strym_plan(d); expect_true(pl$applicable[pl$id == "modified_shannon_entropy"]); expect_true(pl$applicable[pl$id == "genist_network"])
  d$expr_unit <- "normalized"; pl <- strym_plan(d); expect_false(pl$applicable[pl$id == "modified_shannon_entropy"])
  expect_match(pl$note[pl$id == "modified_shannon_entropy"], "normalized")
})
