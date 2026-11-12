write_inputs <- function(d, n = 14, tp = 6, seed = 11) {
  set.seed(seed); genes <- sprintf("G%02d", seq_len(n))
  tc <- matrix(exp(rnorm(n * tp, 3, .5)), n, tp)
  write.csv(data.frame(gene = genes, tc), file.path(d, "tc.csv"), row.names = FALSE)
  write.csv(data.frame(gene = genes), file.path(d, "genes.csv"), row.names = FALSE)
  write.csv(data.frame(gene = genes, matrix(rnorm(n * 8), n, 8)), file.path(d, "expr.csv"), row.names = FALSE)
  genes
}

test_that("GENIST with no transcription factor among the genes returns only isolated nodes", {
  d <- tempfile("e"); dir.create(d); genes <- write_inputs(d)
  write.csv(data.frame(TF = c("X1", "X2")), file.path(d, "tf.csv"), row.names = FALSE)
  quiet(genist_run(file.path(d, "genes.csv"), file.path(d, "tc.csv"), tf_file = file.path(d, "tf.csv"),
                   output_folder = d, output_file = "net.txt"))
  out <- read.delim(file.path(d, "net.txt"))
  expect_equal(nrow(out), length(genes)); expect_true(all(out$width == 0)); expect_true(all(out$source == out$target))
})

test_that("GENIST skips genes missing from the data and says so", {
  d <- tempfile("e"); dir.create(d); genes <- write_inputs(d)
  write.csv(data.frame(gene = c(genes, "NOT_THERE")), file.path(d, "genes.csv"), row.names = FALSE)
  msgs <- testthat::capture_messages(genist_run(file.path(d, "genes.csv"), file.path(d, "tc.csv"),
                                                output_folder = d, output_file = "net.txt"))
  expect_true(any(grepl("missing from the time course", msgs)))
  expect_false("NOT_THERE" %in% read.delim(file.path(d, "net.txt"))$source)
})

test_that("GENIST needs at least 3 time points", {
  d <- tempfile("e"); dir.create(d); write_inputs(d, tp = 2)
  expect_error(quiet(genist_run(file.path(d, "genes.csv"), file.path(d, "tc.csv"), output_folder = d)), "at least 3 time points")
})

test_that("GENIST uses symbols from a name file without a header row, including its first row", {
  d <- tempfile("e"); dir.create(d); genes <- write_inputs(d, n = 8)
  write.table(data.frame(genes, paste0("SYM_", genes)), file.path(d, "sym.csv"), sep = ",", row.names = FALSE,
              col.names = FALSE, quote = FALSE)                              # no header: first row is data
  quiet(genist_run(file.path(d, "genes.csv"), file.path(d, "tc.csv"), symbol_file = file.path(d, "sym.csv"),
                   output_folder = d, output_file = "net.txt"))
  out <- read.delim(file.path(d, "net.txt"))
  expect_true(all(grepl("^SYM_", c(out$source, out$target))))
})

test_that("RTP-STAR without any regulator writes an empty result instead of failing", {
  d <- tempfile("e"); dir.create(d); write_inputs(d)
  write.csv(data.frame(TF = "X1"), file.path(d, "tf.csv"), row.names = FALSE)
  quiet(rtpstar_run(file.path(d, "genes.csv"), file.path(d, "expr.csv"), d, "res.txt", tf_file = file.path(d, "tf.csv"), ntrees = 50))
  expect_equal(length(readLines(file.path(d, "res.txt"))), 0)
})

test_that("RTP-STAR rejects a time course that is too short", {
  d <- tempfile("e"); dir.create(d); write_inputs(d, tp = 2)
  expect_error(quiet(rtpstar_run(file.path(d, "genes.csv"), file.path(d, "expr.csv"), d, time_course_file = file.path(d, "tc.csv"))),
               "at least 3 time points")
})

test_that("RTP-STAR matches genes by name, not by row position, across files", {
  d <- tempfile("e"); dir.create(d); genes <- write_inputs(d, n = 12, tp = 5)
  tc <- read.csv(file.path(d, "tc.csv")); write.csv(tc[rev(seq_len(nrow(tc))), ], file.path(d, "tc_rev.csv"), row.names = FALSE)  # reversed row order
  for (f in c("tc.csv", "tc_rev.csv"))
    quiet(rtpstar_run(file.path(d, "genes.csv"), file.path(d, "expr.csv"), d, paste0("r_", f), time_course_file = file.path(d, f),
                      ntrees = 200, rf_seed = 3))
  expect_equal(readLines(file.path(d, "r_tc.csv")), readLines(file.path(d, "r_tc_rev.csv")))
})

test_that("command-line runners print DONE on success and ERROR (exit 1) on failure", {
  d <- tempfile("e"); dir.create(d); write_inputs(d)
  runner <- file.path(here, "run_genist_cli.R")
  cfg_ok <- file.path(d, "ok.json")
  jsonlite::write_json(list(genes_file = file.path(d, "genes.csv"), time_course_file = file.path(d, "tc.csv"),
                            output_folder = d, output_file = "net.txt"), cfg_ok, auto_unbox = TRUE)
  res <- suppressWarnings(system2("Rscript", c(shQuote(runner), shQuote(cfg_ok)), stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("^DONE: ", res))); expect_null(attr(res, "status"))
  cfg_bad <- file.path(d, "bad.json")
  jsonlite::write_json(list(genes_file = file.path(d, "missing.csv"), time_course_file = file.path(d, "tc.csv"),
                            output_folder = d), cfg_bad, auto_unbox = TRUE)
  res <- suppressWarnings(system2("Rscript", c(shQuote(runner), shQuote(cfg_bad)), stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("^ERROR: .*cannot be found", res))); expect_equal(attr(res, "status"), 1)
})
