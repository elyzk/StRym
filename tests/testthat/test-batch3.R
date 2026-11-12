# Batch 3: chosen FASTQ files, the auto-built sample sheet, the run estimate, finding dropped files, threads.
fake_fq <- function(path, n = 3) {
  rec <- unlist(lapply(seq_len(n), function(i) c(sprintf("@r%d", i), "ACGTACGTAC", "+", "IIIIIIIIII")))
  con <- gzfile(path, "wt"); writeLines(rec, con); close(con)
}
reads_dir <- function(root, files) { r <- file.path(root, "reads"); dir.create(r, showWarnings = FALSE); for (f in files) fake_fq(file.path(r, f)); r }
d0 <- function(...) { d <- strym_default_design(); a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]; d }

# ---- chosen FASTQ files ----
test_that("only the chosen FASTQ files are used, and a missing chosen file is reported", {
  root <- tempfile("b3"); dir.create(root)
  r <- reads_dir(root, c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz", "TRT1_x.fastq.gz", "TRT2_x.fastq.gz", "OTHER_x.fastq.gz"))
  d <- d0(data_kind = "fastq", data_path = r, n_conditions = 2L)
  expect_equal(nrow(strym_quantify_samples(d)), 5)                                       # nothing chosen: every file in the folder
  d$fastq_files <- file.path(r, c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz", "TRT1_x.fastq.gz", "TRT2_x.fastq.gz"))
  s <- strym_quantify_samples(d); expect_equal(basename(s$file), basename(d$fastq_files)); expect_equal(sort(unique(s$group)), c("CTL", "TRT"))
  d$fastq_files <- file.path(r, c("CTL1_x.fastq.gz", "gone.fastq.gz"))
  v <- strym_validate_reads(d); expect_true(any(grepl("Chosen FASTQ file not found: gone.fastq.gz", v$errors)))
  expect_equal(strym_normalize_design(d0(fastq_files = c("a", "b")))$fastq_files, c("a", "b")); expect_length(strym_normalize_design(d0())$fastq_files, 0)
})

test_that("choosing paired-end files keeps both mates and renumbers the replicates", {
  s <- data.frame(group = c("A", "A", "B"), k = c(1L, 2L, 1L), file = c("/x/A1_R1.fastq.gz", "/x/A2_R1.fastq.gz", "/x/B1_R1.fastq.gz"),
                  mate = c("/x/A1_R2.fastq.gz", "/x/A2_R2.fastq.gz", "/x/B1_R2.fastq.gz"), sample = c("A1_R1", "A2_R1", "B1_R1"))
  f <- strym_filter_fastq(s, c("/y/A2_R2.fastq.gz", "/y/B1_R1.fastq.gz"), TRUE)             # choosing either mate selects the pair
  expect_equal(basename(f$file), c("A2_R1.fastq.gz", "B1_R1.fastq.gz")); expect_equal(f$k, c(1L, 1L))
  expect_equal(nrow(strym_filter_fastq(s, character(), TRUE)), 3)
})

# ---- the sample sheet built from the inputs ----
test_that("the sample sheet for raw reads lists the real files and maps groups to the condition names typed", {
  root <- tempfile("b3"); dir.create(root)
  r <- reads_dir(root, c("C2B3_Fe.fastq.gz", "C2B4_Fe.fastq.gz", "C2D2_Fe.fastq.gz", "C2D4_Fe.fastq.gz"))
  d <- d0(data_kind = "fastq", data_path = r, n_conditions = 2L, condition_names = c("B", "D"), timepoint_names = "24hr")
  a <- strym_auto_sample_sheet(d)
  expect_equal(names(a$sheet), c("sample", "condition", "time_point", "replicate", "file"))
  expect_equal(a$sheet$sample, c("C2B3_Fe", "C2B4_Fe", "C2D2_Fe", "C2D4_Fe")); expect_equal(a$sheet$condition, c("B", "B", "D", "D"))
  expect_equal(unique(a$sheet$time_point), "24hr"); expect_equal(a$sheet$replicate, c(1L, 2L, 1L, 2L)); expect_equal(a$sheet$file, paste0(a$sheet$sample, ".fastq.gz"))
  expect_match(a$note, "4 samples found")
  d$fastq_files <- file.path(r, c("C2B3_Fe.fastq.gz", "C2D2_Fe.fastq.gz"))                    # only the chosen files appear
  expect_equal(strym_auto_sample_sheet(d)$sheet$sample, c("C2B3_Fe", "C2D2_Fe"))
  d$sample_sheet <- ""; d$condition_names <- character()                                       # no names typed: the guessed groups
  expect_equal(strym_auto_sample_sheet(d)$sheet$condition, c("C2B", "C2D"))
})

test_that("the sample sheet for a counts table comes from its columns, using the names typed", {
  root <- tempfile("b3"); dir.create(root); demo <- strym_make_demo_data(root)
  d <- d0(data_path = demo$counts, n_conditions = 2L, condition_names = c("control", "heat"), n_timepoints = 3L, timepoint_names = c("0h", "1h", "4h"))
  a <- strym_auto_sample_sheet(d); s <- a$sheet
  expect_equal(nrow(s), 12); expect_false("file" %in% names(s)); expect_equal(as.vector(table(s$condition)), c(6, 6))
  expect_equal(sort(unique(s$time_point)), c("0h", "1h", "4h")); expect_equal(sort(unique(s$replicate)), 1:2)
  expect_equal(s$sample, strym_counts_samples(demo$counts))
  expect_true(all(s$condition[grepl("^heat", s$sample)] == "heat"))
  one <- strym_auto_sample_sheet(modifyList(d, list(n_timepoints = 1L, timepoint_names = character())))$sheet
  expect_equal(unique(one$time_point), "t1")
})

test_that("the sample sheet explains what is missing instead of failing", {
  expect_null(strym_auto_sample_sheet(d0())$sheet); expect_match(strym_auto_sample_sheet(d0())$note, "Choose your input data")
  root <- tempfile("b3"); dir.create(root); e <- file.path(root, "empty"); dir.create(e)
  expect_match(strym_auto_sample_sheet(d0(data_kind = "fastq", data_path = e))$note, "No FASTQ files")
  expect_match(strym_auto_sample_sheet(d0(data_type = "single_cell", data_kind = "10x", data_path = e))$note, "Single cell/single nucleus")
  bad <- file.path(root, "bad.csv"); writeLines("", bad)
  expect_match(strym_auto_sample_sheet(d0(data_path = bad))$note, "could not be read")
})

test_that("a sample sheet whose sample names are not in the counts table is flagged", {
  root <- tempfile("b3"); dir.create(root); demo <- strym_make_demo_data(root)
  ss <- read.csv(demo$sample_sheet); ss$sample[1] <- "typo_sample"; write.csv(ss, file.path(root, "ss.csv"), row.names = FALSE)
  d <- d0(project_name = "p", project_dir = file.path(root, "proj"), data_path = demo$counts, sample_sheet = file.path(root, "ss.csv"), n_conditions = 2L, n_timepoints = 3L)
  v <- strym_validate_design(d); expect_true(any(grepl("1 of the 12 sample names.*not columns of the counts table \\(e.g. typo_sample\\)", v$warnings)))
  expect_false(any(grepl("not columns", strym_validate_design(modifyList(d, list(sample_sheet = demo$sample_sheet)))$warnings)))
})

# ---- the run estimate ----
test_that("durations are written in plain words", {
  expect_equal(strym_fmt_duration(11857), "3 h 18 min"); expect_equal(strym_fmt_duration(1500), "25 min")
  expect_equal(strym_fmt_duration(30), "under a minute"); expect_equal(strym_fmt_duration(3599), "1 h 00 min")
  expect_equal(strym_fmt_duration(75), "1 min"); expect_equal(strym_fmt_duration(75 * 1.1, up = TRUE), "2 min")     # an upper bound never reads lower than the estimate
  expect_equal(strym_fmt_duration(11857 * 1.1, up = TRUE), "3 h 38 min")                                           # 3 h 37.4 min, rounded up
})

test_that("the time estimate reproduces the real 4-sample run (3 h 17 min on 1 thread) and shrinks with more threads", {
  root <- tempfile("b3"); r <- file.path(root, "reads"); dir.create(r, recursive = TRUE)
  sizes <- c(C2B3_Fe = 1482212296, C2B4_Fe = 1232635260, C2D2_Fe = 1778333638, C2D4_Fe = 1634382304)   # the real file sizes (sparse stand-ins)
  for (n in names(sizes)) system2("truncate", c("-s", format(sizes[[n]], scientific = FALSE), file.path(r, paste0(n, ".fastq.gz"))))
  d <- d0(data_kind = "fastq", data_path = r, project_dir = file.path(root, "proj"), threads = 1L, n_conditions = 2L)
  e1 <- strym_estimate_run(d, c("quantify", "deseq2"))
  expect_equal(e1$seconds, 3 * 3600 + 17 * 60 + 15, tolerance = 0.03)                           # real run: 3 h 17 min (+ one quick step)
  expect_equal(e1$reads_m, 79.4, tolerance = 0.03); expect_true(is.finite(e1$need)); expect_true(is.finite(e1$free))
  expect_match(e1$lines, "about 80 million reads from 4 samples with 1 thread")
  e4 <- strym_estimate_run(modifyList(d, list(threads = 4L)), "quantify")
  expect_lt(e4$seconds, e1$seconds / 2); expect_gt(e4$seconds, e1$seconds / 5)
  e0 <- strym_estimate_run(d0(project_dir = file.path(root, "proj")), c("deseq2", "pca"))      # counts-table steps are quick
  expect_lt(e0$seconds, 60); expect_true(is.na(e0$need)); expect_length(e0$lines, 0)
})

# ---- dropped files and threads ----
test_that("dropped files are located by name and size", {
  root <- normalizePath(tempfile("drop")); dir.create(file.path(root, "work", "sub"), recursive = TRUE)
  writeLines("hello", file.path(root, "work", "a.csv")); writeLines("hello world", file.path(root, "b.csv")); dir.create(file.path(root, "work", "folder_x"))
  f <- locate_dropped(c("a.csv", "b.csv", "folder_x"), c(6, 12, 0), c(FALSE, FALSE, TRUE), hints = file.path(root, "work"), max_seconds = 5)
  expect_equal(unname(f), file.path(root, c("work/a.csv", "b.csv", "work/folder_x")))        # the hint folder and its parent are looked at first
  g <- locate_dropped("a.csv", 999, FALSE, hints = file.path(root, "work"), max_seconds = 2)  # same name, other size: not the same file
  expect_true(is.na(g[["a.csv"]]))
  h <- locate_dropped("nope_zzz_9f3a.fastq.gz", 10, FALSE, hints = file.path(root, "work"), max_seconds = 2)
  expect_true(is.na(h[[1]]))
})

test_that("dropped files are found in the usual places on this computer", {
  skip_if_not(dir.exists("/mnt/c/Users"))
  nm <- sprintf("strym_drop_test_%d.txt", sample.int(1e8, 1)); dr <- file.path(path.expand("~"), "strym_drop_dir"); dir.create(dr, showWarnings = FALSE)
  on.exit(unlink(dr, recursive = TRUE)); writeLines("find me", file.path(dr, nm))
  f <- locate_dropped(nm, 8, FALSE, hints = character(), max_seconds = 90)
  expect_equal(unname(f), file.path(dr, nm))
})

test_that("threads typed by the user become a sensible whole number", {
  mx <- parallel::detectCores()
  expect_equal(whole_threads(3), 3L); expect_equal(whole_threads(2.6), 3L); expect_equal(whole_threads("4"), 4L)
  expect_equal(whole_threads(NA), 1L); expect_equal(whole_threads(0), 1L); expect_equal(whole_threads(-5), 1L); expect_equal(whole_threads("abc"), 1L)
  expect_equal(whole_threads(10000), as.integer(mx)); expect_equal(whole_threads(NULL), 1L)
})
