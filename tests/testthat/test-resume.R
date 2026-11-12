# Resuming an interrupted Tux run: finished, intact parts are reused; anything doubtful is redone. Real tools, synthetic reads.
skip_if(!all(check_tools()$found), "Tuxedo tools not on PATH")

tux_fixture <- function() {
  root <- tempfile("rs"); dir.create(root)
  system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(here, "make_synthetic.R")), "--out", shQuote(root), "--layout flat --gz --seed 7"), stdout = FALSE)
  list(root = root, inputs = file.path(root, "inputs"), gtf = file.path(root, "genes.gtf"), fasta = file.path(root, "genome.fa"), out = file.path(root, "tux out"))   # a space in the path on purpose
}
go <- function(fx, resume = TRUE, ...) {
  a <- list(inputs_folder = fx$inputs, gtf_file = fx$gtf, fasta_file = fx$fasta, output_folder = fx$out, read_type = "single", threads = 2,
            adapter_fasta = file.path(here, "IlluminaAdaptorSeq.fasta"), resume = resume); b <- list(...); for (k in names(b)) a[[k]] <- b[[k]]
  testthat::capture_messages(do.call(run_tuxedo, a))
}
n_msg <- function(m, pat) sum(grepl(pat, m))
count_ran <- function(m) c(align = n_msg(m, "^hisat2 \\+ samtools sort"), cufflinks = n_msg(m, "^Cufflinks [^ ]+\\.\\.\\."), merge = n_msg(m, "^Cuffmerge\\.\\.\\."), diff = n_msg(m, "^Cuffdiff\\.\\.\\."))
stopifnot(sum(grepl("^Cufflinks [^ ]+\\.\\.\\.", c("Cufflinks CTL1_Fe...\n", "Cufflinks TRT2_Fe...\n", "Cufflinks successful on all files\n"))) == 2)   # sample names are not all alike; captured messages end in a newline

test_that("a finished run is recognized as complete: nothing is redone, results are untouched", {
  fx <- tux_fixture(); m1 <- go(fx, resume = FALSE)
  expect_equal(count_ran(m1), c(align = 6, cufflinks = 6, merge = 1, diff = 1))
  ge <- file.path(fx$out, "diff_out", "gene_exp.diff"); before <- tools::md5sum(ge); t0 <- file.mtime(ge)
  m2 <- go(fx)
  expect_equal(count_ran(m2), c(align = 0, cufflinks = 0, merge = 0, diff = 0))
  expect_true(any(grepl("Reusing finished work from an earlier run \\(alignments 6/6, transcript assembly 6/6, merge, differential expression\\). Still to do: nothing", m2)))
  expect_equal(n_msg(m2, "Alignment already done"), 6); expect_equal(n_msg(m2, "Assembly already done"), 6)
  expect_equal(unname(tools::md5sum(ge)), unname(before)); expect_equal(file.mtime(ge), t0)
  expect_true(any(grepl("\\[STEP 4/4\\]", m2))); expect_true(any(grepl("^Done\\. Results in", m2)))                  # the step markers StRym's progress bar follows are all there
  expect_true(any(grepl("^\\[SAMPLE 6/6\\]", m2)))
})

test_that("an interruption in Cuffdiff (result tables empty, as after the reboot) redoes only Cuffdiff", {
  fx <- tux_fixture(); go(fx, resume = FALSE)
  ge <- file.path(fx$out, "diff_out", "gene_exp.diff"); ref <- read.delim(ge); ref_rg <- read.delim(file.path(fx$out, "diff_out", "genes.read_group_tracking"))
  for (f in list.files(file.path(fx$out, "diff_out"), full.names = TRUE)) if (grepl("diff$|tracking$|\\.info$", f) && !grepl("bias_params|var_model", f)) writeLines(character(), f)   # all tables 0 bytes
  m <- go(fx)
  expect_equal(count_ran(m), c(align = 0, cufflinks = 0, merge = 0, diff = 1)); expect_true(any(grepl("Still to do: Cuffdiff\\.", m)))
  expect_equal(nrow(read.delim(ge)), nrow(ref)); expect_equal(read.delim(ge)$q_value, ref$q_value)                       # the same, complete answer
  expect_equal(nrow(read.delim(file.path(fx$out, "diff_out", "genes.read_group_tracking"))), nrow(ref_rg))
  # and the half-written folder of the failed attempt is not mixed in: a stray file from it is gone
  writeLines("x", file.path(fx$out, "diff_out", "stray_from_failed_attempt.txt")); for (f in list.files(file.path(fx$out, "diff_out"), full.names = TRUE, pattern = "diff$")) writeLines(character(), f)
  go(fx); expect_false(file.exists(file.path(fx$out, "diff_out", "stray_from_failed_attempt.txt")))
})

test_that("a damaged BAM is detected: that sample is realigned and everything after it is redone", {
  fx <- tux_fixture(); go(fx, resume = FALSE)
  bam <- file.path(fx$out, "Aligned", "hisat2OutSortedCTL_2.bam"); sz <- file.size(bam)
  con <- file(bam, "r+b"); seek(con, sz %/% 2); writeBin(raw(0), con); close(con); system2("truncate", c("-s", sz %/% 2, shQuote(bam)))     # cut in half: no end-of-file marker
  expect_false(system2("samtools", c("quickcheck", shQuote(bam)), stdout = FALSE, stderr = FALSE) == 0)
  st <- tux_resume_status(infer_samples(fx$inputs), fx$out, fx$gtf, fx$fasta); expect_equal(sum(st$align), 5); expect_false(st$merge); expect_false(st$cuffdiff)
  expect_equal(which(!st$align), 2)
  m <- go(fx); expect_equal(count_ran(m), c(align = 1, cufflinks = 1, merge = 1, diff = 1))
  expect_equal(system2("samtools", c("quickcheck", shQuote(bam)), stdout = FALSE, stderr = FALSE), 0)
  expect_equal(count_ran(go(fx)), c(align = 0, cufflinks = 0, merge = 0, diff = 0))                                       # and now it is whole again
})

test_that("a changed reference, changed cleaning settings or different samples are never reused", {
  fx <- tux_fixture(); go(fx, resume = FALSE); smp <- infer_samples(fx$inputs)
  expect_true(all(tux_resume_status(smp, fx$out, fx$gtf, fx$fasta)$cuffdiff))
  st <- tux_resume_status(smp, fx$out, fx$gtf, fx$fasta, quality_thresh = 20); expect_false(any(st$align)); expect_match(st$why, "cleaning settings")
  st <- tux_resume_status(smp, fx$out, fx$gtf, fx$fasta, window_size = 7); expect_false(any(st$align))
  expect_match(tux_resume_text(st), "Nothing to reuse \\(the cleaning settings or the reference differ")
  st <- tux_resume_status(smp[-1, ], fx$out, fx$gtf, fx$fasta); expect_false(any(st$align)); expect_match(st$why, "samples differ")
  st <- tux_resume_status(smp, fx$out, fx$gtf, fx$fasta, paired = TRUE); expect_false(any(st$align))                       # read type changed (paired flag is part of the record)
  # a newer annotation: alignments stay valid (they do not depend on it), assembly onwards is redone
  Sys.setFileTime(fx$gtf, Sys.time() + 5)
  st <- tux_resume_status(smp, fx$out, fx$gtf, fx$fasta); expect_true(all(st$align)); expect_false(any(st$cufflinks)); expect_false(st$merge); expect_false(st$cuffdiff)
  # a different annotation of another size: the record says it is not the same reference
  writeLines(c(readLines(fx$gtf), readLines(fx$gtf)[1]), fx$gtf)
  expect_false(any(tux_resume_status(smp, fx$out, fx$gtf, fx$fasta)$align))
  # a reads file replaced by a newer one: that sample is not reused
  fx2 <- tux_fixture(); go(fx2, resume = FALSE); smp2 <- infer_samples(fx2$inputs)
  Sys.setFileTime(smp2$file[3], Sys.time() + 5); st <- tux_resume_status(smp2, fx2$out, fx2$gtf, fx2$fasta)
  expect_equal(which(!st$align), 3); expect_false(st$merge)
  expect_false(any(tux_resume_status(smp2, file.path(fx2$root, "nowhere"), fx2$gtf, fx2$fasta)$align))                     # no earlier run at all
})

test_that("a run made before the settings record existed can still be resumed with the default settings only", {
  fx <- tux_fixture(); go(fx, resume = FALSE); smp <- infer_samples(fx$inputs)
  unlink(file.path(fx$out, "run_params.json"))                                                                           # as left by the earlier version
  expect_true(tux_resume_status(smp, fx$out, fx$gtf, fx$fasta)$cuffdiff)
  expect_false(any(tux_resume_status(smp, fx$out, fx$gtf, fx$fasta, quality_thresh = 25)$align))
  m <- go(fx); expect_equal(count_ran(m), c(align = 0, cufflinks = 0, merge = 0, diff = 0)); expect_true(file.exists(file.path(fx$out, "run_params.json")))   # the record is written now
})

test_that("without resume an existing output is still refused; StRym's step resumes by default and 'from scratch' replaces everything", {
  fx <- tux_fixture(); go(fx, resume = FALSE)
  expect_error(go(fx, resume = FALSE), "Previous outputs exist")
  d <- strym_default_design(); d$project_name <- "rs"; d$project_dir <- file.path(fx$root, "proj"); d$data_type <- "bulk"; d$data_kind <- "fastq"
  d$data_path <- fx$inputs; d$gtf_file <- fx$gtf; d$fasta_file <- fx$fasta; d$read_type <- "single"; d$threads <- 2L; d$n_conditions <- 2L
  msgs <- testthat::capture_messages(strym_run(d, steps = "quantify", app_dir = here)); expect_equal(count_ran(msgs)[["align"]], 6)
  tux <- file.path(d$project_dir, "01_quantification", "tux")
  # an interruption inside Cuffdiff, then the same Run again
  for (f in list.files(file.path(tux, "diff_out"), full.names = TRUE, pattern = "diff$|tracking$")) writeLines(character(), f)
  unlink(file.path(d$project_dir, "01_quantification", "counts.csv"))
  m2 <- testthat::capture_messages(res2 <- strym_run(d, steps = "quantify", app_dir = here))
  expect_equal(res2$summary$status, "ok"); expect_equal(count_ran(m2), c(align = 0, cufflinks = 0, merge = 0, diff = 1))
  expect_true(file.exists(file.path(d$project_dir, "01_quantification", "counts.csv"))); expect_gt(nrow(read.csv(file.path(d$project_dir, "01_quantification", "counts.csv"))), 5)
  # "start the Tux step from scratch"
  d$params <- list(tux_fresh = TRUE)
  m3 <- testthat::capture_messages(strym_run(d, steps = "quantify", app_dir = here)); expect_equal(count_ran(m3), c(align = 6, cufflinks = 6, merge = 1, diff = 1))
  expect_false(any(grepl("Reusing finished work", m3)))
})

test_that("the resume wording is plain", {
  st <- list(align = c(TRUE, TRUE, TRUE, TRUE), cufflinks = rep(TRUE, 4), merge = TRUE, cuffdiff = FALSE, why = "")
  expect_equal(tux_resume_text(st), "Reusing finished work from an earlier run (alignments 4/4, transcript assembly 4/4, merge). Still to do: Cuffdiff.")
  st2 <- list(align = c(TRUE, FALSE), cufflinks = c(FALSE, FALSE), merge = FALSE, cuffdiff = FALSE, why = "")
  expect_equal(tux_resume_text(st2), "Reusing finished work from an earlier run (alignments 1/2, transcript assembly 0/2). Still to do: alignment, assembly, merge, Cuffdiff.")
  expect_match(tux_resume_text(list(align = FALSE, cufflinks = FALSE, merge = FALSE, cuffdiff = FALSE, why = "")), "^Nothing to reuse: starting from the beginning")
})
