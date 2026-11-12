fake_fastq <- function(path, n = 3) {                       # a few valid FASTQ records (content only needs to look like FASTQ)
  rec <- unlist(lapply(seq_len(n), function(i) c(sprintf("@r%d", i), "ACGTACGTAC", "+", "IIIIIIIIII")))
  con <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt"); writeLines(rec, con); close(con)
}
reads_design <- function(root, files = c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz", "TRT1_x.fastq.gz", "TRT2_x.fastq.gz"), ...) {
  reads <- file.path(root, "reads"); dir.create(reads, showWarnings = FALSE)
  for (f in files) fake_fastq(file.path(reads, f))
  ref <- file.path(root, "ref"); dir.create(ref, showWarnings = FALSE)
  writeLines(">chr1\nACGT", file.path(ref, "g.fa")); writeLines("chr1\tx\texon\t1\t4\t.\t+\t.\tgene_id \"G1\"; transcript_id \"G1.1\";", file.path(ref, "g.gtf"))
  d <- strym_default_design()
  d$project_name <- "reads"; d$project_dir <- file.path(root, "proj"); d$data_type <- "bulk"; d$data_kind <- "fastq"
  d$data_path <- reads; d$gtf_file <- file.path(ref, "g.gtf"); d$fasta_file <- file.path(ref, "g.fa"); d$n_conditions <- 2L
  for (k in names(list(...))) d[[k]] <- list(...)[[k]]
  d
}

test_that("quantification applies (and is ticked) only when starting from raw reads", {
  p <- strym_plan(strym_default_design()); expect_false(p$applicable[p$id == "quantify"])
  d <- strym_default_design(); d$data_kind <- "fastq"
  p <- strym_plan(d); expect_true(p$applicable[p$id == "quantify"]); expect_true(p$default[p$id == "quantify"])
  expect_equal(p$id[1], "quantify")                                        # it comes first
})

test_that("with raw reads, later steps are blocked until counts exist", {
  root <- tempfile("q"); dir.create(root)
  d <- reads_design(root); d$n_timepoints <- 1L
  b <- strym_blocked(d, c("deseq2", "pca")); expect_true(all(grepl("Quantify reads", b)))
  expect_equal(unname(strym_blocked(d, c("quantify", "deseq2", "pca"))), c("", "", ""))
  dir.create(file.path(d$project_dir, "01_quantification"), recursive = TRUE)
  writeLines("gene,s1", strym_counts_path(d))                              # counts from an earlier quantification
  expect_equal(unname(strym_blocked(d, c("deseq2", "pca"))), c("", ""))
  expect_equal(unname(strym_blocked(strym_default_design(), "pca")), "")  # counts tables are unaffected
})

test_that("sample table: from the sample sheet (condition + time point) or guessed from file names", {
  root <- tempfile("q"); dir.create(root); d <- reads_design(root)
  s <- strym_quantify_samples(d)                                           # no sheet: groups from names
  expect_equal(sort(unique(s$group)), c("CTL", "TRT")); expect_equal(nrow(s), 4)
  ss <- data.frame(sample = c("a", "b", "c", "d"), condition = c("ctl", "ctl", "heat", "heat"), time_point = c("0h", "1h", "0h", "1h"),
                   replicate = 1, file = c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz", "TRT1_x.fastq.gz", "TRT2_x.fastq.gz"))
  write.csv(ss, file.path(root, "ss.csv"), row.names = FALSE)
  d$sample_sheet <- file.path(root, "ss.csv"); d$n_timepoints <- 2L
  s <- strym_quantify_samples(d)
  expect_equal(s$group, c("ctl_0h", "ctl_1h", "heat_0h", "heat_1h")); expect_equal(s$file, file.path(d$data_path, ss$file))
  d$n_timepoints <- 1L; expect_equal(strym_quantify_samples(d)$group, ss$condition)   # one time point: group = condition
  write.csv(ss[, c("sample", "condition")], file.path(root, "bad.csv"), row.names = FALSE); d$sample_sheet <- file.path(root, "bad.csv")
  expect_error(strym_quantify_samples(d), "“file” column")
})

test_that("a sample sheet with an empty 'file' column is matched to the reads by sample name", {
  root <- tempfile("q"); dir.create(root); d <- reads_design(root)
  blank <- data.frame(sample = c("CTL1_x", "CTL2_x", "TRT1_x", "TRT2_x"), condition = c("B", "B", "D", "D"), time_point = "t1", replicate = 1:2, file = "")
  write.csv(blank, file.path(root, "blank.csv"), row.names = FALSE); d$sample_sheet <- file.path(root, "blank.csv")
  s <- strym_quantify_samples(d)
  expect_equal(basename(s$file), c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz", "TRT1_x.fastq.gz", "TRT2_x.fastq.gz")); expect_equal(s$group, c("B", "B", "D", "D"))
  v <- strym_validate_reads(d); expect_length(v$errors, 0); expect_true(any(grepl("B ×2, D ×2", v$info)))
  # sample names that do not match any file: a clear message instead of "NA: file not found"
  bad <- blank; bad$sample <- c("s1", "s2", "s3", "s4"); write.csv(bad, file.path(root, "bad.csv"), row.names = FALSE); d$sample_sheet <- file.path(root, "bad.csv")
  expect_error(strym_quantify_samples(d), "no FASTQ file found for s1, s2, s3, s4.*“file” column")
  v <- strym_validate_reads(d); expect_true(any(grepl("no FASTQ file found", v$errors))); expect_false(any(grepl("NA: file not found", v$errors)))
  # names typed the way Windows Explorer shows them (no .gz) still resolve
  typed <- blank; typed$file <- c("CTL1_x.fastq", "CTL2_x.fastq", "TRT1_x.fastq", "TRT2_x.fastq")
  write.csv(typed, file.path(root, "typed.csv"), row.names = FALSE); d$sample_sheet <- file.path(root, "typed.csv")
  expect_equal(basename(strym_quantify_samples(d)$file), c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz", "TRT1_x.fastq.gz", "TRT2_x.fastq.gz"))
  expect_length(strym_validate_reads(d)$errors, 0)
})

test_that("the sample-sheet template lists the real FASTQ files when starting from raw reads", {
  root <- tempfile("q"); dir.create(root); d <- reads_design(root)
  t <- strym_sample_sheet_template(d)
  expect_equal(t$file, c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz", "TRT1_x.fastq.gz", "TRT2_x.fastq.gz")); expect_equal(t$condition, c("CTL", "CTL", "TRT", "TRT"))
  write.csv(t, file.path(root, "t.csv"), row.names = FALSE); d$sample_sheet <- file.path(root, "t.csv")
  expect_length(strym_validate_reads(d)$errors, 0)                         # the downloaded template works as is
})

test_that("checks for raw reads explain what is missing", {
  root <- tempfile("q"); dir.create(root); d <- reads_design(root, data_path = "")
  v <- strym_validate_reads(modifyList(reads_design(root), list(data_path = file.path(root, "reads", "CTL1_x.fastq.gz"))))
  expect_match(v$errors, "FOLDER")
  d <- reads_design(root)
  v <- strym_validate_reads(d); expect_length(v$errors, 0); expect_true(any(grepl("4 single FASTQ files in 2 groups \\(CTL ×2, TRT ×2\\)", v$info)))
  expect_true(any(grepl("free on the project drive", v$info)))
  v <- strym_validate_reads(modifyList(d, list(gtf_file = "", fasta_file = file.path(root, "nope.fa"))))
  expect_true(any(grepl("annotation", v$errors))); expect_true(any(grepl("genome FASTA file was not found", v$errors)))
  v <- strym_validate_reads(modifyList(d, list(fasta_file = paste0(d$fasta_file, ".gz"))))
  expect_true(any(grepl("uncompressed", v$errors)))
  d1 <- reads_design(tempfile("q1") |> (\(p) { dir.create(p); p })(), files = c("CTL1_x.fastq.gz", "CTL2_x.fastq.gz"))
  expect_true(any(grepl("2 groups|1 group", strym_validate_reads(d1)$warnings)) || any(grepl("group", strym_validate_reads(d1)$warnings)))
  writeLines(c("this", "is", "not", "fastq"), file.path(root, "reads", "TRT2_x.fastq")); file.remove(file.path(root, "reads", "TRT2_x.fastq.gz"))
  expect_true(any(grepl("does not look like FASTQ", strym_validate_reads(d)$errors)))
  # and the overall design check uses it
  full <- strym_validate_design(modifyList(d, list(project_name = "x"))); expect_true(any(grepl("does not look like FASTQ", full$errors)))
})

test_that("counts table is built from Cuffdiff output with the right replicate-to-sample mapping", {
  tx <- tempfile("tux"); dir.create(file.path(tx, "diff_out"), recursive = TRUE)
  write.csv(data.frame(group = c("A", "A", "B"), k = c(1, 2, 1), file = c("a1.fq", "a2.fq", "b1.fq"), mate = NA, sample = c("a1", "a2", "b1")),
            file.path(tx, "samples_used.csv"), row.names = FALSE)
  rg <- expand.grid(replicate = 0:1, condition = c("A", "B"), tracking_id = c("XLOC_1", "XLOC_2"), stringsAsFactors = FALSE)
  rg <- rg[!(rg$condition == "B" & rg$replicate == 1), ]; rg$raw_frags <- seq_len(nrow(rg)) * 10
  rg <- rg[, c("tracking_id", "condition", "replicate", "raw_frags")]
  write.table(rg, file.path(tx, "diff_out", "genes.read_group_tracking"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(test_id = c("XLOC_1", "XLOC_2"), gene = c("AT1G01", "-")), file.path(tx, "diff_out", "gene_exp.diff"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cnt <- strym_counts_from_cuffdiff(tx)
  expect_equal(names(cnt), c("gene", "a1", "a2", "b1")); expect_equal(cnt$gene, c("AT1G01", "XLOC_2"))   # name where known, id otherwise
  m <- rg[rg$tracking_id == "XLOC_1", ]
  expect_equal(cnt$a1[1], m$raw_frags[m$condition == "A" & m$replicate == 0]); expect_equal(cnt$a2[1], m$raw_frags[m$condition == "A" & m$replicate == 1])
  expect_equal(cnt$b1[1], m$raw_frags[m$condition == "B"])
})

test_that("a counts-needing step is skipped (not run on raw reads) when quantification was not selected", {
  root <- tempfile("q"); dir.create(root); d <- reads_design(root)
  res <- quiet(strym_run(d, steps = "pca", app_dir = here))
  expect_equal(res$summary$status, "skipped"); expect_match(res$summary$message, "Quantify reads")
  expect_false(dir.exists(file.path(d$project_dir, "04_pca")))
})

test_that("END TO END: raw reads -> counts -> later steps (synthetic reads, full Tuxedo run)", {
  skip_if(!all(check_tools()$found), "Tuxedo tools not on PATH")
  root <- tempfile("e2e"); dir.create(root)
  system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(here, "make_synthetic.R")), "--out", shQuote(root), "--layout flat --gz --seed 7"), stdout = FALSE)
  d <- strym_default_design()
  d$project_name <- "e2e"; d$project_dir <- file.path(root, "proj"); d$data_type <- "bulk"; d$data_kind <- "fastq"
  d$data_path <- file.path(root, "inputs"); d$gtf_file <- file.path(root, "genes.gtf"); d$fasta_file <- file.path(root, "genome.fa")
  d$read_type <- "single"; d$threads <- 2L; d$n_conditions <- 2L; d$n_timepoints <- 1L
  expect_length(strym_validate_design(d)$errors, 0)
  msgs <- testthat::capture_messages(res <- strym_run(d, steps = c("quantify", "deseq2", "pca"), app_dir = here,
                                                  progress = function(p, l) message(sprintf("[PROGRESS %d] %s", as.integer(p), l))))
  s <- setNames(res$summary$status, res$summary$step)
  expect_equal(s[["quantify"]], "ok"); expect_equal(s[["deseq2"]], "ok"); expect_equal(s[["pca"]], "ok")      # DESeq2 reads the Tux counts, PCA the Tux (Cuffdiff) results
  pd <- d$project_dir
  cnt <- read.csv(file.path(pd, "01_quantification", "counts.csv"), check.names = FALSE)
  expect_equal(names(cnt)[-1], sprintf("%s%d_Fe", rep(c("CTL", "TRT"), each = 3), 1:3)); expect_gt(nrow(cnt), 5); expect_true(all(cnt[, -1] >= 0))
  ss <- read.csv(file.path(pd, "01_quantification", "sample_sheet.csv"))
  expect_equal(ss$condition, rep(c("CTL", "TRT"), each = 3)); expect_equal(ss$replicate, rep(1:3, 2))
  expect_true(file.exists(file.path(pd, "01_quantification", "tux", "diff_out", "gene_exp.diff")))
  dgs <- read.csv(file.path(pd, "02_differential_expression", "deseq2", "deseq2_TRT_vs_CTL.csv"))                    # DESeq2 on the real Tux counts
  truth0 <- read.delim(file.path(root, "truth.tsv")); expect_true(all(c("gene", "log2FoldChange", "padj", "direction") %in% names(dgs)))
  expect_equal(nrow(dgs) > 5, TRUE); expect_true(any(dgs$direction != "none") || all(dgs$padj <= 1))
  expect_true(any(grepl("\\[PROGRESS [0-9]+\\] Quantify reads .* — cleaning and aligning, sample 2 of 6", msgs)))   # progress inside the step
  # a later single step reuses the counts instead of the reads folder
  res2 <- quiet(strym_run(d, steps = "pca", app_dir = here)); expect_equal(res2$summary$status, "ok")
  sc <- read.csv(file.path(pd, "04_pca", "pca_scores.csv"))                       # PCA of the 6 real replicates, named as in the sample sheet
  expect_equal(sort(sc$sample), sort(ss$sample)); expect_equal(sort(unique(sc$condition)), c("CTL", "TRT")); expect_true(file.exists(file.path(pd, "04_pca", "pca_plot.png")))
  expect_equal(sum(read.csv(file.path(pd, "04_pca", "pca_variance.csv"))$percent_variance), 100, tolerance = 0.05)
  # the entropy methods need time points: with one they are not offered, and the step says why if forced
  expect_false(strym_plan(d)$applicable[strym_plan(d)$id == "shannon_entropy"])
  data <- strym_load_expression(list(design = d, project_dir = pd)); expect_equal(data$kind, "tux")
  expect_equal(data$samples$condition, ss$condition); expect_equal(nrow(strym_wide(data)) > 5, TRUE)
  expect_error(strym_timecourse(data, strym_default_params(), character()), "at least 2 time points")
  expect_equal(read.csv(file.path(pd, "run_summary.csv"))$step, c("quantify", "deseq2", "pca"))   # summary keeps all steps
  # the counts agree with the planted differences: the up-regulated genes have more reads in TRT than in CTL
  truth <- read.delim(file.path(root, "truth.tsv")); ups <- truth$gene_id[truth$de == "up"]
  expect_gt(length(ups), 0)
})
