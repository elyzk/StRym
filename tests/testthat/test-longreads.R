# Long reads (Oxford Nanopore / PacBio, e.g. Plasmidsaurus): minimap2 + featureCounts, then the same downstream steps.
test_that("long reads are counted per gene and the planted genes are found", {
  skip_if(!all(check_long_tools()$found), "minimap2 / featureCounts not installed")
  e <- file.path(here, "example_data", "protocol", "exp6_longread_tiny"); skip_if(!dir.exists(e), "example data missing")
  proj <- tempfile("long_"); dir.create(proj)
  d <- strym_normalize_design(list(project_name = "t", project_dir = proj, data_type = "bulk", data_kind = "fastq", data_path = file.path(e, "inputs"),
                                   gtf_file = file.path(e, "genes.gtf"), fasta_file = file.path(e, "genome.fa"), read_type = "single", read_platform = "long",
                                   n_conditions = 2, n_timepoints = 1, threads = 2))
  expect_equal(strym_validate_design(d)$errors, character())
  out <- file.path(proj, "01_quantification"); dir.create(out)
  ctx <- list(design = d, data_path = d$data_path, sample_sheet = "", out_dir = out, project_dir = proj, app_dir = here, prereq = list(), log = function(...) invisible(),
              sub_progress = function(...) invisible())
  source(file.path(here, "analysis", "quantify.R")); source(file.path(here, "longreads.R"))
  r <- strym_step_quantify(ctx); expect_equal(r$status, "ok")
  cn <- utils::read.csv(file.path(out, "counts.csv")); expect_equal(nrow(cn), 12); expect_equal(ncol(cn), 7)
  truth <- utils::read.delim(file.path(e, "truth.tsv")); up <- truth$gene_id[truth$de == "up"]
  m <- as.matrix(cn[, -1]); rownames(m) <- cn$gene
  fc <- rowMeans(m[, 4:6]) / rowMeans(m[, 1:3]); expect_true(all(fc[up] > 3))                           # the planted up-regulated genes
  expect_true(file.exists(file.path(out, "sample_sheet.csv")))
})

test_that("short reads set to long-read mode and long reads set to short-read mode are both caught", {
  d <- tempfile("lr_"); dir.create(d)
  fq <- function(f, len) { l <- unlist(lapply(1:6, function(i) c(sprintf("@r%d", i), strrep("A", len), "+", strrep("I", len)))); writeLines(l, f) }
  fq(file.path(d, "a.fastq"), 3000); fq(file.path(d, "b.fastq"), 3000)
  s <- normalize_samples(data.frame(group = c("A", "B"), file = file.path(d, c("a.fastq", "b.fastq"))), FALSE)
  expect_match(paste(validate_samples(s, FALSE, platform = "short")$errors, collapse = " "), "Choose .Long reads")
  expect_equal(validate_samples(s, FALSE, platform = "long")$errors, character())
  fq(file.path(d, "c.fastq"), 100); fq(file.path(d, "e.fastq"), 100)
  s2 <- normalize_samples(data.frame(group = c("A", "B"), file = file.path(d, c("c.fastq", "e.fastq"))), FALSE)
  expect_match(paste(validate_samples(s2, FALSE, platform = "long")$warnings, collapse = " "), "choose Short reads")
})