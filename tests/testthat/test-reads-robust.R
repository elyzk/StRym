# Real sequencer output is untidy. The read checks must catch what would otherwise give silent wrong results.
fq_write <- function(path, n = 20, len = 50, qual = "I", id = function(i) sprintf("@read%d", i), gz = grepl("\\.gz$", path)) {
  l <- unlist(lapply(seq_len(n), function(i) c(id(i), paste(sample(c("A", "C", "G", "T"), len, TRUE), collapse = ""), "+", strrep(qual, len))))
  if (gz) { con <- gzfile(path, "wt"); writeLines(l, con); close(con) } else writeLines(l, path)
  invisible(path)
}
mkdir <- function() { d <- tempfile("reads_"); dir.create(d); d }

test_that("long reads (Nanopore / PacBio / Plasmidsaurus) are refused with a pointer to the counts table", {
  d <- mkdir(); fq_write(file.path(d, "a1.fastq"), n = 5, len = 3000); fq_write(file.path(d, "b1.fastq"), n = 5, len = 3000)
  s <- normalize_samples(data.frame(group = c("A", "B"), file = file.path(d, c("a1.fastq", "b1.fastq"))), FALSE)
  v <- validate_samples(s, FALSE); expect_match(paste(v$errors, collapse = " "), "long.*Nanopore.*counts table")
})

test_that("mates must be readable and belong together", {
  d <- mkdir(); set.seed(1)
  fq_write(file.path(d, "s1_R1.fastq"), id = function(i) sprintf("@r%d/1", i)); fq_write(file.path(d, "s1_R2.fastq"), id = function(i) sprintf("@r%d/2", i))
  fq_write(file.path(d, "s2_R1.fastq"), id = function(i) sprintf("@x%d 1:N:0", i)); fq_write(file.path(d, "s2_R2.fastq"), id = function(i) sprintf("@y%d 2:N:0", i))   # not the same reads
  s <- normalize_samples(data.frame(group = c("A", "B"), file = file.path(d, c("s1_R1.fastq", "s2_R1.fastq")), mate = file.path(d, c("s1_R2.fastq", "s2_R2.fastq"))), TRUE)
  v <- validate_samples(s, TRUE); expect_match(paste(v$errors, collapse = " "), "s2_R1.fastq and s2_R2.fastq do not belong together"); expect_false(grepl("s1_R1", paste(v$errors, collapse = " ")))
  writeLines(character(), file.path(d, "s2_R2.fastq")); v2 <- validate_samples(s, TRUE); expect_match(paste(v2$errors, collapse = " "), "s2_R2.fastq: .*empty")
})

test_that("a name that only starts with A_ is not an A_/B_ mate pair; extensions in any case are found", {
  d <- mkdir(); fq_write(file.path(d, "A_thaliana_c1_R1.fastq")); fq_write(file.path(d, "A_thaliana_c1_R2.fastq")); fq_write(file.path(d, "X.FASTQ.GZ")); fq_write(file.path(d, "Y.Fq"))
  mi <- .mate_info(list.files(d, full.names = TRUE, pattern = "A_thaliana")); expect_equal(mi$first, c(TRUE, FALSE)); expect_equal(basename(mi$mate[1]), "A_thaliana_c1_R2.fastq")
  expect_length(list.files(d, pattern = .FQ_EXT, full.names = TRUE), 4)
})

test_that("Undetermined and index files are not samples; lanes and paired-looking files are warned about", {
  d <- mkdir(); for (f in c("Ctrl-1_S1_L001_R1_001.fastq", "Ctrl-1_S1_L001_R2_001.fastq", "Ctrl-1_S1_L001_I1_001.fastq", "Undetermined_S0_L001_R1_001.fastq")) fq_write(file.path(d, f))
  fq_write(file.path(d, "Ctrl-1_S1_L002_R1_001.fastq")); fq_write(file.path(d, "Ctrl-1_S1_L002_R2_001.fastq"))
  l <- quiet(.list_fastq(d, TRUE)); expect_equal(sort(basename(l$file)), c("Ctrl-1_S1_L001_R1_001.fastq", "Ctrl-1_S1_L002_R1_001.fastq"))
  s <- data.frame(group = c("Ctrl", "Ctrl"), k = 1:2, file = l$file, mate = l$mate, sample = c("a", "b"), stringsAsFactors = FALSE)
  v <- validate_samples(s, TRUE); expect_match(paste(v$warnings, collapse = " "), "split over sequencing lanes.*cat ")
  s1 <- normalize_samples(data.frame(group = c("A", "B"), file = list.files(d, pattern = "^Ctrl-1_S1_L001_R[12]_001.fastq$", full.names = TRUE)), FALSE)
  expect_match(paste(validate_samples(s1, FALSE)$warnings, collapse = " "), "looks like paired-end")
})

test_that("a truncated compressed file is detected before the run, a comma in a path is made safe, the Phred offset is explicit", {
  d <- mkdir(); f <- file.path(d, "ok.fastq.gz"); fq_write(f, n = 200); b <- file.path(d, "cut.fastq.gz"); fq_write(b, n = 2000, len = 100)
  raw <- readBin(b, "raw", file.size(b)); writeBin(raw[seq_len(length(raw) %/% 2)], b)
  expect_equal(check_gz_integrity(c(f, b)), "cut.fastq.gz")
  p <- file.path(d, "Smith,Lab"); dir.create(p); expect_false(identical(.safe(p), p)); .cleanup_links()
  expect_equal(.detect_phred(f), 33L); g <- file.path(d, "p33.fastq"); fq_write(g, n = 10, qual = "I"); expect_equal(.detect_phred(g), 33L)     # all bases Q40 in Phred+33 (what fastq-mcf alone would call +64)
  h <- file.path(d, "real64.fastq"); l <- readLines(f <- fq_write(h, n = 4)); l[4] <- paste0(strrep("h", 25), strrep("B", 25)); writeLines(l, h); expect_equal(.detect_phred(h), 64L)
})

test_that("a sample whose condition is blank stops the sample sheet instead of being dropped", {
  d <- mkdir(); fq_write(file.path(d, "a.fastq")); fq_write(file.path(d, "b.fastq"))
  ss <- file.path(d, "s.csv"); writeLines(c("file,condition", "a.fastq,Ctrl", "b.fastq,"), ss)
  des <- list(data_kind = "fastq", data_path = d, sample_sheet = ss, read_type = "single", n_conditions = 2, n_timepoints = 1)
  expect_error(strym_quantify_samples(des), "no condition")
})