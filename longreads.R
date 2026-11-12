# longreads.R -- counting LONG reads (Oxford Nanopore, PacBio; e.g. Plasmidsaurus) from FASTQ.
#
# The short-read route (tuxedo.R: fastq-mcf -> HISAT2 -> Cufflinks -> Cuffdiff) cannot take these reads: they are 0.5-10 kb long and noisy,
# so quality trimming removes them and HISAT2 is not made for them. Here each file is
#   1. aligned, splice-aware, with minimap2 (-ax splice) and sorted with samtools -> Aligned/<sample>.bam
#   2. counted per gene with featureCounts in long-read mode (-L), from the annotation's exons
# and the result is a genes x samples table of raw counts, exactly what the rest of StRym reads (DESeq2, PCA, trajectories, and the
# entropy steps, whose q-values are then calculated by StRym from the counts: see analysis/derived_diff.R).
#
# Tools: minimap2, samtools, featureCounts (subread). They live in their own WSL environment ("longread", see setup_modules.sh /
# environment-longread.yml) or anywhere on the PATH.  Strandedness: cDNA libraries are not stranded, so reads are counted on either strand.
# (the tool lookup .long_tool() / check_long_tools() lives in tuxedo.R next to the short-read tools, so the design page can check for them)

# samples: data.frame(group, k, file, sample) as from normalize_samples(); returns the counts table (gene + one column per sample)
run_longreads <- function(samples, out_folder, gtf_file, fasta_file, threads = 1, resume = TRUE, log = message) {
  tl <- check_long_tools()
  if (!all(tl$found)) stop("Missing tools for long reads: ", paste(tl$tool[!tl$found], collapse = ", "),
                           ". Install them with setup_modules.sh (see the README, “Long reads”).", call. = FALSE)
  minimap2 <- tl$path[tl$tool == "minimap2"]; samtools <- tl$path[tl$tool == "samtools"]; fcounts <- tl$path[tl$tool == "featureCounts"]
  aligned <- file.path(out_folder, "Aligned"); dir.create(aligned, recursive = TRUE, showWarnings = FALSE)
  fasta_run <- fasta_file; gtf_run <- gtf_file                  # these tools take quoted paths with spaces; no symlink detour (it would put output files in a temp folder)
  threads <- max(1L, as.integer(threads)); sort_threads <- max(1L, min(threads, 4L))
  bams <- character(nrow(samples)); mapped <- numeric(nrow(samples))
  for (i in seq_len(nrow(samples))) {
    s <- samples$sample[i]; f <- samples$file[i]; bam <- file.path(aligned, sprintf("%03d_%s.bam", i, gsub("[^A-Za-z0-9._-]+", "_", s))); bams[i] <- bam      # numbered: "S 1" and "S_1", or "A" and "a", never share a file
    message(sprintf("[SAMPLE %d/%d] %s (%s)", i, nrow(samples), s, samples$group[i]))
    fresh <- resume && file.exists(bam) && file.size(bam) > 0 && !is.na(file.mtime(bam)) && file.mtime(bam) >= file.mtime(f) && file.mtime(bam) >= file.mtime(fasta_file)
    if (fresh) message("Alignment already done: reusing it")
    else {
      tmp <- file.path(aligned, sprintf(".sorttmp_%03d", i))
      cmd <- sprintf("%s -ax splice -t %d --secondary=no %s %s | %s sort -@ %d -m 400M -T %s -O bam -o %s -",
                     shQuote(minimap2), threads, shQuote(fasta_run), shQuote(f), shQuote(samtools), sort_threads, shQuote(tmp), shQuote(bam))
      .run_pipe(cmd, sprintf("minimap2 + samtools sort %s", basename(f)))
      unlink(Sys.glob(paste0(tmp, "*")))
      .run(samtools, c("index", shQuote(bam)), paste("Indexing", basename(bam)))
    }
    fs <- suppressWarnings(system2(samtools, c("flagstat", shQuote(bam)), stdout = TRUE, stderr = FALSE))
    pct <- suppressWarnings(as.numeric(sub(".*\\(([0-9.]+)%.*", "\\1", grep("primary mapped", fs, value = TRUE)[1])))
    if (is.na(pct)) pct <- suppressWarnings(as.numeric(sub(".*\\(([0-9.]+)%.*", "\\1", grep(" mapped \\(", fs, value = TRUE)[1])))
    mapped[i] <- pct
    if (!is.na(pct)) { message(sprintf("%s: %.1f%% of the reads aligned", s, pct))
      if (pct < 30) message(sprintf("WARNING: only %.1f%% of the reads of %s aligned. Wrong genome or annotation for this species, or a contaminated sample?", pct, s)) }
  }
  message("[STEP counting] featureCounts")
  cnt <- file.path(out_folder, "featureCounts.txt")
  .run(fcounts, c("-L", "-T", threads, "--primary", "-a", shQuote(gtf_run), "-t", "exon", "-g", "gene_id", "-o", shQuote(cnt),
                  vapply(bams, function(b) shQuote(b), "")), "featureCounts")
  x <- utils::read.delim(cnt, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  m <- x[, -(1:6), drop = FALSE]; colnames(m) <- samples$sample
  out <- data.frame(gene = x$Geneid, m, check.names = FALSE, stringsAsFactors = FALSE)
  attr(out, "mapped_pct") <- stats::setNames(mapped, samples$sample)
  out
}

# The quantify step for long reads (called from analysis/quantify.R). Same outputs as the short-read route: counts.csv, sample_sheet.csv.
strym_quantify_long <- function(ctx, samples) {
  d <- ctx$design; paired <- FALSE
  if (identical(d$read_type, "paired")) stop("Long reads are single reads: choose “Single-end” as the read type.", call. = FALSE)
  ns <- normalize_samples(samples[, c("group", "file", "mate")], FALSE)
  v <- validate_samples(ns, FALSE, platform = "long"); if (length(v$errors)) stop(paste(c("Input problems:", v$errors), collapse = "\n  "), call. = FALSE)
  for (w in v$warnings) ctx$log(paste("WARNING:", w))
  bad_gz <- check_gz_integrity(ns$file)                                          # a cut-off .gz would otherwise align only the part before the cut and still look fine
  if (length(bad_gz)) stop("These compressed files are damaged or incomplete (a copy that was interrupted?): ", paste(bad_gz, collapse = ", "), ". Copy them again.", call. = FALSE)
  resume <- !isTRUE(strym_normalize_params(d$params)$tux_fresh)
  ctx$log(sprintf("Quantifying %d long-read samples in %d groups: %s", nrow(ns), length(unique(ns$group)), paste(unique(ns$group), collapse = ", ")))
  nrun <- 0L
  follow <- function(m) {
    txt <- conditionMessage(m)
    if (grepl("^\\[SAMPLE [0-9]+/[0-9]+\\]", txt)) { ij <- as.integer(strsplit(sub("^\\[SAMPLE ([0-9]+/[0-9]+)\\].*", "\\1", txt), "/")[[1]])
      ctx$sub_progress(0.9 * (ij[1] - 1) / ij[2], sprintf("aligning sample %d of %d", ij[1], ij[2])) }
    else if (grepl("^\\[STEP counting\\]", txt)) ctx$sub_progress(0.92, "counting reads per gene")
  }
  out <- file.path(ctx$out_dir, "long")
  counts <- withCallingHandlers(run_longreads(ns, out, d$gtf_file, d$fasta_file, d$threads, resume, ctx$log), message = follow)
  tot <- colSums(counts[-1]); med <- stats::median(tot)
  if (any(tot == 0)) stop(sprintf("No reads were counted for: %s. None of their reads aligned to a gene of this annotation (wrong species, genome or annotation file?).", paste(names(tot)[tot == 0], collapse = ", ")), call. = FALSE)
  if (any(tot < 0.05 * med)) ctx$log(sprintf("WARNING: very few reads were counted for %s (under 5%% of the median sample): their results are unreliable.", paste(names(tot)[tot < 0.05 * med], collapse = ", ")))
  utils::write.csv(counts, file.path(ctx$out_dir, "counts.csv"), row.names = FALSE)
  m <- match(normalizePath(ns$file, mustWork = FALSE), normalizePath(samples$file, mustWork = FALSE))
  ss <- data.frame(sample = ns$sample, condition = samples$condition[m], time_point = samples$time_point[m], replicate = samples$replicate[m],
                   file = basename(ns$file), stringsAsFactors = FALSE)
  if (all(is.na(ss$replicate))) ss$replicate <- as.integer(ave(seq_len(nrow(ss)), ss$condition, ss$time_point, FUN = seq_along))
  if (all(!nzchar(ss$time_point))) ss$time_point <- if (length(d$timepoint_names) == 1) d$timepoint_names else "t1"
  utils::write.csv(ss, file.path(ctx$out_dir, "sample_sheet.csv"), row.names = FALSE)
  mp <- attr(counts, "mapped_pct"); if (length(mp)) utils::write.csv(data.frame(sample = names(mp), percent_aligned = mp), file.path(ctx$out_dir, "alignment_rates.csv"), row.names = FALSE)
  list(status = "ok", files = c("counts.csv", "sample_sheet.csv", "alignment_rates.csv", "long/featureCounts.txt"),
       message = sprintf("%s genes × %d long-read samples counted (minimap2 + featureCounts); q-values for the entropy steps are calculated from these counts.", format(nrow(counts), big.mark = ","), ncol(counts) - 1),
       data_path = file.path(ctx$out_dir, "counts.csv"), sample_sheet = file.path(ctx$out_dir, "sample_sheet.csv"))
}
