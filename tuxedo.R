# tuxedo.R -- R port of the TUX (Tuxedo) tab of TuxNet.
#
# Pipeline (same tools and flags as the MATLAB app):
#   1. Clean + align, one sample at a time:
#        ea-utils fastq-mcf (or fastp, e.g. on macOS) -> hisat2 --dta-cufflinks | samtools sort (BAM)
#   2. cufflinks (per sample)   3. cuffmerge   4. cuffdiff
#
# Designed for very large data on external drives:
#   * inputs are only ever READ; all outputs go to a separate `output_folder`
#   * no raw SAM is written (hisat2 is piped into samtools sort -> BAM)
#   * each sample's cleaned reads are deleted right after alignment
#   * free disk space is estimated and checked before starting
#
# The external tools are Linux/macOS binaries (that's why the Windows standalone of
# TuxNet has no TUX tab). Run this in Linux/WSL with the tools on PATH (or under
# `tool_dir`, laid out like the TuxNet install: ea-utils/clipper/fastq-mcf,
# hisat2/{hisat2,hisat2-build}, samtools/samtools, cufflinks/{cufflinks,cuffmerge,cuffdiff}).
#
# Samples: everything runs from a "sample sheet" data.frame(group, file[, mate]).
# You can pass your own (`samples=`) or let infer_samples() build one from the inputs folder:
#   * flat folder of FASTQ files  -> group inferred from file names (C2B1_x, C2B2_x -> "C2B")
#   * folder with one sub-folder per condition -> group = sub-folder name (TuxNet layout)
# FASTQ files may end .fastq/.fq, optionally .gz. Paired-end mates: A_x/B_x, x_R1/x_R2, x_1/x_2.
#
# Output layout (inside output_folder): Cleaned/ Aligned/ Cufflinks/ merged_asm/ diff_out/

.OUTPUT_DIRS <- c("Cleaned", "Aligned", "Cufflinks", "merged_asm", "diff_out", "rtpstar_inputs")
.FQ_EXT <- "\\.([Ff][Aa][Ss][Tt][Qq]|[Ff][Qq])(\\.[Gg][Zz])?$"

# Resolve a tool path: <tool_dir>/<rel> if tool_dir is given, else bare name on PATH.
.tool <- function(tool_dir, rel) {
  if (is.null(tool_dir)) basename(rel) else file.path(tool_dir, rel)
}

# Run a command; stop with a clear message on non-zero exit.
.run <- function(cmd, args, step) {
  message(step, "...")
  status <- system2(cmd, args = args, stdout = "", stderr = "")
  if (status != 0) .tool_failed(step, status)
  invisible(status)
}
.tool_failed <- function(step, status)
  stop(step, " failed (exit status ", status, ").", if (status %in% c(137, 139)) " The process was killed, most likely because the computer ran out of memory: use fewer threads (Experimental Design page) or a smaller genome." else "", call. = FALSE)

# Run a shell pipeline (bash, pipefail) so a failure in ANY stage is an error.
.run_pipe <- function(cmdline, step) {
  message(step, "...")
  status <- system2("bash", c("-o", "pipefail", "-c", shQuote(cmdline)), stdout = "", stderr = "")
  if (status != 0) .tool_failed(step, status)
  invisible(status)
}

# ---- space-safe paths ------------------------------------------------------
# Some tools (notably the hisat2 launcher script) re-split their arguments on whitespace, so a
# path like "/mnt/d/TEST FOR TUX/x.fastq.gz" breaks them even when quoted. Every path we hand to a
# tool therefore goes through .safe(): unchanged if it has only ordinary characters, otherwise a
# symlink with a plain name in a per-run temp dir. The links point at the real files/folders
# (which are never moved or modified) and are removed again by .cleanup_links().
.links <- new.env()
.links$dir <- NULL; .links$made <- character(); .links$n <- 0L

.safe <- function(p) {
  if (is.null(p) || is.na(p) || !nzchar(p) || !grepl("[^A-Za-z0-9_./-]", p)) return(p)
  if (is.null(.links$dir)) {
    .links$dir <- file.path(tempdir(), "tuxnet_links"); dir.create(.links$dir, showWarnings = FALSE)
  }
  .links$n <- .links$n + 1L
  nm <- gsub("[^A-Za-z0-9_.-]", "_", basename(p))
  link <- file.path(.links$dir, sprintf("l%03d_%s", .links$n, nm))
  if (!file.symlink(normalizePath(p), link)) stop("Could not create a link for path: ", p)
  .links$made <- c(.links$made, link)
  link
}

# Remove only the links we made (rm -f never follows a symlink), then the temp dir if empty.
.cleanup_links <- function() {
  for (l in .links$made) system2("rm", c("-f", shQuote(l)))
  if (!is.null(.links$dir)) system2("rmdir", c(shQuote(.links$dir)), stderr = FALSE)
  .links$made <- character(); .links$dir <- NULL
}

TOOL_PATHS <- c("ea-utils/clipper/fastq-mcf", "hisat2/hisat2", "hisat2/hisat2-build",
                "samtools/samtools", "cufflinks/cufflinks", "cufflinks/cuffmerge",
                "cufflinks/cuffdiff")

# Read cleaning: "fastq-mcf" (ea-utils, what TuxNet uses; Linux and Windows) or "fastp" (also on macOS). "auto" = fastq-mcf when it is installed, else fastp.
.find_tool <- function(tool_dir, rel) {
  if (is.null(tool_dir)) return(unname(Sys.which(basename(rel))))
  if (file.exists(file.path(tool_dir, rel))) file.path(tool_dir, rel) else unname(Sys.which(basename(rel)))
}
pick_cleaner <- function(cleaner = "auto", tool_dir = NULL) {
  cleaner <- match.arg(cleaner, c("auto", "fastq-mcf", "fastp"))
  if (cleaner != "auto") return(cleaner)
  if (nzchar(.find_tool(tool_dir, "ea-utils/clipper/fastq-mcf"))) "fastq-mcf"
  else if (nzchar(.find_tool(tool_dir, "fastp"))) "fastp" else "fastq-mcf"
}

# Returns a data.frame(tool, found, path). With tool_dir = NULL, searches PATH. The list holds the read cleaner that will be used (see pick_cleaner).
check_tools <- function(tool_dir = NULL, cleaner = "auto") {
  cl <- pick_cleaner(cleaner, tool_dir)
  rels <- if (cl == "fastp") c("fastp", TOOL_PATHS[-1]) else TOOL_PATHS
  path <- vapply(rels, function(rel) {
    if (rel == "fastp") .find_tool(tool_dir, rel)
    else if (is.null(tool_dir)) unname(Sys.which(basename(rel)))
    else if (file.exists(file.path(tool_dir, rel))) file.path(tool_dir, rel) else ""
  }, character(1))
  data.frame(tool = basename(rels), found = nzchar(path), path = unname(path),
             stringsAsFactors = FALSE)
}

# ---- samples -------------------------------------------------------------
.stem <- function(f) sub(.FQ_EXT, "", basename(f))

# For each file: is it a first mate (A_, _R1, _1) and where is its mate?
.mate_info <- function(f) {
  b <- basename(f); d <- dirname(f)
  mate <- rep(NA_character_, length(f)); first <- rep(FALSE, length(f))
  a <- grepl("^A_", b) & file.exists(file.path(d, paste0("B_", substring(b, 3))))          # a name that merely starts with A_ (A_thaliana_...) is not a mate
  mate[a] <- file.path(d[a], paste0("B_", substring(b[a], 3))); first[a] <- TRUE
  r <- !a & grepl("[_.]R?1(_001)?\\.(fastq|fq)(\\.gz)?$", b, ignore.case = TRUE)
  mate[r] <- file.path(d[r], sub("([_.]R?)1((_001)?\\.(fastq|fq)(\\.gz)?)$", "\\12\\2", b[r], ignore.case = TRUE))
  first[r] <- TRUE
  data.frame(first = first, mate = mate, stringsAsFactors = FALSE)
}

# "C2B1_Fe" -> "C2B";  "WT_rep2_R1" -> "WT";  "E4_Fe" -> "E"
.infer_group <- function(stem) {
  s <- sub("^[AB]_", "", sub("_R?[12](_001)?$", "", stem))
  s <- sub("[_.].*$", "", s)
  g <- sub("[-_. ]+$", "", sub("[0-9]+$", "", s))                                     # Control-1 -> Control (not "Control-")
  ifelse(nzchar(g), g, s)
}

.list_fastq <- function(d, paired) {
  f <- sort(list.files(d, pattern = .FQ_EXT, full.names = TRUE))
  skip <- grepl("^Undetermined", basename(f), ignore.case = TRUE)                            # reads that belong to no sample
  if (any(grepl("[_.]R[12](_001)?\\.", basename(f), ignore.case = TRUE))) skip <- skip | grepl("[_.]I[12](_001)?\\.", basename(f), ignore.case = TRUE)   # index reads of a bcl2fastq run
  if (any(skip)) { message("Note: not used as samples: ", paste(utils::head(basename(f)[skip], 6), collapse = ", "), if (sum(skip) > 6) ", ..." else ""); f <- f[!skip] }
  if (!length(f)) return(data.frame(file = character(), mate = character(), stringsAsFactors = FALSE))
  if (!paired) return(data.frame(file = f, mate = NA_character_, stringsAsFactors = FALSE))
  mi <- .mate_info(f)
  data.frame(file = f[mi$first], mate = mi$mate[mi$first], stringsAsFactors = FALSE)
}

# One unique name per sample. The file name without its extension is the name, but the same file name can occur in several condition folders
# (ctrl/rep1.fastq.gz and heat/rep1.fastq.gz): those samples get their group in front, so their counts and conditions can never be mixed up.
.unique_sample_names <- function(stem, group) {
  dup <- stem %in% stem[duplicated(stem)]
  stem[dup] <- paste(group[dup], stem[dup], sep = "_")
  if (anyDuplicated(stem)) stem <- make.unique(stem, sep = "_")
  stem
}

# Clean a sample sheet: trim/sanitize group names, order groups by first appearance,
# number replicates within each group. Returns data.frame(group, k, file, mate, sample).
normalize_samples <- function(samples, paired) {
  samples <- as.data.frame(samples, stringsAsFactors = FALSE)
  if (!all(c("group", "file") %in% names(samples)))
    stop("Sample sheet needs columns 'group' and 'file'.")
  if (!"mate" %in% names(samples)) samples$mate <- rep(NA_character_, nrow(samples))
  samples <- samples[nzchar(trimws(samples$group)) & nzchar(trimws(samples$file)), , drop = FALSE]
  if (!nrow(samples)) return(data.frame(group = character(), k = integer(), file = character(),
                                        mate = character(), sample = character()))
  samples$group <- gsub("[^A-Za-z0-9._-]+", "_", trimws(samples$group))
  samples$file <- trimws(samples$file)
  samples$group <- factor(samples$group, levels = unique(samples$group))
  samples <- samples[order(samples$group), , drop = FALSE]
  samples$k <- as.integer(ave(seq_len(nrow(samples)), samples$group, FUN = seq_along))
  samples$group <- as.character(samples$group)
  samples$sample <- .unique_sample_names(.stem(samples$file), samples$group)
  rownames(samples) <- NULL
  samples[, c("group", "k", "file", "mate", "sample")]
}

# Build a sample sheet from an inputs folder (flat folder or one sub-folder per group).
infer_samples <- function(inputs_folder, paired = FALSE) {
  flat <- .list_fastq(inputs_folder, paired)
  subs <- setdiff(list.dirs(inputs_folder, full.names = FALSE, recursive = FALSE), c(.OUTPUT_DIRS, "tuxedo_results"))
  rows <- lapply(subs, function(g) {                                               # one sub-folder per condition; files at any depth below it belong to that condition
    x <- .list_fastq(file.path(inputs_folder, g), paired)
    if (!nrow(x)) {                                                                  # deeper layouts (condition/sample/reads.fastq.gz): look below
      deep <- list.dirs(file.path(inputs_folder, g), full.names = TRUE, recursive = TRUE)[-1]
      x <- do.call(rbind, lapply(deep, function(dd) .list_fastq(dd, paired)))
    }
    if (!is.null(x) && nrow(x)) cbind(group = g, x, stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  if (nrow(flat) && !is.null(rows) && nrow(rows)) {                                # both: sub-folders are the design, the loose files are said to be ignored
    message("Note: files directly in the folder are not used because the condition folders hold the samples: ", paste(utils::head(basename(flat$file), 4), collapse = ", "), if (nrow(flat) > 4) ", ..." else "")
    return(normalize_samples(rows, paired))
  }
  if (nrow(flat)) {
    flat$group <- .infer_group(.stem(flat$file))
    return(normalize_samples(flat, paired))
  }
  if (is.null(rows)) return(normalize_samples(data.frame(group = character(), file = character()), paired))
  normalize_samples(rows, paired)
}

# Cheap FASTQ sanity check (first record only; fast even on huge files). Results are cached
# per file (path + size + mtime) because the GUI re-validates on every edit.
.fq_cache <- new.env()
.fastq_problem <- function(path) {
  key <- paste(path, file.size(path), format(file.mtime(path), "%s"))
  if (!is.na(file.size(path)) && exists(key, envir = .fq_cache, inherits = FALSE))
    return(get(key, envir = .fq_cache))
  res <- .fastq_problem_uncached(path)
  if (!is.na(file.size(path))) assign(key, res, envir = .fq_cache)
  res
}
.fastq_problem_uncached <- function(path) {
  con <- tryCatch(if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt"),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(con)) return("cannot be opened (is the drive still connected?)")
  on.exit(close(con))
  l <- tryCatch(readLines(con, n = 4), error = function(e) character())
  if (length(l) < 4) return("file is empty, truncated or unreadable")
  if (!startsWith(l[1], "@") || !startsWith(l[3], "+") || nchar(l[2]) != nchar(l[4]))
    return("does not look like FASTQ")
  ""
}

# The first reads of a FASTQ(.gz): names, lengths and quality characters (a few hundred reads tell short from long reads, Phred+33 from +64,
# and whether two mate files belong together).
.fq_head <- function(path, n = 100) {
  key <- paste("head", path, file.size(path), format(file.mtime(path), "%s"))
  if (exists(key, envir = .fq_cache, inherits = FALSE)) return(get(key, envir = .fq_cache))
  res <- list(ids = character(), len = integer(), qmin = NA_integer_, qmax = NA_integer_)
  con <- tryCatch(if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt"), error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(con)) {
    l <- tryCatch(readLines(con, n = 4 * n, warn = FALSE), error = function(e) character(), warning = function(w) character()); close(con)
    if (length(l) >= 4) { k <- seq(1, length(l) - 3, by = 4)
      res$ids <- sub("[/ ].*$", "", sub("^@", "", l[k])); res$len <- nchar(l[k + 1])
      q <- tryCatch(unlist(lapply(l[k + 3], utf8ToInt)), error = function(e) integer()); if (length(q)) { res$qmin <- min(q); res$qmax <- max(q) } }
  }
  assign(key, res, envir = .fq_cache); res
}
# fastq-mcf guesses the quality offset from the data, and guesses +64 for a small file whose bases are all high quality (every read is then
# thrown away). Modern data is Phred+33; +64 only when the characters say so (none below ";" and some above "J").
.detect_phred <- function(path) { h <- .fq_head(path); if (!is.na(h$qmin) && h$qmin >= 64 && h$qmax >= 100) 64L else 33L }       # +64: nothing below "@" and values up to "h" (Q40); a "K" is a normal +33 base
# Names of the compressed files that are damaged or cut short (a copy that was interrupted): their first reads look fine, the rest is missing.
check_gz_integrity <- function(files) {
  files <- unique(files[!is.na(files) & grepl("\\.gz$", files, ignore.case = TRUE) & file.exists(files)]); bad <- character()
  for (f in files) { st <- suppressWarnings(system2("gzip", c("-t", shQuote(f)), stdout = FALSE, stderr = FALSE)); if (!identical(as.integer(st), 0L)) bad <- c(bad, basename(f)) }
  bad
}

# Validate a (normalized) sample sheet. Returns list(errors, warnings).
validate_samples <- function(samples, paired, platform = "short") {
  errs <- character(); warns <- character()
  if (!nrow(samples))
    return(list(errors = paste0("No FASTQ files found (files must end in .fastq/.fq, ",
                                "optionally .gz).", if (paired)
                                " Paired-end mates must be named A_x/B_x, x_R1/x_R2 or x_1/x_2." else ""),
                warnings = warns))
  ng <- length(unique(samples$group))
  if (ng < 2) errs <- c(errs, sprintf(
    "Need at least 2 groups (conditions) to compare; found %d. Edit the Group column.", ng))
  single <- names(which(table(samples$group) == 1))
  if (length(single)) warns <- c(warns, paste0(
    "Only one replicate in group(s): ", paste(single, collapse = ", "),
    ". Cuffdiff statistics are weak without replicates."))
  if (anyDuplicated(samples$file)) errs <- c(errs, "The same file is listed more than once.")
  if (!paired) {                                                                                   # single-end chosen, but the files come in pairs
    mi <- .mate_info(samples$file); pr <- mi$first & !is.na(mi$mate) & file.exists(mi$mate)
    if (any(pr)) warns <- c(warns, sprintf("%d of these files have a matching mate (for example %s and %s): this looks like paired-end data. Choose Paired-end as the read type, otherwise every mate is counted as a separate sample (if these are single-end replicates that happen to end in _1 and _2, ignore this).", sum(pr), basename(samples$file[pr][1]), basename(mi$mate[pr][1])))
  }
  lane <- grepl("_L00[0-9]", basename(samples$file))
  if (any(lane) && anyDuplicated(sub("_L00[0-9]", "", .stem(samples$file)[lane])))
    warns <- c(warns, "Several files look like one sample split over sequencing lanes (_L001, _L002 ...): each file would be counted as its own replicate. Join the lanes of each sample first, for example  cat S1_L001_R1_001.fastq.gz S1_L002_R1_001.fastq.gz > S1_R1.fastq.gz  (compressed files can simply be joined this way), and use the joined files.")
  ex <- utils::head(samples$file[file.exists(samples$file)], 6)
  if (length(ex)) { len <- unlist(lapply(ex, function(x) .fq_head(x, 50)$len))
    if (length(len) && identical(platform, "short") && stats::median(len) > 600)
      errs <- c(errs, sprintf("These reads are long (median %s bases): this looks like Oxford Nanopore or PacBio data, for example from Plasmidsaurus. The short-read pipeline (HISAT2, Cufflinks) would throw these reads away or misalign them. Choose \u201cLong reads\u201d under Sequencing platform (reads are then aligned with minimap2 and counted with featureCounts), or use the counts table your provider supplies (Experimental Design: \u201cA counts table\u201d).", format(round(stats::median(len)), big.mark = ",")))
    if (length(len) && identical(platform, "long") && stats::median(len) < 300)
      warns <- c(warns, sprintf("These reads are short (median %s bases) although the sequencing platform is set to Long reads: for Illumina-type reads choose Short reads.", format(round(stats::median(len)), big.mark = ","))) }
  for (i in seq_len(nrow(samples))) {
    f <- samples$file[i]
    if (!file.exists(f)) { errs <- c(errs, sprintf("%s: file not found (drive disconnected?).", basename(f))); next }
    p <- .fastq_problem(f)
    if (nzchar(p)) errs <- c(errs, sprintf("%s: %s", basename(f), p))
    if (paired) {
      m <- samples$mate[i]
      if (is.na(m) || !file.exists(m))
        errs <- c(errs, sprintf("%s: mate file %s is missing.", basename(f),
                                if (is.na(m)) "(none)" else basename(m)))
      else {                                                                                       # the mate must be a readable FASTQ that belongs to this file
        pm <- .fastq_problem(m)
        if (nzchar(pm)) errs <- c(errs, sprintf("%s: %s", basename(m), pm))
        else { h1 <- .fq_head(f, 50)$ids; h2 <- .fq_head(m, 50)$ids; n <- min(length(h1), length(h2))
               if (n > 0 && !identical(h1[seq_len(n)], h2[seq_len(n)]))
                 warns <- c(warns, sprintf("%s and %s may not belong together: their read names differ (re-sorted, from different runs, or wrongly paired). If you know they do, you can run anyway.", basename(f), basename(m))) }
      }
    }
  }
  list(errors = errs, warnings = warns)
}

# ---- disk space ----------------------------------------------------------
# Free bytes on the filesystem holding `path` (NA if it can't be determined).
free_space_bytes <- function(path) {
  while (!file.exists(path) && dirname(path) != path) path <- dirname(path)
  out <- tryCatch(suppressWarnings(system2("df", c("-Pk", shQuote(path)), stdout = TRUE, stderr = FALSE)),
                  error = function(e) character())
  if (length(out) < 2) return(NA_real_)
  as.numeric(strsplit(trimws(out[length(out)]), "\\s+")[[1]][4]) * 1024
}

# Rough disk requirement for the run, in bytes. Uncompressed-FASTQ-equivalent size F;
# gz files count as 4x. Peak ~= sorted BAMs for all samples (~0.4 F) + Cufflinks/
# Cuffdiff output + working set for one sample (cleaned reads ~1 F_i). Padded.
estimate_space <- function(samples, paired) {
  one <- function(f) { ok <- !is.na(f) & file.exists(f); out <- rep(0, length(f)); out[ok] <- file.size(f[ok]) * ifelse(grepl("\\.gz$", f[ok]), 4, 1); out }
  per_sample <- one(samples$file) + if (paired) one(samples$mate) else 0              # each sample on its own: a missing file counts as 0, nothing is shifted
  per_sample <- per_sample[is.finite(per_sample)]
  if (!length(per_sample) || sum(per_sample) == 0) return(list(total = 0, largest_sample = 0, need = 0))
  total <- sum(per_sample); big <- max(per_sample)
  list(total = total, largest_sample = big, need = 0.5 * total + 1.5 * big + 1e9)
}

fmt_bytes <- function(b) {
  if (is.null(b) || is.na(b)) return("unknown")
  u <- c("B", "KB", "MB", "GB", "TB"); i <- if (b > 0) min(4, floor(log(b, 1024))) else 0
  sprintf("%.1f %s", b / 1024^i, u[i + 1])
}

# Check the output folder: must be creatable/writable and free of earlier outputs.
# Returns list(errors, warnings).
check_output_folder <- function(out_folder) {
  errs <- character(); warns <- character()
  if (!nzchar(out_folder)) return(list(errors = "Choose an output folder.", warnings = character()))
  if (dir.exists(out_folder)) {
    existing <- .OUTPUT_DIRS[dir.exists(file.path(out_folder, .OUTPUT_DIRS))]
    if (length(existing))
      warns <- c(warns, paste0("Previous outputs exist in the output folder (",
                               paste(existing, collapse = ", "),
                               "). Delete them or choose a different output folder."))
    probe <- file.path(out_folder, paste0(".tux_write_test_", Sys.getpid()))
    if (!isTRUE(tryCatch(suppressWarnings(file.create(probe)), error = function(e) FALSE)))
      errs <- c(errs, "Output folder is not writable.")
    else unlink(probe)
  } else if (!dir.exists(dirname(out_folder))) {
    errs <- c(errs, "Output folder's parent does not exist (is the drive mounted?).")
  }
  list(errors = errs, warnings = warns)
}

.bam_name <- function(g, k) sprintf("hisat2OutSorted%s_%d.bam", g, k)

# ---- resuming an interrupted run -----------------------------------------
# A run can be cut short (a reboot, a closed laptop, a crash). tux_resume_status() looks at what an earlier run left in `out_folder` and
# says which parts can be reused. A part is reused ONLY when it is demonstrably intact and made from the same inputs:
#   * the same samples (same FASTQ files, groups and replicate numbers) as the earlier run, and the same cleaning settings and reference
#     (recorded in run_params.json; a run from before that file existed counts as made with the default cleaning settings);
#   * alignment:  the BAM passes `samtools quickcheck` and is newer than its FASTQ file(s) and the genome;
#   * Cufflinks:  its transcripts.gtf / genes.fpkm_tracking exist, are not empty, and are newer than the BAM and the reference;
#   * Cuffmerge:  merged.gtf is newer than every assembly;     * Cuffdiff:  its result tables exist and are newer than all of the above.
# Each part also needs the parts before it, so anything uncertain is simply redone. Returns list(align, cufflinks, merge, cuffdiff, why).
tux_resume_status <- function(samples, out_folder, gtf_file, fasta_file, quality_thresh = 30, min_seq_length = 30, window_size = 5,
                              paired = FALSE, tool_dir = NULL, cleaner = "fastq-mcf") {
  n <- nrow(samples)
  none <- function(why = "") list(align = rep(FALSE, n), cufflinks = rep(FALSE, n), merge = FALSE, cuffdiff = FALSE, why = why)
  if (!dir.exists(out_folder)) return(none("no earlier run found"))
  ts <- function(f) { m <- suppressWarnings(as.numeric(file.mtime(f))); if (length(m) != 1) NA_real_ else m }
  ok_file <- function(f) file.exists(f) && isTRUE(file.size(f) > 0)
  newest <- function(x) { x <- x[!is.na(x)]; if (length(x)) max(x) else -Inf }
  prev <- file.path(out_folder, "samples_used.csv")
  if (!ok_file(prev)) return(none("no earlier run found"))
  old <- tryCatch(utils::read.csv(prev, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(old) || !all(c("group", "k", "file") %in% names(old))) return(none("the earlier run's sample list is unreadable"))
  key <- function(g, k, f) paste(g, k, normalizePath(f, mustWork = FALSE))
  if (!setequal(key(old$group, old$k, old$file), key(samples$group, samples$k, samples$file))) return(none("the samples differ from the earlier run"))
  now <- c(q = quality_thresh, l = min_seq_length, w = window_size, paired = as.integer(paired), gtf = file.size(gtf_file), fasta = file.size(fasta_file))
  rec_f <- file.path(out_folder, "run_params.json")
  if (file.exists(rec_f) && requireNamespace("jsonlite", quietly = TRUE)) {
    rec <- tryCatch(jsonlite::fromJSON(rec_f), error = function(e) NULL)
    if (is.null(rec) || !identical(as.character(unlist(rec[names(now)])), as.character(now))) return(none("the cleaning settings or the reference differ from the earlier run"))
    if (!identical(if (is.null(rec$cl)) "fastq-mcf" else as.character(rec$cl), cleaner)) return(none("the read cleaner differs from the earlier run"))
  } else if (!(quality_thresh == 30 && min_seq_length == 30 && window_size == 5))
    return(none("the cleaning settings differ from the defaults used by earlier runs"))

  ref_t <- newest(c(ts(gtf_file), ts(fasta_file)))
  samtools <- .tool(tool_dir, "samtools/samtools")
  bam_ok <- function(bam) ok_file(bam) && isTRUE(suppressWarnings(system2(samtools, c("quickcheck", shQuote(bam)), stdout = FALSE, stderr = FALSE)) == 0)
  bams <- file.path(out_folder, "Aligned", .bam_name(samples$group, samples$k))
  cdirs <- file.path(out_folder, "Cufflinks", sprintf("CufflinksOut%s_%d", samples$group, samples$k))
  align <- cuff <- logical(n)
  for (i in seq_len(n)) {
    inputs_t <- newest(c(ts(samples$file[i]), if (paired) ts(samples$mate[i]), ts(fasta_file)))
    align[i] <- bam_ok(bams[i]) && isTRUE(ts(bams[i]) >= inputs_t)
    tg <- file.path(cdirs[i], "transcripts.gtf")
    cuff[i] <- align[i] && ok_file(tg) && ok_file(file.path(cdirs[i], "genes.fpkm_tracking")) && isTRUE(ts(tg) >= newest(c(ts(bams[i]), ref_t)))
  }
  mf <- file.path(out_folder, "merged_asm", "merged.gtf"); asm <- file.path(out_folder, "Cufflinks", "assemblies.txt")
  merge <- all(cuff) && ok_file(mf) && ok_file(asm) && isTRUE(ts(mf) >= newest(c(vapply(file.path(cdirs, "transcripts.gtf"), ts, 1), ts(asm))))
  ge <- file.path(out_folder, "diff_out", "gene_exp.diff"); rg <- file.path(out_folder, "diff_out", "genes.read_group_tracking")
  cd <- merge && ok_file(ge) && ok_file(rg) && isTRUE(ts(ge) >= newest(c(ts(mf), vapply(bams, ts, 1))))
  list(align = align, cufflinks = cuff, merge = merge, cuffdiff = cd, why = "")
}

# One plain-language line for a resume status (shown in the run log and in the confirmation window).
tux_resume_text <- function(st) {
  n <- length(st$align)
  if (!any(st$align) && !any(st$cufflinks) && !st$merge && !st$cuffdiff) return(if (nzchar(st$why)) paste0("Nothing to reuse (", st$why, "): starting from the beginning.") else "Nothing to reuse: starting from the beginning.")
  done <- c(sprintf("alignments %d/%d", sum(st$align), n), sprintf("transcript assembly %d/%d", sum(st$cufflinks), n), if (st$merge) "merge", if (st$cuffdiff) "differential expression")
  todo <- c(if (!all(st$align)) "alignment", if (!all(st$cufflinks)) "assembly", if (!st$merge) "merge", if (!st$cuffdiff) "Cuffdiff")
  sprintf("Reusing finished work from an earlier run (%s). Still to do: %s.", paste(done, collapse = ", "), if (length(todo)) paste(todo, collapse = ", ") else "nothing (the run is complete)")
}

# fastp counterpart of the fastq-mcf call: -q / -l / -w of TuxNet become a sliding-window quality trim from the 3' end (window -w, mean quality -q),
# a minimum length (-l) and the same adapter file. Pairs stay in step (a read whose mate is dropped is dropped too). The reports go to a temporary
# folder. fastp reads Phred+33; Phred+64 files are converted (--phred64).
fastp_args <- function(f, out, adapter_fasta, q, l, w, threads, phred = 33L, mate = NULL, mate_out = NULL) {
  rep_dir <- tempfile("fastp_report"); dir.create(rep_dir)
  a <- c("-i", shQuote(f), "-o", shQuote(out))
  if (!is.null(mate)) a <- c(a, "-I", shQuote(mate), "-O", shQuote(mate_out))
  c(a, "--adapter_fasta", shQuote(adapter_fasta), "--cut_tail", "--cut_window_size", min(max(as.integer(w), 1L), 1000L), "--cut_mean_quality", q,
    "--length_required", l, "--thread", min(max(as.integer(threads), 1L), 16L), if (phred == 64L) "--phred64",
    "--json", shQuote(file.path(rep_dir, "fastp.json")), "--html", shQuote(file.path(rep_dir, "fastp.html")))
}

# ---- 1. Clean + align, one sample at a time ------------------------------
clean_and_align <- function(samples, out_folder, fasta_file, adapter_fasta,
                            quality_thresh, min_seq_length, window_size,
                            threads, paired, tool_dir, keep_cleaned = FALSE, skip = rep(FALSE, nrow(samples)), cleaner = "auto") {
  cleaner <- pick_cleaner(cleaner, tool_dir)
  cleaned_dir <- file.path(out_folder, "Cleaned"); dir.create(cleaned_dir, showWarnings = FALSE)
  aligned_dir <- file.path(out_folder, "Aligned"); dir.create(aligned_dir, showWarnings = FALSE)
  fastq_mcf    <- if (cleaner == "fastp") .find_tool(tool_dir, "fastp") else .tool(tool_dir, "ea-utils/clipper/fastq-mcf")   # (variable name kept: it is "the cleaner")
  message("Read cleaning with ", cleaner)
  hisat2       <- .tool(tool_dir, "hisat2/hisat2")
  hisat2_build <- .tool(tool_dir, "hisat2/hisat2-build")
  samtools     <- .tool(tool_dir, "samtools/samtools")

  index <- file.path(aligned_dir, "index")
  if (!all(skip)) .run(hisat2_build, c("-p", threads, shQuote(fasta_file), shQuote(index)), "hisat2-build")   # the index is only needed to align

  common0 <- c("-q", quality_thresh, "-l", min_seq_length, "-w", window_size)
  sort_threads <- max(1L, min(as.integer(threads), 4L))                           # samtools sort keeps -m per thread in memory: capped so a long thread list cannot exhaust the RAM
  for (i in seq_len(nrow(samples))) {
    common <- c(common0, "-P", .detect_phred(samples$file[i]))
    g <- samples$group[i]; k <- samples$k[i]; f <- samples$file_run[i]
    message(sprintf("[SAMPLE %d/%d] %s (%s)", i, nrow(samples), samples$sample[i], g))
    if (skip[i]) { message("Alignment already done: reusing it"); next }
    if (paired) {
      a_out <- file.path(cleaned_dir, sprintf("A_Cleaned%s_%d.fastq", g, k))
      b_out <- file.path(cleaned_dir, sprintf("B_Cleaned%s_%d.fastq", g, k))
      if (cleaner == "fastp")
        .run(fastq_mcf, fastp_args(f, a_out, adapter_fasta, quality_thresh, min_seq_length, window_size, threads, .detect_phred(samples$file[i]), mate = samples$mate_run[i], mate_out = b_out),
             paste("Cleaning", basename(f)))
      else .run(fastq_mcf, c(shQuote(adapter_fasta), shQuote(f), shQuote(samples$mate_run[i]), common,
                        "-o", shQuote(a_out), "-o", shQuote(b_out)),
           paste("Cleaning", basename(f)))
      reads <- paste("-1", shQuote(a_out), "-2", shQuote(b_out)); cleaned <- c(a_out, b_out)
    } else {
      a_out <- file.path(cleaned_dir, sprintf("Cleaned%s_%d.fastq", g, k))
      if (cleaner == "fastp")
        .run(fastq_mcf, fastp_args(f, a_out, adapter_fasta, quality_thresh, min_seq_length, window_size, threads, .detect_phred(samples$file[i])),
             paste("Cleaning", basename(f)))
      else .run(fastq_mcf, c(shQuote(adapter_fasta), shQuote(f), common, "-o", shQuote(a_out)),
           paste("Cleaning", basename(f)))
      reads <- paste("-U", shQuote(a_out)); cleaned <- a_out
    }
    bam <- file.path(aligned_dir, .bam_name(g, k))
    tmp <- file.path(aligned_dir, sprintf(".sorttmp_%s_%d", g, k))
    cmd <- sprintf("%s --dta-cufflinks -p %s -x %s %s | %s sort -@ %s -m 400M -T %s -O bam -o %s -",
                   shQuote(hisat2), threads, shQuote(index), reads,
                   shQuote(samtools), sort_threads, shQuote(tmp), shQuote(bam))
    .run_pipe(cmd, sprintf("hisat2 + samtools sort %s", basename(f)))
    unlink(Sys.glob(paste0(tmp, "*")))
    # a sample that lost all its reads, or none of which aligned, is caught HERE, not after hours of Cufflinks
    fs <- suppressWarnings(system2(samtools, c("flagstat", shQuote(bam)), stdout = TRUE, stderr = FALSE))
    tot <- suppressWarnings(as.numeric(sub(" .*", "", fs[1]))); pct <- suppressWarnings(as.numeric(sub(".*\\(([0-9.]+)%.*", "\\1", grep(" mapped \\(", fs, value = TRUE)[1])))
    if (!is.na(tot) && tot == 0)
      stop(sprintf("No reads of %s were left after cleaning: the quality, length and adapter filters removed all of them (reads shorter than the minimum length, long reads, a wrong quality setting?).", samples$sample[i]), call. = FALSE)
    if (!is.na(pct)) {
      message(sprintf("%s: %.1f%% of the reads aligned", samples$sample[i], pct))
      if (pct == 0) stop(sprintf("None of the reads of %s aligned to the genome. Wrong species or genome file, or not RNA-seq reads?", samples$sample[i]), call. = FALSE)
      if (pct < 30) message(sprintf("WARNING: only %.1f%% of the reads of %s aligned (wrong genome or annotation for this species, or a contaminated sample?).", pct, samples$sample[i]))
    }
    if (!keep_cleaned) unlink(cleaned)                # free space immediately
  }
  if (!keep_cleaned) unlink(cleaned_dir, recursive = TRUE)
  message("Cleaning and alignment successful on all files")
}

# ---- 2. Cufflinks --------------------------------------------------------
# Look INSIDE the annotation / genome files so a wrong file (e.g. a sample sheet) is caught before hours of alignment.
# Returns "" if the file looks right, otherwise a plain-language problem.
annotation_problem <- function(path) {
  if (grepl("\\.(gz|zip)$", path, ignore.case = TRUE)) return("The annotation must be uncompressed (unzip it first).")
  ln <- tryCatch(readLines(path, n = 200, warn = FALSE), error = function(e) character())
  ln <- ln[nzchar(ln) & !startsWith(ln, "#")]
  if (!length(ln)) return("The annotation file is empty (or only comments).")
  f <- strsplit(ln[1], "\t", fixed = TRUE)[[1]]
  if (length(f) < 9 || is.na(suppressWarnings(as.integer(f[4]))) || is.na(suppressWarnings(as.integer(f[5]))))
    return(sprintf("This does not look like a GTF/GFF annotation (expected 9 tab-separated columns with start/end positions). Is it the right file? It starts with: %s", substr(ln[1], 1, 60)))
  ""
}
genome_problem <- function(path) {
  ln <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character())
  if (!length(ln) || !startsWith(ln[1], ">")) return("This does not look like a genome FASTA (the first line should start with '>').")
  ""
}

run_cufflinks <- function(samples, out_folder, gtf_file, fasta_file, threads, tool_dir, skip = rep(FALSE, nrow(samples))) {
  cl_dir <- file.path(out_folder, "Cufflinks"); dir.create(cl_dir, showWarnings = FALSE)
  cufflinks <- .tool(tool_dir, "cufflinks/cufflinks")
  assemblies <- character()
  for (i in seq_len(nrow(samples))) {
    g <- samples$group[i]; k <- samples$k[i]
    message(sprintf("[SAMPLE %d/%d] %s (%s)", i, nrow(samples), samples$sample[i], g))
    out <- file.path(cl_dir, sprintf("CufflinksOut%s_%d", g, k))
    if (skip[i]) message("Assembly already done: reusing it")
    else .run(cufflinks, c("-q", "--no-update-check", "-p", threads, "-G", shQuote(gtf_file), "-b", shQuote(fasta_file),
                           "-o", shQuote(out), shQuote(file.path(out_folder, "Aligned", .bam_name(g, k)))),
              sprintf("Cufflinks %s", samples$sample[i]))
    assemblies <- c(assemblies, file.path(out, "transcripts.gtf"))
  }
  writeLines(assemblies, file.path(cl_dir, "assemblies.txt"))
  message("Cufflinks successful on all files")
}

# ---- 3. Cuffmerge --------------------------------------------------------
run_cuffmerge <- function(out_folder, gtf_file, fasta_file, threads, tool_dir) {
  cuffmerge <- .tool(tool_dir, "cufflinks/cuffmerge")
  work <- tempfile(".cuffmerge-", tmpdir = out_folder); dir.create(work)
  old_wd <- setwd(work)
  on.exit({ setwd(old_wd); unlink(work, recursive = TRUE) }, add = TRUE)
  .run(cuffmerge, c("-g", shQuote(gtf_file), "-s", shQuote(fasta_file), "-p", threads,
                    "-o", shQuote(file.path(out_folder, "merged_asm")),
                    shQuote(file.path(out_folder, "Cufflinks", "assemblies.txt"))),
       "Cuffmerge")
  message("Cuffmerge successful")
}

# ---- 4. Cuffdiff ---------------------------------------------------------
run_cuffdiff <- function(samples, out_folder, fasta_file, threads, tool_dir) {
  cuffdiff <- .tool(tool_dir, "cufflinks/cuffdiff")
  groups <- unique(samples$group)
  # One comma-separated list of replicate BAMs per group, groups as separate args.
  bam_groups <- vapply(groups, function(g) {
    k <- samples$k[samples$group == g]
    paste(shQuote(file.path(out_folder, "Aligned", .bam_name(g, k))), collapse = ",")
  }, character(1))
  .run(cuffdiff, c("-q", "--no-update-check", "-o", shQuote(file.path(out_folder, "diff_out")),
                   "-b", shQuote(fasta_file), "-p", threads,
                   "-L", paste(groups, collapse = ","),
                   "-u", shQuote(file.path(out_folder, "merged_asm", "merged.gtf")),
                   unname(bam_groups)),
       "Cuffdiff")
  message("Cuffdiff successful")
}

# ---- Full pipeline (equivalent of pressing "Run" on the TUX tab) --------
#' @param inputs_folder  folder with FASTQ files (flat, or one sub-folder per group); read-only
#' @param gtf_file       reference annotation (.gtf/.gff)
#' @param fasta_file     reference genome (.fa/.fasta)
#' @param output_folder  where all results/intermediates go (default: inputs_folder/tuxedo_results)
#' @param read_type      "single" or "paired"
#' @param samples        optional sample sheet data.frame(group, file[, mate]); default: infer
#' @param cleaner  read cleaning program: \
#' @param min_seq_length fastq-mcf -l (default 30)
#' @param window_size    fastq-mcf -w (default 5)
#' @param threads        threads for all tools (default 1)
#' @param adapter_fasta  adapter file (IlluminaAdaptorSeq.fasta from TuxNet-MATLAB)
#' @param tool_dir       dir holding the tool sub-folders, or NULL to use PATH
#' @param keep_cleaned   keep cleaned FASTQs (default FALSE: deleted after alignment)
#' @param skip_space_check  run even if estimated free space looks too small
#' @param resume         reuse the finished, intact parts of an earlier run in output_folder (see tux_resume_status); default FALSE
run_tuxedo <- function(inputs_folder, gtf_file, fasta_file,
                       output_folder = file.path(inputs_folder, "tuxedo_results"),
                       read_type = c("single", "paired"), samples = NULL,
                       quality_thresh = 30, min_seq_length = 30, window_size = 5,
                       threads = 1,
                       adapter_fasta = "IlluminaAdaptorSeq.fasta",
                       tool_dir = NULL, keep_cleaned = FALSE, skip_space_check = FALSE, resume = FALSE, cleaner = "auto") {
  read_type <- match.arg(read_type)
  cleaner <- pick_cleaner(cleaner, tool_dir)
  paired <- read_type == "paired"
  output_folder <- normalizePath(output_folder, mustWork = FALSE)

  if (!nzchar(gtf_file) || !file.exists(gtf_file))     stop("Annotation (GTF/GFF) file cannot be found.")
  if (!nzchar(fasta_file) || !file.exists(fasta_file)) stop("Genome FASTA file cannot be found.")
  if (!file.exists(adapter_fasta)) stop("Adapter file not found: ", adapter_fasta)
  p <- annotation_problem(gtf_file); if (nzchar(p)) stop(p, call. = FALSE)
  p <- genome_problem(fasta_file);   if (nzchar(p)) stop(p, call. = FALSE)
  if (is.null(samples)) {
    if (!nzchar(inputs_folder) || !dir.exists(inputs_folder)) stop("Inputs folder cannot be found.")
    samples <- infer_samples(normalizePath(inputs_folder), paired)
  } else samples <- normalize_samples(samples, paired)

  oc <- check_output_folder(output_folder)
  if (length(oc$errors) || (!resume && length(oc$warnings))) stop(paste(c(oc$errors, oc$warnings), collapse = " "))   # resuming: earlier outputs are expected
  tools <- check_tools(tool_dir, cleaner)
  if (!all(tools$found)) stop("Missing tools: ", paste(tools$tool[!tools$found], collapse = ", "),
                              if (!cleaner %in% tools$tool[tools$found] && cleaner == "fastq-mcf" && !nzchar(.find_tool(tool_dir, "fastp")))
                                " (install fastp in the tuxnet environment, or use fastq-mcf)" else "")
  v <- validate_samples(samples, paired)
  if (length(v$errors)) stop(paste(c("Input problems:", v$errors), collapse = "\n  "))
  for (w in v$warnings) message("WARNING: ", w)
  message("Checking that the compressed files are complete...")
  bad_gz <- check_gz_integrity(c(samples$file, if (paired) samples$mate))
  if (length(bad_gz)) stop("These compressed files are damaged or incomplete (a copy that was interrupted?): ", paste(bad_gz, collapse = ", "), ". Copy them again.", call. = FALSE)

  message(sprintf("%d samples in %d groups: %s", nrow(samples), length(unique(samples$group)),
                  paste(sprintf("%s(%d)", unique(samples$group), as.integer(table(samples$group)[unique(samples$group)])),
                        collapse = " ")))
  est <- estimate_space(samples, paired)
  free <- free_space_bytes(output_folder)
  message(sprintf("Input size ~%s (uncompressed equivalent); estimated space needed ~%s; free on output drive: %s",
                  fmt_bytes(est$total), fmt_bytes(est$need), fmt_bytes(free)))
  if (!skip_space_check && !is.na(free) && free < est$need)
    stop(sprintf("Not enough free space on the output drive: need about %s, have %s. ",
                 fmt_bytes(est$need), fmt_bytes(free)),
         "Choose an output folder on a bigger drive.")

  st <- if (resume) tux_resume_status(samples, output_folder, gtf_file, fasta_file, quality_thresh, min_seq_length, window_size, paired, tool_dir, cleaner)
        else list(align = rep(FALSE, nrow(samples)), cufflinks = rep(FALSE, nrow(samples)), merge = FALSE, cuffdiff = FALSE, why = "")
  if (resume) message(tux_resume_text(st))
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  write.csv(samples, file.path(output_folder, "samples_used.csv"), row.names = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE))                          # what this run was made with (checked by a later resume)
    jsonlite::write_json(list(q = quality_thresh, l = min_seq_length, w = window_size, paired = as.integer(paired),
                              gtf = file.size(gtf_file), fasta = file.size(fasta_file), cl = cleaner), file.path(output_folder, "run_params.json"), auto_unbox = TRUE)

  # From here on every path given to a tool is space-free (see .safe).
  on.exit(.cleanup_links(), add = TRUE)
  out_run <- .safe(output_folder)
  samples$file_run <- vapply(samples$file, .safe, character(1), USE.NAMES = FALSE)
  samples$mate_run <- if (paired) vapply(samples$mate, .safe, character(1), USE.NAMES = FALSE) else NA_character_
  gtf_run <- .safe(gtf_file); fasta_run <- .safe(fasta_file); adapter_run <- .safe(adapter_fasta)

  message("[STEP 1/4] Cleaning and aligning reads")
  clean_and_align(samples, out_run, fasta_run, adapter_run, quality_thresh,
                  min_seq_length, window_size, threads, paired, tool_dir, keep_cleaned, skip = st$align, cleaner = cleaner)
  message("[STEP 2/4] Assembling transcripts")
  run_cufflinks(samples, out_run, gtf_run, fasta_run, threads, tool_dir, skip = st$cufflinks)
  message("[STEP 3/4] Merging assemblies")
  if (st$merge) message("Cuffmerge already done: reusing it") else { unlink(file.path(output_folder, "merged_asm"), recursive = TRUE); run_cuffmerge(out_run, gtf_run, fasta_run, threads, tool_dir) }
  message("[STEP 4/4] Differential expression")
  if (st$cuffdiff) message("Cuffdiff already done: reusing it")
  else { unlink(file.path(output_folder, "diff_out"), recursive = TRUE); run_cuffdiff(samples, out_run, fasta_run, threads, tool_dir) }   # never mix with a half-written earlier attempt

  tux_write_rtpstar_inputs(output_folder)                       # the genes file and the expression file RTP-STAR needs (never stops the run)
  message("Done. Results in ", file.path(output_folder, "diff_out"))
  invisible(file.path(output_folder, "diff_out"))
}

# Remove the output folders of a previous (possibly failed) run so it can be re-run.
# Only touches the five known output folders inside out_folder.
clear_outputs <- function(out_folder) {
  dirs <- file.path(out_folder, .OUTPUT_DIRS)
  dirs <- dirs[dir.exists(dirs)]
  unlink(dirs, recursive = TRUE)
  invisible(basename(dirs))
}

# ---- tools for long reads (see longreads.R) -----------------------------------------------------------------------------
# The "longread" environment sits next to the main one. Check an explicit StRym root, micromamba's active root,
# the locations used by the Linux/WSL, macOS and Docker installers, and the running R environment's sibling.
.long_env_bins <- function() {
  roots <- c(Sys.getenv("STRYM_MAMBA_ROOT", ""), Sys.getenv("MAMBA_ROOT_PREFIX", ""),
             path.expand(c("~/mamba", "~/micromamba")), "/root/mamba", "/opt/mamba",
             file.path(dirname(dirname(dirname(R.home("bin")))), ""))
  roots <- roots[nzchar(roots)]
  file.path(sub("/+$", "", roots), "envs", "longread", "bin")
}
.LONG_ENV_BIN <- "/root/mamba/envs/longread/bin"

.long_tool <- function(name) {
  cand <- c(Sys.getenv(paste0("STRYM_", toupper(name)), ""), file.path(c(.LONG_ENV_BIN, .long_env_bins()), name), unname(Sys.which(name)))
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (length(cand)) cand[1] else ""
}
check_long_tools <- function() {
  t <- c("minimap2", "samtools", "featureCounts"); p <- vapply(t, .long_tool, "")
  data.frame(tool = t, found = nzchar(p), path = unname(p), stringsAsFactors = FALSE)
}

# ---- Inputs for RTP-STAR from the Tux results -----------------------------------------------------------------------------
# RTP-STAR (random-forest network) needs two files the researcher would otherwise build by hand:
#   genes file       the differentially expressed genes (first column, with a header)
#   expression file  gene names in column 1, then one column per SAMPLE (replicates in separate columns)
# A Tux run holds both: gene_exp.diff has the group-to-group tests, genes.read_group_tracking the FPKM of every gene in every replicate.
# They are written to <output folder>/rtpstar_inputs/ at the end of a Tux run, and can be built again with other cutoffs in the RTP-STAR tab.
# Genes are named as in the annotation (gene_short_name, like the GENIST files); when fewer than half of the loci have a name the Tux ids are used instead.
# Nothing in the Tux folders is changed.
#
# tux_rtpstar_inputs(): returns list(genes_file, expression_file, readme, n_genes, n_pass, n_significant, capped, n_expression_genes, samples, summary, warnings)
# and stops with a plain message when the results cannot give a useful pair of files (no replicate-level table, fewer than 3 samples, fewer than 2 genes).
tux_rtpstar_inputs <- function(diff_dir, output_folder, samples_used = file.path(dirname(diff_dir), "samples_used.csv"),
                               q_cutoff = 0.05, min_log2fc = 1, min_fpkm = 1, max_genes = 500) {
  f_rg <- file.path(diff_dir, "genes.read_group_tracking"); f_de <- file.path(diff_dir, "gene_exp.diff"); f_tr <- file.path(diff_dir, "genes.fpkm_tracking")
  have <- function(f) file.exists(f) && !is.na(file.size(f)) && file.size(f) > 0
  if (!have(f_de)) stop("gene_exp.diff (the Cuffdiff comparisons) was not found in ", diff_dir, ".", call. = FALSE)
  if (!have(f_rg)) stop("genes.read_group_tracking (the expression of every replicate) was not found in ", diff_dir, ": the expression file for RTP-STAR needs it.", call. = FALSE)
  if (!is.numeric(q_cutoff) || length(q_cutoff) != 1 || is.na(q_cutoff) || q_cutoff <= 0 || q_cutoff > 1) stop("The q-value cutoff must be between 0 and 1.", call. = FALSE)
  if (length(min_log2fc) != 1 || is.na(min_log2fc) || min_log2fc < 0) stop("The minimum fold change must be 0 or more.", call. = FALSE)
  rd <- function(f) utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")

  # expression of every gene in every replicate (a failed locus counts as 0, as everywhere in StRym)
  rg <- rd(f_rg)
  if (!all(c("tracking_id", "condition", "replicate", "FPKM") %in% names(rg))) stop("genes.read_group_tracking does not have the expected columns.", call. = FALSE)
  key <- paste(rg$condition, rg$replicate)
  nm <- rep(NA_character_, nrow(rg))
  if (have(samples_used)) {
    used <- utils::read.csv(samples_used, stringsAsFactors = FALSE)
    if (all(c("group", "k", "sample") %in% names(used))) nm <- as.character(used$sample)[match(key, paste(used$group, used$k - 1))]      # Cuffdiff numbers replicates 0, 1, ...
  }
  nm[is.na(nm)] <- paste0(rg$condition, "_rep", suppressWarnings(as.integer(rg$replicate)) + 1L)[is.na(nm)]
  ord <- unique(data.frame(cond = rg$condition, rep = suppressWarnings(as.integer(rg$replicate)), sample = nm, stringsAsFactors = FALSE))
  ord <- ord[order(factor(ord$cond, levels = unique(rg$condition)), ord$rep), , drop = FALSE]
  samples <- make.unique(ord$sample, sep = "_"); names(samples) <- paste(ord$cond, ord$rep)
  if (length(samples) < 3) stop(sprintf("RTP-STAR needs several samples; the Tux results have %d. With so few samples a network cannot be inferred.", length(samples)), call. = FALSE)
  fp <- suppressWarnings(as.numeric(rg$FPKM)); fp[is.na(fp)] <- 0
  ids <- unique(as.character(rg$tracking_id))
  M <- matrix(0, length(ids), length(samples), dimnames = list(ids, unname(samples)))
  M[cbind(match(as.character(rg$tracking_id), ids), match(key, names(samples)))] <- fp

  # gene names
  de <- rd(f_de)
  need <- c("sample_1", "sample_2", "status", "log2(fold_change)", "q_value")
  if (!all(need %in% names(de))) stop("gene_exp.diff does not have the expected columns.", call. = FALSE)
  lab <- stats::setNames(rep(NA_character_, length(ids)), ids)
  if (have(f_tr)) { tr <- rd(f_tr)
    if (all(c("tracking_id", "gene_short_name") %in% names(tr))) lab[] <- trimws(as.character(tr$gene_short_name))[match(ids, tr$tracking_id)] }
  if (all(is.na(lab)) && all(c("test_id", "gene") %in% names(de))) lab[] <- trimws(as.character(de$gene))[match(ids, de$test_id)]
  lab[lab %in% c("-", "", "NA")] <- NA
  use_ids <- mean(!is.na(lab)) < 0.5                                               # a genome without gene names: the Tux ids are the names
  label <- if (use_ids) stats::setNames(ids, ids) else lab
  keep <- !is.na(label); M <- M[keep, , drop = FALSE]; lbl <- unname(label[keep])
  o <- order(-rowMeans(M)); M <- M[o, , drop = FALSE]; lbl <- lbl[o]                  # one row per gene name: its most expressed locus
  first <- !duplicated(lbl); M <- M[first, , drop = FALSE]; lbl <- lbl[first]
  o <- order(lbl, method = "radix"); M <- M[o, , drop = FALSE]; lbl <- lbl[o]; rownames(M) <- lbl
  if (nrow(M) < 2) stop("Fewer than 2 genes with a name in the Tux results.", call. = FALSE)

  # the differentially expressed genes: significant, and a large change, between any two of the groups
  lf <- suppressWarnings(as.numeric(de[["log2(fold_change)"]])); q <- suppressWarnings(as.numeric(de$q_value))
  dl <- if ("test_id" %in% names(de)) unname(label[match(as.character(de$test_id), ids)]) else rep(NA_character_, nrow(de))
  if ("gene" %in% names(de)) { g2 <- trimws(as.character(de$gene)); g2[g2 %in% c("-", "", "NA")] <- NA; dl[is.na(dl) & !use_ids] <- g2[is.na(dl) & !use_ids] }
  istrue <- function(x) !is.na(x) & x
  hit <- de$status == "OK" & istrue(q <= q_cutoff) & istrue(abs(lf) >= min_log2fc) & !is.na(dl)
  strength <- if (any(hit)) tapply(pmin(abs(lf[hit]), 20), dl[hit], max) else numeric()     # the strongest change of each gene (Inf capped)
  n_sig <- length(strength)
  strength <- strength[names(strength) %in% lbl]
  if (length(strength)) strength <- strength[apply(M[names(strength), , drop = FALSE], 1, function(r) any(r >= min_fpkm))]
  n_pass <- length(strength)
  capped <- !is.na(max_genes) && max_genes > 0 && n_pass > max_genes
  if (capped) strength <- strength[order(-strength, names(strength), method = "radix")][seq_len(max_genes)]
  genes <- if (length(strength)) sort(names(strength), method = "radix") else character()
  if (length(genes) < 2) stop(sprintf("Only %d gene(s) are differentially expressed with these cutoffs (q <= %s, |log2 fold change| >= %s, FPKM >= %s somewhere): too few for a network. Relax the cutoffs.",
                                      length(genes), format(q_cutoff), format(min_log2fc), format(min_fpkm)), call. = FALSE)

  # write: the expression file holds the chosen genes and every other gene that is expressed (RTP-STAR may use more genes than the genes file)
  rows <- unique(c(genes, rownames(M)[apply(M, 1, max) >= min_fpkm])); rows <- sort(rows, method = "radix")
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  genes_file <- file.path(output_folder, "rtpstar_genes_from_tux.csv"); expr_file <- file.path(output_folder, "rtpstar_expression_from_tux.csv")
  readme <- file.path(output_folder, "README_rtpstar_inputs.txt")
  utils::write.csv(data.frame(gene = genes), genes_file, row.names = FALSE)
  utils::write.csv(data.frame(gene = rows, as.data.frame(round(M[rows, , drop = FALSE], 4), check.names = FALSE), check.names = FALSE), expr_file, row.names = FALSE)
  per_group <- table(ord$cond); warn <- character()
  if (any(per_group < 2)) warn <- c(warn, sprintf("%s has a single replicate: RTP-STAR works best with replicates in separate columns.", paste(names(per_group)[per_group < 2], collapse = ", ")))
  if (length(samples) < 6) warn <- c(warn, sprintf("Only %d samples: the random forests have little to learn from; treat the network with care.", length(samples)))
  summary <- sprintf("%s differentially expressed genes (q <= %s, |log2 fold change| >= %s between two groups, FPKM >= %s somewhere)%s; expression of %s genes in %d samples (%s).",
                     format(length(genes), big.mark = ","), format(q_cutoff), format(min_log2fc), format(min_fpkm),
                     if (capped) sprintf("; kept the %s strongest of %s", format(max_genes, big.mark = ","), format(n_pass, big.mark = ",")) else "",
                     format(length(rows), big.mark = ","), length(samples), paste(unname(samples), collapse = ", "))
  writeLines(c("Inputs for RTP-STAR, made from your Tux results", "",
               "rtpstar_genes_from_tux.csv        -> RTP-STAR tab, 'Genes file'",
               "rtpstar_expression_from_tux.csv   -> RTP-STAR tab, 'Gene expression file'", "",
               summary, if (use_ids) "Genes are named by their Tux ids (the annotation has no gene names for most loci)." else "Genes are named as in the annotation (gene_short_name).",
               if (length(warn)) c("", paste("Note:", warn)), "",
               "Optional in RTP-STAR: a transcription factor list (without it every gene may regulate every other and the run takes longer), a gene-name file,",
               "a clustering file and a time course file. For a time course file with the time points in order, use the GENIST tab's 'Start from your Tux results'.",
               "Nothing in the Tux folders was changed. To use other cutoffs, build the files again in the RTP-STAR tab (Start from your Tux results)."), readme)
  list(genes_file = genes_file, expression_file = expr_file, readme = readme, n_genes = length(genes), n_pass = n_pass, n_significant = n_sig, capped = capped,
       n_expression_genes = length(rows), samples = unname(samples), summary = summary, warnings = warn)
}

# At the end of a Tux run: write the RTP-STAR inputs next to the results (<output folder>/rtpstar_inputs/). When the results cannot give them
# (too few samples, too few differentially expressed genes ...) a note says why instead. Never raises an error: the Tux run itself is done.
tux_write_rtpstar_inputs <- function(output_folder, ...) {
  dir_out <- file.path(output_folder, "rtpstar_inputs"); unlink(dir_out, recursive = TRUE)           # never leave the files of an earlier run beside new results
  r <- tryCatch(tux_rtpstar_inputs(file.path(output_folder, "diff_out"), dir_out, ...), error = function(e) e)
  if (inherits(r, "error")) {
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    try(writeLines(c("No inputs for RTP-STAR were made from these Tux results.", "", conditionMessage(r), "",
                     "You can still prepare them yourself, or build them with other cutoffs in the RTP-STAR tab (Start from your Tux results)."),
                   file.path(dir_out, "rtpstar_inputs_NOT_made.txt")), silent = TRUE)
    message("Note: no RTP-STAR inputs were made: ", conditionMessage(r))
    return(invisible(NULL))
  }
  message("RTP-STAR inputs: ", r$summary, " Saved in ", dir_out)
  for (w in r$warnings) message("Note: ", w)
  invisible(r)
}
