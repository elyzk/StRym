# analysis/common.R -- shared code for the analysis steps that work on Tux (Cuffdiff) results.
#
# The baseline scripts (modified Shannon entropy, Shannon entropy, PCA, trajectories) were written for one experiment:
# one folder of files per genotype x time point, "heat" vs "control" in the condition names. Here they are generalized
# to what StRym's Tux step actually produces: ONE Cuffdiff run over all groups, where a group is
#   <condition>            (one time point)    or    <condition>_<time point>   (several time points).
# The experiment is described by the sample sheet (sample, condition, time_point, replicate); the treated-vs-control
# pairs are worked out from the condition names (or set by the researcher in "Analysis settings").
#
# Data model used by every step ("expression data"):
#   $kind     "tux" (FPKM + Cuffdiff comparisons from a Tux run) or "counts" (an expression table, e.g. counts -> CPM, or FPKM as given)
#   $unit     "FPKM" | "CPM" | "FPKM/TPM"
#   $samples  data.frame(sample, condition, time_point, replicate, group)
#   $rep      data.frame(gene, sample, value, bad)   one row per gene x sample ("bad" = Cuffdiff status FAIL/HIDATA)
#   $diff     Cuffdiff's gene_exp.diff. From a Tux run, or (kind "counts") from comparison files the researcher supplied next to the
#             table (design$diff_path); NULL when there are none. The steps that need q-values test for this, not for $kind.
#   $gene_names  named character vector: gene id -> readable name (NA when unknown), or NULL

# (the analysis settings -- strym_default_params(), strym_normalize_params() -- live in pipeline.R)

# ---- reading Tux (Cuffdiff) results ---------------------------------------------------------------------------------
strym_tux_dir <- function(project_dir) file.path(project_dir, "01_quantification", "tux")
strym_have_tux <- function(project_dir) {
  td <- strym_tux_dir(project_dir)
  f <- file.path(td, c(file.path("diff_out", "genes.read_group_tracking"), file.path("diff_out", "gene_exp.diff"), "samples_used.csv"))
  all(file.exists(f)) && all(file.size(f) > 0)
}

strym_load_tux <- function(project_dir) {
  td <- strym_tux_dir(project_dir)
  if (!strym_have_tux(project_dir))
    stop("No Tux (Cuffdiff) results found in this project. Start from raw reads and run “Quantify reads” first.", call. = FALSE)
  used <- utils::read.csv(file.path(td, "samples_used.csv"), stringsAsFactors = FALSE)
  rg <- utils::read.delim(file.path(td, "diff_out", "genes.read_group_tracking"), stringsAsFactors = FALSE)
  de <- utils::read.delim(file.path(td, "diff_out", "gene_exp.diff"), stringsAsFactors = FALSE)
  rg$sample <- used$sample[match(paste(rg$condition, rg$replicate), paste(used$group, used$k - 1))]   # Cuffdiff numbers replicates 0, 1, ...
  if (anyNA(rg$sample)) stop("Could not match Cuffdiff replicates to the samples used.", call. = FALSE)
  sf <- file.path(project_dir, "01_quantification", "sample_sheet.csv")
  sheet <- if (file.exists(sf)) utils::read.csv(sf, stringsAsFactors = FALSE) else NULL
  samples <- data.frame(sample = used$sample, group = used$group, replicate = used$k, stringsAsFactors = FALSE)
  m <- if (!is.null(sheet)) match(samples$sample, sheet$sample) else rep(NA_integer_, nrow(samples))
  samples$condition <- if (!is.null(sheet)) as.character(sheet$condition[m]) else samples$group
  samples$time_point <- if (!is.null(sheet)) as.character(sheet$time_point[m]) else ""
  samples$condition[is.na(samples$condition)] <- samples$group[is.na(samples$condition)]
  samples$time_point[is.na(samples$time_point)] <- ""
  nm <- if ("gene" %in% names(de)) setNames(as.character(de$gene), de$test_id) else NULL
  if (!is.null(nm)) { nm <- nm[!duplicated(names(nm))]; nm[nm %in% c("-", "")] <- NA }
  list(kind = "tux", unit = "FPKM", samples = samples[, c("sample", "condition", "time_point", "replicate", "group")],
       rep = data.frame(gene = as.character(rg$tracking_id), sample = rg$sample, value = suppressWarnings(as.numeric(rg$FPKM)),
                        bad = rg$status %in% c("FAIL", "HIDATA"), stringsAsFactors = FALSE),
       diff = de, gene_names = nm)
}

# ---- reading a counts table (genes x samples) -------------------------------------------------------------------------
# Raw table: genes x samples matrix as given, plus each sample's condition / time point from the sample sheet (without one every
# sample is its own condition). Returns list(m, samples).
#
# Columns that are NOT samples are left out: text columns (gene_name, description ...) and annotation columns such as the Chr / Start / End /
# Strand / Length that featureCounts writes in front of its counts, unless the sample sheet names one as a sample. (Otherwise a text column becomes
# a fake all-zero sample and Start / End / Length become samples with enormous "counts".)
.ANNOT_NAMES <- paste0("^(chr|chrom|chromosome|contig|start|end|stop|strand|length|width|size|gene[ ._-]?(name|symbol|id|type|biotype|description|desc|short[ ._-]?name)|",
                       "symbol|name|description|desc|biotype|type|source|feature|annotation|locus|accession|ensembl|entrez|product|note|notes|",
                       "(gene|transcript|effective|exon|cds|eff|tx)[ ._-]?len(gth)?|num[ ._-]?reads?|tpm|fpkm)$")   # lengths and the other unit's columns of Salmon / featureCounts tables
strym_count_columns <- function(tab, listed = character()) {
  cols <- seq_along(tab)[-1]; nm <- names(tab)[cols]
  texty <- vapply(tab[cols], function(v) { s <- trimws(as.character(v)); s <- s[!is.na(s) & nzchar(s) & !s %in% c("NA", "N/A", "n/a", "NaN", "-", ".", "null", "NULL")]
                                           !(length(s) > 0 && mean(!is.na(.to_num(s))) >= 0.5) }, logical(1))      # numbers as spreadsheets write them (1,234 / 1,5) count as numbers
  annot <- grepl(.ANNOT_NAMES, trimws(nm), ignore.case = TRUE) & !(tolower(trimws(nm)) %in% tolower(trimws(listed)))
  list(keep = cols[!texty & !annot], dropped = nm[texty | annot])
}

.counts_raw <- function(data_path, sample_sheet = "") {
  tab <- read_table_any(data_path)
  ss <- NULL; listed <- character()
  if (nzchar(sample_sheet) && file.exists(sample_sheet)) {
    ss <- strym_read_sample_sheet(sample_sheet)
    if (inherits(ss, "error")) stop("Sample sheet: ", conditionMessage(ss), call. = FALSE)
    sc0 <- .find_col(ss, "^sample")
    if (is.na(sc0)) stop("The sample sheet needs a column named “sample” with the column names of the table (it has: ", paste(utils::head(names(ss), 6), collapse = ", "), ").", call. = FALSE)
    listed <- trimws(as.character(ss[[sc0]]))
  }
  cc <- strym_count_columns(tab, listed)
  if (length(cc$dropped)) message("Note: left out ", length(cc$dropped), " column(s) of the counts table that are not samples: ", paste(head(cc$dropped, 8), collapse = ", "), if (length(cc$dropped) > 8) ", ..." else "", ".")
  if (!length(cc$keep)) stop("The counts table has no numeric sample columns (every column after the gene names is text or annotation).", call. = FALSE)
  # A readable-name column (gene_short_name, gene_name, symbol ...) is not data, but it is kept to label the genes in the results.
  namecol <- NULL
  cand <- cc$dropped[grepl("^(gene[ ._-]?(name|symbol|short[ ._-]?name)|symbol|name)$", trimws(cc$dropped), ignore.case = TRUE)]
  if (length(cand)) namecol <- trimws(as.character(tab[[cand[1]]]))
  tab <- tab[c(1, cc$keep)]
  genes <- trimws(as.character(tab[[1]]))
  blank <- is.na(genes) | !nzchar(genes) | genes == "NA"; dup <- !blank & duplicated(genes); keep <- !blank & !dup
  if (any(!keep)) message(sprintf("Note: %d row(s) of the table were left out: %d without a gene id, %d with a gene id already listed (the first one is kept).", sum(!keep), sum(blank), sum(dup)))
  dropped <- list(empty = character(), unlisted = character(), missing = character())
  gnames <- NULL
  if (!is.null(namecol) && length(namecol) == length(genes)) {
    nmk <- namecol[keep]; gk <- genes[keep]; ok <- !is.na(nmk) & nzchar(nmk) & !nmk %in% c("-", "NA") & nmk != gk
    if (any(ok)) gnames <- stats::setNames(nmk[ok], gk[ok])
  }
  m <- sapply(tab[-1], .to_num, simplify = TRUE)
  if (!is.matrix(m)) m <- matrix(m, nrow = nrow(tab), dimnames = list(NULL, names(tab)[-1]))
  rownames(m) <- genes; m <- m[keep, , drop = FALSE]
  dropped$rows_blank <- sum(blank); dropped$rows_dup <- sum(dup)
  if (any(is.infinite(m))) { w <- which(is.infinite(m), arr.ind = TRUE)                              # Inf has no meaning in a table of counts or FPKM: like an empty cell it is read as 0, and said so
    message(sprintf("Note: %d cell(s) hold Inf / -Inf and were read as 0 (first: gene %s, sample %s).", nrow(w), rownames(m)[w[1, 1]], colnames(m)[w[1, 2]]))
    dropped$inf_cells <- nrow(w); m[is.infinite(m)] <- 0 }
  na <- is.na(m); dropped$na_cells <- sum(na)
  if (any(na)) { message(sprintf("Note: %s empty or non-numeric cell(s) (in %d sample%s) were read as 0.", format(sum(na), big.mark = ","), sum(colSums(na) > 0), if (sum(colSums(na) > 0) == 1) "" else "s")); m[na] <- 0 }
  empty <- colSums(m) == 0                                              # a library with nothing in it carries no information (and breaks the models)
  if (any(empty)) { message("Note: sample(s) with no counts at all were left out: ", paste(colnames(m)[empty], collapse = ", "), "."); dropped$empty <- colnames(m)[empty]; m <- m[, !empty, drop = FALSE]
    if (!ncol(m)) stop("Every sample of the table is empty (all zeros).", call. = FALSE) }
  smp <- colnames(m); cond <- smp; tp <- rep("", length(smp)); rp <- seq_along(smp)
  if (!is.null(ss)) {
    sc <- .find_col(ss, "^sample"); cc <- .find_col(ss, "^condition"); tc <- .find_col(ss, "^time"); rc <- .find_col(ss, "^replicate")
    nm_sheet <- trimws(as.character(ss[[sc]])); i <- match(smp, nm_sheet)
    if (!anyDuplicated(tolower(smp)) && !anyDuplicated(tolower(nm_sheet))) { j <- match(tolower(smp), tolower(nm_sheet)); i[is.na(i)] <- j[is.na(i)] }   # same name written in another case (only when that is unambiguous)
    if (all(is.na(i))) stop(sprintf("None of the table's sample columns (%s ...) is listed in the sample sheet (%s ...). The names must match exactly.",
                                    paste(utils::head(smp, 3), collapse = ", "), paste(utils::head(nm_sheet, 3), collapse = ", ")), call. = FALSE)
    if (anyNA(i)) {                                                                      # table columns the sheet does not describe are not part of the experiment
      message(sprintf("Note: %d column(s) of the table are not in the sample sheet and were left out: %s.", sum(is.na(i)), paste(utils::head(smp[is.na(i)], 6), collapse = ", ")))
      dropped$unlisted <- smp[is.na(i)]
      m <- m[, !is.na(i), drop = FALSE]; i <- i[!is.na(i)]; smp <- colnames(m); tp <- rep("", length(smp)); rp <- seq_along(smp)
    }
    miss <- setdiff(tolower(nm_sheet[nzchar(nm_sheet)]), tolower(smp)); dropped$missing <- miss
    if (length(miss)) message(sprintf("Note: %d sample(s) of the sample sheet are not in the table and were ignored: %s.", length(miss), paste(utils::head(miss, 6), collapse = ", ")))
    cond <- if (!is.na(cc)) trimws(as.character(ss[[cc]]))[i] else smp
    if (any(is.na(cond) | !nzchar(cond))) stop("The sample sheet has an empty condition for: ", paste(utils::head(smp[is.na(cond) | !nzchar(cond)], 6), collapse = ", "), ".", call. = FALSE)
    if (!is.na(tc)) tp <- trimws(as.character(ss[[tc]]))[i]; tp[is.na(tp)] <- ""
    if (!is.na(rc)) { r <- suppressWarnings(as.integer(ss[[rc]]))[i]; rp <- ifelse(is.na(r), seq_along(smp), r) }
  }
  samples <- data.frame(sample = smp, condition = cond, time_point = tp, replicate = rp, group = ifelse(nzchar(tp), paste(cond, tp, sep = "_"), cond),
                        stringsAsFactors = FALSE)
  list(m = m, samples = samples, gene_names = gnames, dropped = dropped)
}

# Text -> numbers, tolerant of what spreadsheets export: blanks, NA, "-", thousands separators (1,234.5) and decimal commas (1,5).
.to_num <- function(v) {
  if (is.numeric(v)) return(as.numeric(v))
  s <- trimws(as.character(v)); s[s %in% c("", "NA", "N/A", "n/a", "NaN", "-", ".", "null", "NULL")] <- NA
  x <- suppressWarnings(as.numeric(s)); bad <- is.na(x) & !is.na(s)
  if (any(bad)) {
    s2 <- s[bad]
    dec_style <- any(grepl("^-?[0-9]+,([0-9]{1,2}|[0-9]{4,})$", s2))              # a comma followed by 1-2 or 4+ digits can only be a decimal comma: this column uses them
    th <- !dec_style & grepl("^-?[1-9][0-9]{0,2}(,[0-9]{3})+(\\.[0-9]+)?$", s2); dc <- !th & grepl("^-?[0-9]+,[0-9]+$", s2)       # "0,125" is never a thousands separator
    s2[th] <- gsub(",", "", s2[th], fixed = TRUE); s2[dc] <- sub(",", ".", s2[dc], fixed = TRUE)
    x[bad] <- suppressWarnings(as.numeric(s2))
  }
  x
}

# Expression table -> expression data. `unit` says what the numbers are:
#   "counts"      raw counts: library-size normalized to counts per million so samples are comparable
#   "normalized"  already normalized (FPKM / TPM / CPM ...): used as they are (re-scaling them would distort them)
#   "auto"        whole numbers = counts, anything fractional = already normalized
# `diff_path` (optional) = results of a differential-expression run done elsewhere (Cuffdiff gene_exp.diff files: one file, several,
# a folder or a .zip). They give the q-values the time-course steps (Shannon entropy, modified Shannon entropy) need, so a
# table of FPKM values plus those files works exactly like a Tux run done inside StRym.
strym_load_counts <- function(data_path, sample_sheet = "", unit = "auto", diff_path = character(), log = message) {
  raw <- .counts_raw(data_path, sample_sheet); m <- raw$m
  unit <- tolower(trimws(as.character(unit)[1])); if (is.na(unit) || !unit %in% c("auto", "counts", "normalized")) unit <- "auto"
  if (any(m < 0)) { w <- which(m < 0, arr.ind = TRUE)[1, ]
    stop(sprintf("The table has negative values (first: gene %s, sample %s). StRym needs counts or non-negative expression values (FPKM / TPM / CPM), not log2 or z-scored values.",
                 rownames(m)[w[1]], colnames(m)[w[2]]), call. = FALSE) }
  whole <- isTRUE(all(m == floor(m)))
  # A table of counts with a few stray decimals (a "1,5" typed by hand, a rounding artefact) is still a table of counts: they are rounded and said so.
  # A table where decimals are the rule (FPKM / TPM / CPM) stays "already normalized".
  if (!whole && identical(unit, "auto") && mean(m != floor(m)) < 0.001) {
    w <- which(m != floor(m), arr.ind = TRUE)
    log(sprintf("Note: %d cell(s) hold decimals in a table of whole numbers (first: gene %s, sample %s). The table is taken as raw counts and they are rounded.", nrow(w), rownames(m)[w[1, 1]], colnames(m)[w[1, 2]]))
    m <- round(m); whole <- TRUE; raw$dropped$decimal_cells <- nrow(w)
  }
  as_counts <- identical(unit, "counts") || (identical(unit, "auto") && whole)
  if (identical(unit, "auto") && !whole) {                         # decimals: taken as already normalized, but expected counts (RSEM) look the same
    cs <- colSums(m); log(sprintf("The values have decimals, so they are used as they are (already normalized). If this is a table of expected counts (RSEM, Salmon, kallisto) choose “Raw counts” in the design%s.",
                                  if (max(cs) > 3 * min(cs)) sprintf(" (the column totals differ %.0f-fold, which normalized values usually do not)", max(cs) / min(cs)) else ""))
  }
  if (as_counts) {
    lib <- colSums(m)
    val <- m; for (j in seq_len(ncol(m))) val[, j] <- m[, j] * (1e6 / lib[j]); ulab <- "CPM"       # counts per million, one column at a time (no second full copy)
  } else { val <- m; ulab <- "FPKM/TPM" }
  # the expression values are kept as ONE genes x samples matrix (a long gene x sample table costs ~20x more memory)
  out <- list(kind = "counts", unit = ulab, whole_numbers = whole, samples = raw$samples, counts = if (as_counts) m else NULL,    # counts: the raw table, kept for DESeq2-derived q-values
              wide = val, rep = NULL, diff = NULL, gene_names = raw$gene_names, dropped = raw$dropped)
  diff_path <- as.character(unlist(diff_path)); diff_path <- diff_path[!is.na(diff_path) & nzchar(trimws(diff_path))]
  if (length(diff_path)) {
    ld <- strym_load_diff(diff_path, raw$samples, rownames(val), log = log, values = if (as_counts) NULL else val)
    out$diff <- ld$diff
    nm <- ld$gene_names
    if (!is.null(raw$gene_names)) nm <- c(raw$gene_names, nm[setdiff(names(nm), names(raw$gene_names))])   # names in your table win
    out$gene_names <- nm
  }
  out
}

# ---- differential-expression results made elsewhere (Cuffdiff gene_exp.diff) ---------------------------------------------
# Whatever the experiment, the input is "some Cuffdiff comparison files": one file with every comparison, or one per genotype /
# time point (as they come out of separate runs), loose, in folders, or zipped. What matters is only which two groups each file
# compares (its sample_1 / sample_2 columns) and which genes (test_id / gene_id / gene). Both are matched to the sample sheet,
# whatever the naming (T1_Pana_Control, Pana_Control_T1 or just "Control" in a folder called Pana_T1 all resolve to the same group).
.DIFF_SKIP <- "^(isoform|tss_group|cds|splicing|promoters)([_.]|$)|^\\._|\\.(png|jpe?g|svg|pdf|html?|rds|bam|bai|log|md)(\\.gz)?$"      # other Cuffdiff comparisons (isoform, cds.diff, promoters.diff ...), pictures, notes
.diff_tok <- function(x) { t <- tolower(unlist(strsplit(as.character(x), "[^A-Za-z0-9]+"))); t[nzchar(t)] }

# Expand files / folders / zip archives into data.frame(file, context): `context` = words from the folder and file names that can
# tell which genotype / time point a file belongs to when its own sample labels do not.
#
# However they are packed, they are found: a file, a folder (any depth), a .zip, a folder of folders holding zips, a zip holding folders
# holding zips ... (unpacked once and remembered, so the page's checks and every step do not unpack them again).
.DIFF_UNZIPPED <- new.env()
.unzip_cached <- function(z) {
  info <- file.info(z); key <- paste(normalizePath(z, mustWork = FALSE), info$size, as.numeric(info$mtime))
  d <- .DIFF_UNZIPPED[[key]]
  if (is.null(d) || !dir.exists(d)) {
    d <- tempfile("strym_zip_"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    tryCatch(utils::unzip(z, exdir = d), error = function(e) stop("Could not open ", basename(z), " (", conditionMessage(e), ").", call. = FALSE),
             warning = function(w) stop("Could not open ", basename(z), " (", conditionMessage(w), ").", call. = FALSE))
    .DIFF_UNZIPPED[[key]] <- d
  }
  d
}
.looks_like_diff <- function(f) {                  # a comparison table by its header, whatever the file is called
  h <- tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) "", warning = function(w) "")
  length(h) == 1 && grepl("sample_1", h, fixed = TRUE) && grepl("sample_2", h, fixed = TRUE)
}
# Walk `root`: comparison files, and zips (opened and walked in turn). `ctx` = words from the names above this point.
.walk_diff <- function(root, ctx, depth = 0, loose = FALSE) {
  root <- sub("[/\\\\]+$", "", root)
  fs <- list.files(root, recursive = TRUE, full.names = TRUE)
  fs <- fs[!grepl("__MACOSX", fs, fixed = TRUE) & !startsWith(basename(fs), "._")]
  res <- list()
  for (f in fs) {
    rel <- substring(f, nchar(root) + 2L)
    if (grepl("\\.zip$", f, ignore.case = TRUE)) {
      if (depth < 5) res[[length(res) + 1]] <- .walk_diff(.unzip_cached(f), paste(ctx, sub("\\.zip$", "", rel, ignore.case = TRUE)), depth + 1, loose)
    } else if (!grepl(.DIFF_SKIP, basename(f), ignore.case = TRUE) &&
               (if (loose) grepl("\\.(diff|tsv|csv|txt)(\\.gz)?$", f, ignore.case = TRUE) && .looks_like_diff(f) else grepl("diff", basename(f), ignore.case = TRUE)))
      res[[length(res) + 1]] <- data.frame(file = f, context = paste(ctx, rel), stringsAsFactors = FALSE)
  }
  if (length(res)) do.call(rbind, res) else data.frame(file = character(), context = character(), stringsAsFactors = FALSE)
}

strym_diff_files <- function(paths) {
  paths <- trimws(as.character(unlist(paths))); paths <- sub("[/\\\\]+$", "", paths[nzchar(paths)])
  out <- list()
  for (p in paths) {
    if (!file.exists(p)) stop("Differential-expression results not found: ", p, call. = FALSE)
    one <- function(loose) {
      if (dir.exists(p)) .walk_diff(p, basename(p), 0, loose)
      else if (grepl("\\.zip$", p, ignore.case = TRUE)) .walk_diff(.unzip_cached(p), sub("\\.zip$", "", basename(p), ignore.case = TRUE), 1, loose)
      else data.frame(file = p, context = paste(basename(dirname(p)), basename(p)), stringsAsFactors = FALSE)   # a file named by hand is taken as it is
    }
    r <- one(FALSE); if (!nrow(r)) r <- one(TRUE)                  # none named "...diff...": look at the headers instead
    out[[length(out) + 1]] <- r
  }
  if (length(out)) do.call(rbind, out) else data.frame(file = character(), context = character(), stringsAsFactors = FALSE)
}

# How many comparison files a set of paths holds (used by the page's checks).
strym_count_diff_files <- function(paths) nrow(strym_diff_files(paths))

# One comparison file -> data.frame(id_test, id_gene, id_name, sample_1, sample_2, status, q_value, log2fc). Column names are matched
# loosely (Cuffdiff's own: test_id, gene_id, gene, sample_1, sample_2, status, q_value, log2(fold_change); q_value may also be padj / fdr).
.read_diff_file <- function(f) {
  first <- readLines(f, n = 1, warn = FALSE)
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(",", first, fixed = TRUE)) "," else ""
  d <- utils::read.table(f, header = TRUE, sep = sep, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE,
                         strip.white = TRUE, fileEncoding = "UTF-8-BOM")
  nm <- tolower(trimws(names(d)))
  pick <- function(...) { i <- match(c(...), nm); i <- i[!is.na(i)]; if (length(i)) i[1] else NA_integer_ }
  s1 <- pick("sample_1"); s2 <- pick("sample_2"); q <- pick("q_value", "qvalue", "q.value", "padj", "fdr")
  if (anyNA(c(s1, s2, q)))
    stop(basename(f), " is not a Cuffdiff comparison file: it needs the columns sample_1, sample_2 and q_value (it has: ", paste(utils::head(names(d), 8), collapse = ", "), ").", call. = FALSE)
  col <- function(i) if (is.na(i)) rep(NA_character_, nrow(d)) else trimws(as.character(d[[i]]))
  st <- pick("status"); lf <- pick("log2(fold_change)", "log2fc", "log2foldchange"); v1 <- pick("value_1"); v2 <- pick("value_2")
  num <- function(i) if (is.na(i)) rep(NA_real_, nrow(d)) else suppressWarnings(as.numeric(d[[i]]))
  data.frame(value_1 = num(v1), value_2 = num(v2),id_test = col(pick("test_id")), id_gene = col(pick("gene_id")), id_name = col(pick("gene", "gene_name", "gene_short_name")),
             sample_1 = col(s1), sample_2 = col(s2), status = if (is.na(st)) rep("OK", nrow(d)) else col(st),
             q_value = suppressWarnings(as.numeric(d[[q]])), log2fc = if (is.na(lf)) rep(NA_real_, nrow(d)) else suppressWarnings(as.numeric(d[[lf]])),
             stringsAsFactors = FALSE)
}

# The sample-sheet group (condition_timepoint, or the condition) each Cuffdiff sample label stands for; NA = no match.
# A label matches a group when it has the same words (in any order) as the group's condition + time point; a label that is only part
# of that (e.g. "Control") is completed with the words of the file / folder name it came from (`context`).
strym_resolve_diff_groups <- function(labels, samples, context = "") {
  g <- samples[!duplicated(samples$group), c("group", "condition", "time_point"), drop = FALSE]
  gt <- lapply(seq_len(nrow(g)), function(i) unique(c(.diff_tok(g$condition[i]), .diff_tok(g$time_point[i]))))
  cx <- setdiff(.diff_tok(context), c("gene", "exp", "diff", "out", "zip", "csv", "tsv", "txt"))
  vapply(labels, function(lab) {
    i <- which(tolower(g$group) == tolower(lab)); if (length(i) == 1) return(g$group[i])
    L <- unique(.diff_tok(lab)); if (!length(L)) return(NA_character_)
    i <- which(vapply(gt, function(G) setequal(G, L), TRUE)); if (length(i) == 1) return(g$group[i])
    i <- which(vapply(gt, function(G) all(L %in% G) && all(G %in% c(L, cx)), TRUE))
    if (length(i) > 1) { n <- lengths(gt[i]); i <- i[n == max(n)] }
    if (length(i) == 1) g$group[i] else NA_character_
  }, "", USE.NAMES = TRUE)
}

# Last resort for a label whose NAME matches nothing (the run called a genotype "Dujia", the table calls it "PI567357"): Cuffdiff writes the
# mean expression of each group (value_1 / value_2), so a label is the sample-sheet group whose mean expression in the table is the same
# numbers. `values` = genes x samples matrix of the table (only meaningful when it holds FPKM, as Cuffdiff does). Returns the group or NA.
.resolve_label_by_values <- function(d, lab, values, samples) {
  if (is.null(values) || !nrow(d)) return(NA_character_)
  ids <- c("id_test", "id_gene", "id_name")
  ov <- vapply(ids, function(k) { x <- utils::head(d[[k]][!is.na(d[[k]])], 2000); if (length(x)) mean(x %in% rownames(values)) else 0 }, 0)
  if (max(ov) < 0.5) return(NA_character_)
  idc <- ids[which.max(ov)]
  v <- c(d$value_1[d$sample_1 == lab], d$value_2[d$sample_2 == lab]); id <- c(d[[idc]][d$sample_1 == lab], d[[idc]][d$sample_2 == lab])
  k <- which(!is.na(v) & v > 0 & !is.na(match(id, rownames(values)))); if (length(k) < 20) return(NA_character_)
  k <- k[round(seq(1, length(k), length.out = min(length(k), 3000)))]
  v <- v[k]; r <- match(id[k], rownames(values))
  groups <- unique(samples$group); frac <- vapply(groups, function(g) {
    cols <- intersect(samples$sample[samples$group == g], colnames(values)); if (!length(cols)) return(0)
    gm <- rowMeans(values[r, cols, drop = FALSE]); mean(abs(gm - v) <= 1e-4 * pmax(abs(gm), abs(v)) + 1e-9)
  }, 0)
  o <- order(-frac); if (frac[o[1]] >= 0.9 && (length(o) < 2 || frac[o[2]] < 0.9)) groups[o[1]] else NA_character_
}

# Read every comparison, translate their labels to the sample sheet's groups and their gene ids to the table's. Returns
# list(diff = data.frame(test_id, sample_1, sample_2, status, q_value, log2fc), gene_names = named vector or NULL).
.diff_loaded <- new.env()           # the parsed comparisons, remembered for the rest of the run (each entropy step would otherwise read the same files again)
strym_load_diff <- function(paths, samples, genes, log = message, values = NULL) {
  files <- strym_diff_files(paths)
  fi <- file.info(files$file); key <- paste(c(files$file, fi$size, as.numeric(fi$mtime), samples$sample, samples$group, length(genes), genes[1], genes[length(genes)],
                                              if (!is.null(values)) format(sum(values), digits = 12)), collapse = "|")
  if (!is.null(hit <- .diff_loaded[[key]])) return(hit)
  res <- .strym_load_diff_uncached(files, samples, genes, log, values); assign(key, res, envir = .diff_loaded); res
}
.strym_load_diff_uncached <- function(files, samples, genes, log, values) {
  if (!nrow(files)) stop("No Cuffdiff comparison files (gene_exp.diff) were found in the differential-expression results you gave.", call. = FALSE)
  parts <- vector("list", nrow(files)); unresolved <- character(); used <- 0L; by_value <- character()
  for (i in seq_len(nrow(files))) {
    d <- tryCatch(.read_diff_file(files$file[i]), error = function(e) e)
    if (inherits(d, "error")) {                                                    # a file that only looks like a comparison file (its name has "diff" in it) is skipped, not fatal
      if (grepl("gene_exp", basename(files$file[i]), fixed = TRUE) || nrow(files) == 1) stop(conditionMessage(d), call. = FALSE)
      log(sprintf("Skipped %s: %s", basename(files$file[i]), conditionMessage(d))); next
    }
    labs <- unique(c(d$sample_1, d$sample_2))
    map <- strym_resolve_diff_groups(labs, samples, files$context[i])
    for (lab in labs[is.na(map[labs])]) {                                          # the name told nothing: try the numbers
      g <- .resolve_label_by_values(d, lab, values, samples)
      if (!is.na(g) && !g %in% map) { map[lab] <- g; by_value <- c(by_value, paste0(lab, " = ", g)) }
    }
    d$sample_1 <- unname(map[d$sample_1]); d$sample_2 <- unname(map[d$sample_2])
    unresolved <- c(unresolved, labs[is.na(map[labs])])
    ok <- !is.na(d$sample_1) & !is.na(d$sample_2) & d$sample_1 != d$sample_2
    if (any(ok)) { used <- used + 1L; parts[[i]] <- d[ok, , drop = FALSE] }
  }
  de <- do.call(rbind, parts[!vapply(parts, is.null, TRUE)])
  if (is.null(de) || !nrow(de))
    stop("None of the comparisons in your differential-expression results match a condition / time point of the sample sheet. The results compare: ",
         paste(utils::head(unique(unresolved), 6), collapse = ", "), "; the sample sheet has: ", paste(utils::head(unique(samples$group), 8), collapse = ", "), ".", call. = FALSE)
  de <- de[!duplicated(de[c("id_test", "id_gene", "id_name", "sample_1", "sample_2")]), , drop = FALSE]   # the same results given twice (a zip and its unpacked copy)
  ids <- c("id_test", "id_gene", "id_name")
  ov <- vapply(ids, function(k) { x <- utils::head(de[[k]][!is.na(de[[k]])], 20000); if (length(x)) mean(x %in% genes) else 0 }, 0)
  if (max(ov) < 0.05)
    stop("The gene ids in your differential-expression results do not match the gene ids of the expression table (results: ",
         paste(utils::head(unique(de$id_test), 3), collapse = ", "), "; table: ", paste(utils::head(genes, 3), collapse = ", "), ").", call. = FALSE)
  best <- ids[which.max(ov)]; de$test_id <- de[[best]]
  nm <- NULL
  if (best != "id_name" && any(!is.na(de$id_name))) {
    k <- !is.na(de$id_name) & nzchar(de$id_name) & !de$id_name %in% c("-", "NA") & de$id_name != de$test_id & !duplicated(de$test_id)
    if (any(k)) nm <- stats::setNames(de$id_name[k], de$test_id[k])
  }
  de$status[is.na(de$status) | !nzchar(de$status)] <- "OK"
  log(sprintf("Differential-expression results: %d comparison file(s) used, %s comparisons between %d groups (genes matched by %s).",
              used, format(nrow(de), big.mark = ","), length(unique(c(de$sample_1, de$sample_2))), sub("^id_", "", best)))
  if (length(by_value)) log(sprintf("Matched %d group label(s) of the results to your samples by their expression values (their names differ): %s.",
                                    length(by_value), paste(utils::head(by_value, 4), collapse = "; ")))
  unresolved <- setdiff(unresolved, sub(" = .*$", "", by_value))
  if (length(unresolved)) log(sprintf("Left out comparisons of groups that are not in the sample sheet: %s.", paste(utils::head(unique(unresolved), 6), collapse = ", ")))
  list(diff = de[, c("test_id", "sample_1", "sample_2", "status", "q_value", "log2fc"), drop = FALSE], gene_names = nm)
}

# What a step receives: Tux results if the project has them (raw reads), otherwise the counts table.
strym_load_expression <- function(ctx) {
  d <- ctx$design
  if (identical(d$data_kind, "fastq") && identical(d$read_platform, "long")) {                      # long reads: the counts table the quantify step wrote
    cp <- strym_counts_path(d); sh <- file.path(ctx$project_dir, "01_quantification", "sample_sheet.csv")
    if (!file.exists(cp)) stop("No counts found in this project. Run \u201cQuantify reads\u201d first.", call. = FALSE)
    return(strym_load_counts(cp, if (file.exists(sh)) sh else ctx$sample_sheet, unit = "counts", log = if (is.function(ctx$log)) ctx$log else message))
  }
  if (identical(d$data_kind, "fastq") && strym_have_tux(ctx$project_dir)) return(strym_load_tux(ctx$project_dir))
  if (identical(d$data_kind, "fastq")) stop("No Tux (Cuffdiff) results found in this project. Run “Quantify reads” first.", call. = FALSE)
  if (!identical(d$data_kind, "counts")) stop("This step reads a counts table or Tux results.", call. = FALSE)
  strym_load_counts(ctx$data_path, ctx$sample_sheet, unit = d$expr_unit, diff_path = d$diff_path, log = if (is.function(ctx$log)) ctx$log else message)
}

# Samples left out by the researcher (e.g. outliers seen in the PCA).
strym_drop_samples <- function(data, params) {
  drop <- strym_split_list(params$drop_reps)
  unknown <- setdiff(drop, data$samples$sample)
  if (length(unknown)) warning("Not a sample name, ignored: ", paste(unknown, collapse = ", "), call. = FALSE)
  data$samples <- data$samples[!data$samples$sample %in% drop, , drop = FALSE]
  if (!is.null(data$wide)) data$wide <- data$wide[, data$samples$sample, drop = FALSE]
  if (!is.null(data$rep)) data$rep <- data$rep[data$rep$sample %in% data$samples$sample, , drop = FALSE]
  data
}

# Gene ids of the expression data, and its values in long form (gene, sample, value) for the given genes only. The data is held as a
# genes x samples matrix ($wide) when it comes from a table (a long table of every gene in every sample costs ~20x the memory), or as
# the long table ($rep) when it comes from a Tux run.
strym_gene_ids <- function(data) if (!is.null(data$wide)) rownames(data$wide) else unique(data$rep$gene)
strym_long <- function(data, genes) {
  if (is.null(data$wide)) return(data$rep[data$rep$gene %in% genes, c("gene", "sample", "value")])
  g <- intersect(genes, rownames(data$wide)); m <- data$wide[g, , drop = FALSE]
  data.frame(gene = rep(rownames(m), times = ncol(m)), sample = rep(colnames(m), each = nrow(m)), value = as.vector(m), stringsAsFactors = FALSE)
}

# genes x samples matrix of the expression values. na_bad = TRUE sets FAIL / HIDATA values to NA (as the entropy scripts do).
strym_wide <- function(data, na_bad = FALSE) {
  if (!is.null(data$wide)) { m <- data$wide[, data$samples$sample, drop = FALSE]; return(m[order(rownames(m), method = "radix"), , drop = FALSE]) }
  r <- data$rep; genes <- sort(unique(r$gene)); smp <- data$samples$sample
  v <- r$value; if (na_bad) v[r$bad] <- NA
  i <- match(r$gene, genes); j <- match(r$sample, smp); ok <- !is.na(i) & !is.na(j)
  m <- matrix(NA_real_, length(genes), length(smp), dimnames = list(genes, smp))
  m[cbind(i[ok], j[ok])] <- v[ok]
  m
}

# Time points in the order of the experiment: the order typed in the design, else by the number in the label (0h, 1h, 24h), else as found.
strym_order_timepoints <- function(tp, names = character()) {
  tp <- unique(as.character(tp)); tp <- tp[nzchar(tp)]
  if (length(names) && all(tp %in% names)) return(names[names %in% tp])
  v <- .tp_value(tp)
  if (!anyNA(v) && !anyDuplicated(v)) return(tp[order(v)])
  tp
}

# A time label as a number (in hours when it has a unit): "30min" 0.5, "1h" 1, "2d" 48, "T3" 3, "Day10" 10, "-1h" -1. NA when it has no number.
# A "-" is a minus sign only in front of a bare number ("-1h"); in "T-2" it separates the letter from the number.
.tp_value <- function(x) {
  x <- trimws(as.character(x)); out <- rep(NA_real_, length(x))
  mult <- function(u) {
    u <- tolower(sub("(pi|ai|pt|pf|ps|pa|pg)$", "", tolower(u)))                       # 24hpi, 3dpi, 6hai ...
    if (!nzchar(u)) 1 else if (u %in% c("s", "sec", "secs", "second", "seconds")) 1 / 3600 else if (u %in% c("m", "min", "mins", "minute", "minutes")) 1 / 60
    else if (u %in% c("h", "hr", "hrs", "hour", "hours")) 1 else if (u %in% c("d", "day", "days")) 24 else if (u %in% c("w", "wk", "wks", "week", "weeks")) 168 else NA_real_
  }
  for (i in seq_along(x)) {
    r <- regmatches(x[i], regexec("^(-?)([0-9]+(?:\\.[0-9]+)?)\\s*([A-Za-z]*)$", x[i], perl = TRUE))[[1]]
    if (!length(r)) { r2 <- regmatches(x[i], regexec("^[A-Za-z_ ]+[-_ ]?([0-9]+(?:\\.[0-9]+)?)\\s*([A-Za-z]*)$", x[i], perl = TRUE))[[1]]
                      if (length(r2)) r <- c(r2[1], "", r2[2], r2[3]) }
    if (length(r)) { mu <- mult(r[4]); if (!is.na(mu)) out[i] <- (if (r[2] == "-") -1 else 1) * as.numeric(r[3]) * mu }
  }
  out
}

# ---- treated vs control pairs ------------------------------------------------------------------------------------------
# Returns data.frame(treated, control). Order: as given by the researcher, else the order of the conditions.
strym_pairs <- function(samples, params) {
  conds <- unique(samples$condition)
  pr <- strsplit(strym_split_list(params$pairs), "\\s*[:=]\\s*")
  if (length(pr)) {
    if (any(lengths(pr) != 2)) stop("Pairs must look like “Heat:Control, Heat2:Control2” (treated:control).", call. = FALSE)
    out <- data.frame(treated = vapply(pr, `[`, "", 1), control = vapply(pr, `[`, "", 2), stringsAsFactors = FALSE)
    bad <- setdiff(c(out$treated, out$control), conds)
    if (length(bad)) stop("Unknown condition in the pairs: ", paste(bad, collapse = ", "), ". The conditions are: ", paste(conds, collapse = ", "), ".", call. = FALSE)
    return(out)
  }
  ctrl <- strym_split_list(params$control)
  if (length(ctrl)) {
    bad <- setdiff(ctrl, conds)
    if (length(bad)) stop("The control condition “", paste(bad, collapse = ", "), "” is not one of the conditions: ", paste(conds, collapse = ", "), ".", call. = FALSE)
  } else ctrl <- conds[tryCatch(grepl(params$control_pattern, conds, ignore.case = TRUE),
                                error = function(e) stop("The control pattern in Analysis settings is not a valid pattern (", conditionMessage(e), ").", call. = FALSE))]
  if (!length(ctrl)) stop("StRym could not tell which condition is the control (looked for names like control, mock, untreated, WT). ",
                          "Set it in Analysis settings. The conditions are: ", paste(conds, collapse = ", "), ".", call. = FALSE)
  treated <- setdiff(conds, ctrl)
  if (!length(treated)) stop("Every condition looks like a control: nothing to compare.", call. = FALSE)
  if (length(ctrl) == 1) return(data.frame(treated = treated, control = ctrl, stringsAsFactors = FALSE))
  tok <- function(x) tolower(unlist(strsplit(x, "[^A-Za-z0-9]+")))               # several controls: pair by the name they share (Pana_Heat <-> Pana_Control)
  pick <- vapply(treated, function(t) {
    ov <- vapply(ctrl, function(c) length(intersect(tok(t), tok(c))), 1L)
    if (max(ov) == 0 || sum(ov == max(ov)) > 1) NA_character_ else ctrl[which.max(ov)]
  }, "")
  if (anyNA(pick)) stop("There are several control conditions (", paste(ctrl, collapse = ", "), ") and StRym cannot tell which treated condition goes with which. ",
                        "Set the pairs in Analysis settings, e.g. “", treated[1], ":", ctrl[1], "”.", call. = FALSE)
  data.frame(treated = treated, control = unname(pick), stringsAsFactors = FALSE)
}

# Problems a time-course analysis would run into, found BEFORE the run from the sample sheet alone (shown in the page's checks).
# `sheet` has a "condition" column and, for time courses, "time_point". Returns list(pairs = data.frame or NULL, problems = character()).
strym_check_comparisons <- function(sheet, params, tp_names = character()) {
  out <- list(pairs = NULL, problems = character())
  if (!is.data.frame(sheet) || !"condition" %in% names(sheet)) return(out)
  S <- data.frame(condition = as.character(sheet$condition), time_point = if ("time_point" %in% names(sheet)) as.character(sheet$time_point) else "", stringsAsFactors = FALSE)
  pr <- tryCatch(strym_pairs(S, params), error = function(e) { out$problems <<- c(out$problems, conditionMessage(e)); NULL })
  out$pairs <- pr
  tps <- strym_order_timepoints(S$time_point, tp_names)
  if (length(tps) < 2) out$problems <- c(out$problems, sprintf("The sample sheet has %d time point%s: the time-course analyses need at least 2.", length(tps), if (length(tps) == 1) "" else "s"))
  else if (!is.null(pr)) for (cn in unique(c(pr$treated, pr$control))) {
    miss <- setdiff(tps, S$time_point[S$condition == cn])
    if (length(miss)) out$problems <- c(out$problems, sprintf("Condition “%s” has no samples at time point%s %s: that time point will be left out of the time-course analyses (every compared condition needs it).", cn, if (length(miss) == 1) "" else "s", paste(miss, collapse = ", ")))
  }
  out
}

# ---- the time course: mean expression, fold change and q-value per time point, for every treated/control pair -----------
# `tp_names` = time point names typed in the design (gives the order). Mirrors the data preparation of the baseline scripts.
strym_timecourse <- function(data, params, tp_names = character(), log = message) {
  S <- data$samples
  if (is.null(data$diff)) stop("This step needs the Cuffdiff comparisons (q-values): start from raw reads so StRym runs Tux.", call. = FALSE)
  tps <- strym_order_timepoints(S$time_point, tp_names)
  if (length(tps) < 2) stop("This analysis needs at least 2 time points; the sample sheet has ", length(tps), ".", call. = FALSE)
  pairs <- strym_pairs(S, params)
  # A time point where a compared condition has no samples (a failed sample, a time point only some conditions have) cannot be used for any
  # pair: it is left out, and the others are analyzed, instead of the whole step failing.
  has <- function(cond, t) any(S$condition == cond & S$time_point == t)
  all_tps <- tps
  usable <- vapply(tps, function(t) all(vapply(seq_len(nrow(pairs)), function(i) has(pairs$treated[i], t) && has(pairs$control[i], t), TRUE)), TRUE)
  if (!all(usable)) {
    log(sprintf("Left out time point%s %s: a compared condition has no samples there.", if (sum(!usable) == 1) "" else "s", paste(tps[!usable], collapse = ", ")))
    tps <- tps[usable]
    if (length(tps) < 2) stop(sprintf("Only %d time point%s has samples for every compared condition (left out: %s), but this analysis needs at least 2. Check the sample sheet.",
                                      length(tps), if (length(tps) == 1) "" else "s", paste(strym_order_timepoints(S$time_point, tp_names)[!usable], collapse = ", ")), call. = FALSE)
  }
  M <- strym_wide(data, na_bad = TRUE)
  genes <- rownames(M); diff <- data$diff
  gmean <- function(cond, tp) {
    cols <- S$sample[S$condition == cond & S$time_point == tp]
    if (!length(cols)) stop("No samples for condition “", cond, "” at time point “", tp, "”: every condition needs every time point.", call. = FALSE)
    v <- suppressWarnings(rowMeans(M[, cols, drop = FALSE], na.rm = TRUE)); v[is.na(v)] <- 0; v     # failed / all-NA = 0, as in the baseline scripts
  }
  grp <- function(cond, tp) unique(S$group[S$condition == cond & S$time_point == tp])[1]
  heat <- list(); ctrl <- list(); lfc <- list(); qval <- list()
  for (i in seq_len(nrow(pairs))) {
    tr <- pairs$treated[i]; ct <- pairs$control[i]
    h <- do.call(cbind, lapply(tps, function(t) gmean(tr, t))); c0 <- do.call(cbind, lapply(tps, function(t) gmean(ct, t)))
    dimnames(h) <- dimnames(c0) <- list(genes, tps)
    q <- matrix(1, length(genes), length(tps), dimnames = list(genes, tps))
    for (k in seq_along(tps)) {
      g1 <- grp(ct, tps[k]); g2 <- grp(tr, tps[k])
      d <- diff[(diff$sample_1 == g1 & diff$sample_2 == g2) | (diff$sample_1 == g2 & diff$sample_2 == g1), , drop = FALSE]
      if (!nrow(d)) stop("Cuffdiff has no comparison of “", g1, "” with “", g2, "” (found: ", paste(unique(paste(diff$sample_1, diff$sample_2, sep = " vs ")), collapse = "; "), ").", call. = FALSE)
      qq <- ifelse(d$status == "OK", d$q_value, 1); qq[is.na(qq)] <- 1
      keep <- d$test_id %in% genes
      q[d$test_id[keep], k] <- qq[keep]
    }
    heat[[i]] <- h; ctrl[[i]] <- c0; lfc[[i]] <- log2((h + 1) / (c0 + 1)); qval[[i]] <- q
  }
  top <- apply(do.call(cbind, c(heat, ctrl)), 1, max)
  universe <- genes[top >= params$min_fpkm]
  if (!length(universe)) stop(sprintf("No gene reaches %s >= %g in any sample (the highest value is %.3g). Lower “Minimum mean ... in some sample” in Analysis settings; values this small can also mean the table is on a log scale.",
                                      data$unit, params$min_fpkm, max(top, na.rm = TRUE)), call. = FALSE)
  log(sprintf("%d treated-vs-control pair%s (%s); %d time points (%s); %d genes reach %s >= %g somewhere.", nrow(pairs), if (nrow(pairs) == 1) "" else "s",
              paste(paste0(pairs$treated, " vs ", pairs$control), collapse = "; "), length(tps), paste(tps, collapse = ", "), length(universe), data$unit, params$min_fpkm))
  list(pairs = pairs, tps = tps, baseline_ok = identical(tps[1], all_tps[1]), heat = heat, ctrl = ctrl, lfc = lfc, qval = qval, universe = universe, genes = genes, samples = S,
       gene_names = data$gene_names, unit = data$unit)
}

# ---- results table shared by Shannon entropy and MSE ---------------------------------------------------------------------
.calls_table <- function(tc, i, calls, extra, main_min) {
  tps <- tc$tps; g <- rownames(calls); if (is.null(g)) g <- character()           # no gene called: an empty table, not an error
  nm <- if (!is.null(tc$gene_names)) unname(tc$gene_names[g]) else rep(NA_character_, length(g))
  up <- rowSums(calls == 1) > 0; dn <- rowSums(calls == -1) > 0
  base <- data.frame(gene = g, gene_name = nm, comparison = rep(paste0(tc$pairs$treated[i], "_vs_", tc$pairs$control[i]), length(g)), stringsAsFactors = FALSE)
  if (!is.null(extra)) base <- cbind(base, extra)
  res <- cbind(base, data.frame(n_sig_timepoints = rowSums(calls != 0), main_set = rowSums(calls != 0) >= main_min,
                                direction = ifelse(up & dn, "mixed", ifelse(up, "up", "down")), stringsAsFactors = FALSE),
               setNames(as.data.frame(calls), paste0("call_", tps)),
               setNames(as.data.frame(round(tc$lfc[[i]][g, , drop = FALSE], 3)), paste0("log2FC_", tps)),
               setNames(as.data.frame(signif(tc$qval[[i]][g, , drop = FALSE], 3)), paste0("q_", tps)))
  rownames(res) <- NULL
  res[order(!res$main_set, -res$n_sig_timepoints), , drop = FALSE]
}

# ---- Shannon entropy (Schug et al. 2005) of the mean time profile ----------------------------------------------------------
strym_entropy <- function(mat) {                                  # 1 = same level at every time point; lower = changes over time
  p <- mat / rowSums(mat)
  H <- -rowSums(ifelse(p > 0, p * log2(p), 0))
  ifelse(rowSums(mat) > 0, H / log2(ncol(mat)), 1)
}

# The name under which a pair's results are kept (and its file is named): the treated condition, or "<treated>_vs_<control>" when several pairs
# share the same treated condition (Heat vs Control and Heat vs Mock must not overwrite each other).
.pair_keys <- function(pairs) {
  shared <- duplicated(pairs$treated) | duplicated(pairs$treated, fromLast = TRUE)
  ifelse(shared, paste0(pairs$treated, "_vs_", pairs$control), pairs$treated)
}

# Returns list(results = list per pair, main = data.frame, excluded = n per pair).
strym_run_shannon <- function(tc, params) {
  U <- tc$universe; P <- seq_len(nrow(tc$pairs)); res <- list(); excl <- integer(0)
  if (!length(U)) return(list(results = list(), main = NULL, excluded = excl))
  H_ctrl <- do.call(cbind, lapply(P, function(i) strym_entropy(tc$ctrl[[i]][U, , drop = FALSE])))
  drift_i <- do.call(cbind, lapply(P, function(i) { x <- log2(tc$ctrl[[i]][U, , drop = FALSE] + 1); apply(abs(x - apply(x, 1, median)), 1, max) }))
  drift <- apply(ifelse(H_ctrl < params$h_control, drift_i, 0), 1, max); names(drift) <- U      # 0 = stable in every control
  for (i in P) {
    lfc <- tc$lfc[[i]]; sig <- abs(lfc) >= params$min_log2fc & tc$qval[[i]] < params$q_cutoff
    if (params$baseline_first && !isFALSE(tc$baseline_ok)) sig[, 1] <- FALSE                        # the baseline is never a call (unless the real baseline time point was left out: then no time point is blanked)
    heat_effect <- apply(ifelse(sig, abs(lfc), 0), 1, max)
    excluded <- U[drift[U] > 0 & heat_effect[U] < params$rescue_ratio * drift[U]]
    genes <- setdiff(U, excluded)
    H_heat <- if (length(genes)) strym_entropy(tc$heat[[i]][genes, , drop = FALSE]) else numeric(0)
    genes <- genes[H_heat < params$h_heat]
    calls <- sign(lfc[genes, , drop = FALSE]) * sig[genes, , drop = FALSE]
    calls <- calls[rowSums(calls != 0) > 0, , drop = FALSE]
    extra <- data.frame(H_heat = round(H_heat[rownames(calls)], 4), H_control = round(H_ctrl[match(rownames(calls), U), i], 4))
    res[[.pair_keys(tc$pairs)[i]]] <- .calls_table(tc, i, calls, extra, params$min_sig_tp); excl[i] <- length(excluded)
  }
  list(results = res, main = do.call(rbind, lapply(res, function(r) r[r$main_set, , drop = FALSE])), excluded = excl)
}

# ---- Modified Shannon entropy (MSE / ROKU-style outliers) ------------------------------------------------------------------
# For every gene: which time points are outliers (+1 high, -1 low, 0 normal)?
strym_mse_outliers <- function(mat, s_max = 3) {
  L <- ncol(mat)
  keep <- rowSums(mat != 0) > L / 2 + 0.5 & apply(mat, 1, sd) > 0
  mat <- mat[keep, , drop = FALSE]
  s_eff <- min(s_max, L - 2)                                      # at least 2 points must remain after removing outliers
  if (!nrow(mat) || s_eff < 1) return(matrix(0L, nrow(mat), L, dimnames = dimnames(mat)))
  pen <- function(n, s) sqrt(2) * s * log10(factorial(n)) / n
  out <- t(apply(mat, 1, function(x) {
    z <- as.numeric(scale(x)); ord <- order(z); zs <- z[ord]
    models <- list(integer(L)); U <- L * log10(sd(zs))            # model: no outliers
    for (s in seq_len(s_eff)) {
      lo <- integer(L); lo[ord[1:s]] <- -1L                       # the s lowest are outliers
      hi <- integer(L); hi[ord[(L - s + 1):L]] <- 1L              # the s highest are outliers
      models <- c(models, list(lo, hi))
      U <- c(U, (L - s) * log10(sd(zs[(s + 1):L])) + pen(L - s, s),
                (L - s) * log10(sd(zs[1:(L - s)])) + pen(L - s, s))
    }
    if (L >= 4) {
      both <- integer(L); both[ord[1]] <- -1L; both[ord[L]] <- 1L # one low + one high
      models <- c(models, list(both))
      U <- c(U, (L - 2) * log10(sd(zs[2:(L - 1)])) + pen(L - 2, 2))
    }
    U[is.nan(U) | U == Inf] <- Inf                                  # -Inf (identical remaining points = a perfect fit) stays the best model
    models[[which.min(U)]]
  }))
  dimnames(out) <- dimnames(mat)
  out
}

strym_run_mse <- function(tc, params) {
  U <- tc$universe; P <- seq_len(nrow(tc$pairs)); res <- list(); excl <- integer(0)
  if (!length(U)) return(list(results = list(), main = NULL, excluded = excl))
  drift <- apply(do.call(cbind, lapply(P, function(i) {          # control: genes that change over time without treatment
    o <- strym_mse_outliers(tc$ctrl[[i]][U, , drop = FALSE], params$s_max)
    x <- log2(tc$ctrl[[i]][rownames(o), , drop = FALSE] + 1)
    d <- setNames(numeric(length(U)), U)
    for (j in which(rowSums(o != 0) > 0)) {
      shift <- abs(x[j, ] - median(x[j, o[j, ] == 0]))[o[j, ] != 0]
      if (isTRUE(any(shift >= params$ctrl_log2fc, na.rm = TRUE))) d[rownames(o)[j]] <- max(shift, na.rm = TRUE)
    }
    d
  })), 1, max)                                                   # a gene varying in any control counts for every pair
  names(drift) <- U
  for (i in P) {
    lfc <- tc$lfc[[i]]; sig <- abs(lfc) >= params$min_log2fc & tc$qval[[i]] < params$q_cutoff
    if (params$baseline_first && !isFALSE(tc$baseline_ok)) sig[, 1] <- FALSE
    heat_effect <- apply(ifelse(sig, abs(lfc), 0), 1, max)
    excluded <- U[drift[U] > 0 & heat_effect[U] < params$rescue_ratio * drift[U]]
    genes <- setdiff(U, excluded)
    o  <- strym_mse_outliers(tc$heat[[i]][genes, , drop = FALSE], params$s_max)
    fc <- lfc[rownames(o), , drop = FALSE]
    calls <- o * ((o == 1 & fc >= params$min_log2fc) | (o == -1 & fc <= -params$min_log2fc))
    flip <- rowSums(calls != 0) == 0 & rowSums(o != 0) > 0          # sustained response: MSE flags the few unstressed points; the response is the rest
    if (any(flip)) calls[flip, ] <- (o[flip, , drop = FALSE] == 0) * ((fc[flip, , drop = FALSE] >= params$min_log2fc) - (fc[flip, , drop = FALSE] <= -params$min_log2fc))
    calls <- calls * sig[rownames(o), , drop = FALSE]
    calls <- calls[rowSums(calls != 0) > 0, , drop = FALSE]
    res[[.pair_keys(tc$pairs)[i]]] <- .calls_table(tc, i, calls, NULL, params$min_sig_tp); excl[i] <- length(excluded)
  }
  list(results = res, main = do.call(rbind, lapply(res, function(r) r[r$main_set, , drop = FALSE])), excluded = excl)
}

# Write the per-pair tables, the main-set table, the pairs used, and one overview picture. Returns the file names written.
strym_write_calls <- function(out_dir, prefix, run, tc, label, params = NULL) {
  files <- character()
  write_one <- function(df, f) { utils::write.csv(df, file.path(out_dir, f), row.names = FALSE); files <<- c(files, f) }
  for (nm in names(run$results)) write_one(run$results[[nm]], sprintf("%s_%s.csv", prefix, gsub("[^A-Za-z0-9._-]+", "_", nm)))
  main <- run$main; if (is.null(main)) main <- data.frame()
  write_one(main, sprintf("%s_main_genes.csv", prefix))
  write_one(data.frame(treated = tc$pairs$treated, control = tc$pairs$control, time_points = paste(tc$tps, collapse = ", "),
                       genes_removed_for_control_drift = run$excluded[seq_len(nrow(tc$pairs))],
                       genes_called = vapply(run$results, nrow, 1L), main_set = vapply(run$results, function(r) sum(r$main_set), 1L)),
            sprintf("%s_comparisons_used.csv", prefix))
  allres <- do.call(rbind, run$results)
  if (!is.null(allres) && nrow(allres) && requireNamespace("ggplot2", quietly = TRUE)) {
    cc <- grep("^call_", names(allres), value = TRUE)
    long <- do.call(rbind, lapply(cc, function(cn) data.frame(comparison = allres$comparison, time_point = sub("^call_", "", cn),
                                                              call = ifelse(allres[[cn]] > 0, "up", ifelse(allres[[cn]] < 0, "down", "none")), stringsAsFactors = FALSE)))
    long <- long[long$call != "none", , drop = FALSE]
    if (nrow(long)) {
      long$time_point <- factor(long$time_point, levels = tc$tps); long$call <- factor(long$call, levels = c("up", "down"))
      p <- ggplot2::ggplot(long, ggplot2::aes(time_point, fill = call)) + ggplot2::geom_bar(position = "dodge") +
        ggplot2::facet_wrap(~comparison) +
        ggplot2::scale_fill_manual(values = c(up = "#6F7D1C", down = "#4156A1"), drop = FALSE) +
        ggplot2::labs(x = "Time point", y = "Genes called", fill = NULL, title = paste0(label, ": genes called at each time point")) + ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10))                                    # a long title is no longer cut off in a narrow figure
      f <- sprintf("%s_calls_per_timepoint.png", prefix)
      ggplot2::ggsave(file.path(out_dir, f), p, width = 3 + 2.2 * min(4, nrow(tc$pairs)), height = 4, dpi = 150); files <- c(files, f)
    }
  }
  if (is.list(params) && isTRUE(params$entropy_heatmaps))                              # heatmaps are extras: a problem drawing them must never fail the step
    files <- c(files, tryCatch(strym_entropy_heatmaps(out_dir, prefix, run, tc, params, label),
                               error = function(e) { message("Note: the heatmaps could not be drawn (", conditionMessage(e), ")."); character() }))
  files
}

# ---- heatmaps of the genes an entropy step calls (pheatmap) ------------------------------------------------------------------
# Quick-look figures, not publication figures. Per comparison, two heatmaps of the called genes (the main set first, then the rest, at most
# params$heatmap_max_genes of them):
#   <prefix>_heatmap_log2FC_<comparison>.png       log2 fold change treated vs control at every time point (time points in order, not clustered)
#   <prefix>_heatmap_expression_<comparison>.png   expression of the treated and of the control profile, each gene scaled to its own mean (z-score)
# With several comparisons one more, <prefix>_heatmap_log2FC_all_comparisons.png, puts them side by side.
# The size follows the data: the figure grows with the number of genes and columns, gene names are printed only while they stay readable,
# and every awkward case (one gene, two time points, identical values, no variation, missing values) is drawn instead of failing.
.hm_clip <- function(x, lo, hi) max(lo, min(hi, x))
.hm_draw <- function(mat, file, title, col_info = NULL, row_info = NULL, gaps_col = NULL, diverging = TRUE, cluster_rows = TRUE, lim = NULL) {
  mat <- as.matrix(mat); mat[!is.finite(mat)] <- 0
  nr <- nrow(mat); nc <- ncol(mat); if (!nr || !nc) return(FALSE)
  if (is.null(lim)) lim <- max(1, as.numeric(stats::quantile(abs(mat), 0.98, na.rm = TRUE)))
  if (!is.finite(lim) || lim <= 0) lim <- 1
  mat <- pmax(pmin(mat, lim), -lim)                                                       # one outlier must not wash out the colours
  pal <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101)
  brk <- seq(-lim, lim, length.out = 102)                                                 # fixed, so a flat matrix still has valid breaks
  show_rn <- nr <= 80
  w <- .hm_clip(2.6 + 0.42 * nc + (if (show_rn) 2.2 else 0.4) + (if (!is.null(row_info)) 0.9 else 0), 5, 24)
  h <- .hm_clip(2.4 + 0.2 * nr + (if (!is.null(col_info)) 0.8 else 0.3), 4, 28)
  fr <- .hm_clip(11 - nr / 8, 4, 10); fc <- .hm_clip(11 - nc / 4, 5, 10)
  do_cluster <- cluster_rows && nr >= 2 && nr <= 3000 && any(apply(mat, 1, stats::sd) > 0)
  args <- list(mat, color = pal, breaks = brk, cluster_rows = do_cluster, cluster_cols = FALSE, show_rownames = show_rn,
               fontsize_row = fr, fontsize_col = fc, fontsize = 9, main = title, border_color = if (nr * nc <= 600) "grey85" else NA,
               annotation_col = col_info, annotation_row = row_info, gaps_col = gaps_col, filename = file, width = w, height = h)
  ac <- list()                                                                            # fixed colours for the side bars, so up / down / treated / control read the same everywhere
  if (!is.null(row_info) && "direction" %in% names(row_info)) ac$direction <- c(up = "#B2182B", down = "#2166AC", mixed = "#999999")
  if (!is.null(col_info) && "group" %in% names(col_info)) ac$group <- c(treated = "#B2182B", control = "#2166AC")
  if (length(ac)) args$annotation_colors <- ac
  args <- args[!vapply(args, is.null, TRUE)]
  if ("angle_col" %in% names(formals(pheatmap::pheatmap))) args$angle_col <- if (nc > 12) "90" else "45"
  ok <- tryCatch({ suppressMessages(do.call(pheatmap::pheatmap, args)); TRUE }, error = function(e) e)
  if (inherits(ok, "error") && do_cluster) {                                              # clustering can fail on odd data: draw without it
    args$cluster_rows <- FALSE
    ok <- tryCatch({ suppressMessages(do.call(pheatmap::pheatmap, args)); TRUE }, error = function(e) e)
  }
  if (inherits(ok, "error")) { message("Note: heatmap ", basename(file), " not drawn (", conditionMessage(ok), ")."); return(FALSE) }
  isTRUE(ok) && file.exists(file)
}

strym_entropy_heatmaps <- function(out_dir, prefix, run, tc, params, label) {
  if (!requireNamespace("pheatmap", quietly = TRUE) || !requireNamespace("RColorBrewer", quietly = TRUE)) return(character())
  if (!length(run$results)) return(character())
  maxg <- max(2L, as.integer(if (is.null(params$heatmap_max_genes)) 100 else params$heatmap_max_genes))
  keys <- .pair_keys(tc$pairs); files <- character(); tps <- tc$tps
  fname <- function(kind, key) sprintf("%s_heatmap_%s_%s.png", prefix, kind, gsub("[^A-Za-z0-9._-]+", "_", key))
  labs <- function(ids) { nm <- if (is.null(tc$gene_names)) NULL else unname(tc$gene_names[ids])
    l <- if (is.null(nm)) ids else ifelse(is.na(nm) | !nzchar(nm) | nm %in% c("-", "NA"), ids, sub("[,;].*$", "", nm)); make.unique(l, sep = "_") }
  top_genes <- function(r) { r <- r[order(!r$main_set, -r$n_sig_timepoints, -apply(abs(as.matrix(r[, grep("^log2FC_", names(r)), drop = FALSE])), 1, max), method = "radix"), , drop = FALSE]
    utils::head(r$gene, maxg) }
  rowinfo <- function(r, g) { x <- r[match(g, r$gene), , drop = FALSE]; df <- data.frame(direction = factor(x$direction, levels = c("up", "down", "mixed")),
                                                                                      time_points = x$n_sig_timepoints, row.names = labs(g)); df }
  dir_col <- list(direction = c(up = "#B2182B", down = "#2166AC", mixed = "#999999"))
  one_wide <- list()
  for (i in seq_len(nrow(tc$pairs))) {
    key <- keys[i]; r <- run$results[[key]]; if (is.null(r) || !nrow(r)) next
    g <- top_genes(r); if (!length(g)) next
    L <- tc$lfc[[i]][g, , drop = FALSE]; rownames(L) <- labs(g)
    ttl <- sprintf("%s: %s vs %s (%d gene%s%s)", label, tc$pairs$treated[i], tc$pairs$control[i], length(g), if (length(g) == 1) "" else "s",
                   if (nrow(r) > length(g)) sprintf(" of %d called", nrow(r)) else "")
    ri <- rowinfo(r, g); ann_col <- NULL
    if (.hm_draw(L, file.path(out_dir, fname("log2FC", key)), paste0(ttl, "\nlog2 fold change, treated vs control"), row_info = ri, lim = NULL,
                 cluster_rows = TRUE)) files <- c(files, fname("log2FC", key))
    # expression of treated and control, row-scaled (each gene: mean 0, sd 1 over the columns shown)
    H <- log2(tc$heat[[i]][g, , drop = FALSE] + 1); C <- log2(tc$ctrl[[i]][g, , drop = FALSE] + 1)
    E <- cbind(H, C); colnames(E) <- c(paste0(tc$pairs$treated[i], "_", tps), paste0(tc$pairs$control[i], "_", tps))
    mu <- rowMeans(E); sdv <- apply(E, 1, stats::sd); Z <- (E - mu) / ifelse(sdv > 0, sdv, 1); rownames(Z) <- labs(g)
    ann_col <- data.frame(group = factor(rep(c("treated", "control"), each = length(tps)), levels = c("treated", "control")), row.names = colnames(Z))
    if (.hm_draw(Z, file.path(out_dir, fname("expression", key)), paste0(ttl, "\nlog2 expression, scaled per gene (z-score)"), col_info = ann_col, row_info = ri,
                 gaps_col = length(tps), lim = 2.5, cluster_rows = TRUE)) files <- c(files, fname("expression", key))
    one_wide[[key]] <- list(genes = g, score = r$n_sig_timepoints[match(g, r$gene)], lfc = tc$lfc[[i]])
  }
  if (length(one_wide) > 1) {                                                              # all comparisons side by side
    sc <- tapply(unlist(lapply(one_wide, `[[`, "score")), unlist(lapply(one_wide, `[[`, "genes")), sum)
    g <- names(sort(sc, decreasing = TRUE, method = "radix"))[seq_len(min(maxg, length(sc)))]
    cols <- list(); gap <- integer(); n <- 0L
    for (key in names(one_wide)) { m <- one_wide[[key]]$lfc[g, , drop = FALSE]; colnames(m) <- paste0(key, "_", tps); cols[[key]] <- m; n <- n + ncol(m); gap <- c(gap, n) }
    A <- do.call(cbind, cols); rownames(A) <- labs(g)
    ann <- data.frame(comparison = factor(rep(names(cols), each = length(tps)), levels = names(cols)), row.names = colnames(A))
    if (.hm_draw(A, file.path(out_dir, fname("log2FC", "all_comparisons")), sprintf("%s: all comparisons (%d genes), log2 fold change", label, length(g)),
                 col_info = ann, gaps_col = utils::head(gap, -1))) files <- c(files, fname("log2FC", "all_comparisons"))
  }
  files
}
