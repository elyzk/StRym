# analysis/singlecell.R -- helpers for single cell / single nucleus data that the design page and the Seurat step share.
# (Loaded together with common.R and registry.R by strym_load_registry().)

.TENX_DIRS <- c("raw_matrix", "filtered_matrix", "filtered_feature_bc_matrix", "raw_feature_bc_matrix", "matrix", "outs")

# The 10x matrix folders under `path` -> data.frame(sample, dir).
.find_10x <- function(path) {
  has <- function(d) length(list.files(d, pattern = "^matrix\\.mtx(\\.gz)?$")) > 0
  # The sample's folder: climb out of the matrix folders Cell Ranger and others use (Sample/outs/filtered_feature_bc_matrix -> Sample).
  sample_dir <- function(d) { while (basename(d) %in% .TENX_DIRS && dirname(d) != d) d <- dirname(d); d }
  if (has(path)) return(data.frame(sample = basename(sample_dir(path)), dir = path, stringsAsFactors = FALSE))
  m <- list.files(path, pattern = "^matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (!length(m)) return(data.frame(sample = character(), dir = character()))
  base <- normalizePath(path, winslash = "/"); mm <- normalizePath(m, winslash = "/")
  rel <- ifelse(startsWith(mm, base), substring(mm, nchar(base) + 2L), mm)                  # (plain string work: a path may contain ( ) + [ ] . ...)
  m <- m[lengths(strsplit(rel, "/", fixed = TRUE)) <= 4]                                      # at most 3 folders deep (Sample/outs/filtered_feature_bc_matrix/matrix.mtx.gz)
  if (!length(m)) return(data.frame(sample = character(), dir = character()))
  d <- dirname(m); sd <- unname(vapply(d, sample_dir, "")); out <- data.frame(sample = basename(sd), dir = d, sdir = sd, stringsAsFactors = FALSE, row.names = NULL)
  # One sample folder may hold the raw AND the filtered matrix (Cell Ranger writes both): use the filtered one (the raw one has every empty droplet).
  keep <- vapply(seq_len(nrow(out)), function(i) { sib <- which(out$sdir == out$sdir[i])
    length(sib) == 1 || { f <- sib[grepl("^filtered", basename(out$dir[sib]))]; if (length(f)) i == f[1] else i == sib[1] } }, logical(1))
  out <- out[keep, c("sample", "dir"), drop = FALSE]; rownames(out) <- NULL
  out
}

# The sample label inside a cell barcode (NA when there is none). A barcode is a run of at least 10 letters A/C/G/T/N (any case), optionally followed
# by the library number cellranger aggr adds (-1, -2, ...). The label sits in front of it or behind it, joined by one of  - _ . : # | + :
#   SAMN15315800-AAACACGACGCTCTTC -> "SAMN15315800"       Ctrl_AAACCTGAGAAACCAT -> "Ctrl"       Sample1#AAACCTGAGAAACCAT-1 -> "Sample1"
#   AAACCTGAGAAACCAT-1_SampleA -> "SampleA"               aaacctgagaaaccat.B2 -> "B2"
#   AAACCTGAGAAACCAT-2 (an aggr library number) -> "2"     AAACCTGAGAAACCAT -> NA (no label at all)
.BC_RE <- "[ACGTNacgtn]{10,}"; .SEP_RE <- "[-_.:#|+]"
strym_barcode_labels <- function(barcodes) {
  bc <- as.character(barcodes); ok <- !is.na(bc)
  has_tail <- ok & grepl(paste0(.BC_RE, "(-[0-9]+)?$"), bc, perl = TRUE)                  # barcode at the end: the label (if any) is in front
  pre <- ifelse(has_tail, sub(paste0(.SEP_RE, "?", .BC_RE, "(-[0-9]+)?$"), "", bc, perl = TRUE), "")
  lib <- ifelse(has_tail & grepl("-[0-9]+$", bc), sub("^.*-([0-9]+)$", "\\1", bc), "")
  head_re <- paste0("^", .BC_RE, "(?:-[0-9]+)?", .SEP_RE, "(.+)$")                          # barcode first: the label comes after it
  head_ok <- ok & !has_tail & grepl(head_re, bc, perl = TRUE)
  lab <- ifelse(nzchar(pre), pre, lib)
  lab[head_ok] <- sub(head_re, "\\1", bc[head_ok], perl = TRUE)
  lab[!ok | !nzchar(lab)] <- NA_character_
  lab
}

# `barcodes` without the sample `label` (and the separator next to it), whichever side it was on.
strym_strip_label <- function(barcodes, label) {
  bc <- as.character(barcodes); out <- bc; pre <- startsWith(bc, label); suf <- !pre & endsWith(bc, label)
  out[pre] <- sub(paste0("^", .SEP_RE), "", substring(bc[pre], nchar(label) + 1L))
  out[suf] <- sub(paste0(.SEP_RE, "$"), "", substr(bc[suf], 1L, nchar(bc[suf]) - nchar(label)))
  out
}

# The distinct sample labels (in order of first appearance) when ONE barcode list holds several samples; character(0) when it holds one
# sample, when some barcodes have no label, or when the "labels" are nearly one per cell (so they are not sample names).
strym_split_labels <- function(barcodes) {
  lab <- strym_barcode_labels(barcodes)
  if (!length(lab) || anyNA(lab)) return(character())
  u <- unique(lab)
  if (length(u) < 2 || length(u) > length(lab) / 5) return(character())
  u
}

# The samples of a 10x data folder: one per matrix folder, or - for ONE matrix whose barcodes carry sample labels - one per label.
# Returns list(samples, how = "folders" | "barcodes" | "none").
strym_10x_samples <- function(path) {
  f <- .find_10x(path)
  if (!nrow(f)) return(list(samples = character(), how = "none"))
  if (nrow(f) > 1) return(list(samples = f$sample, how = "folders"))
  bf <- list.files(f$dir[1], pattern = "^barcodes\\.tsv(\\.gz)?$", full.names = TRUE)[1]
  labs <- character()
  if (!is.na(bf)) {
    con <- if (grepl("\\.gz$", bf)) gzfile(bf, "rt") else file(bf, "rt")
    bc <- tryCatch(readLines(con), error = function(e) character(), finally = close(con))
    labs <- strym_split_labels(bc[nzchar(bc)])
  }
  if (length(labs)) list(samples = labs, how = "barcodes") else list(samples = f$sample[1], how = "folders")
}
