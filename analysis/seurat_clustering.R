# analysis/seurat_clustering.R -- single cell / single nucleus RNA-seq with Seurat v5 (from the baseline notebooks
# Control-1.Rmd and Control_vs_FCP_int_SeuratV5.Rmd, generalized).
#
# Input  : the data chosen on the design page --
#            * a 10x matrix folder (barcodes / features / matrix .mtx, gz or plain), or a folder of SAMPLE folders that each hold
#              one (e.g. Control-1/raw_matrix, FCP1-1/raw_matrix ...): every sample is read and merged;
#            * a Seurat object (.rds).
#          A sample's condition comes from the sample sheet (columns sample, condition) or, without one, from its name with the
#          replicate number removed (Control-1 -> Control).
# Steps  : as in the notebooks -- Seurat object per sample (min.cells / min.features) -> UMI filter -> merge -> LogNormalize ->
#          variable features -> scale -> PCA -> number of PCs from the elbow rule (the smaller of: PC where the cumulative variance
#          passes 90% with < 5% each, and the last PC whose drop to the next is > 0.1%) -> integration across samples (CCA by default; RPCA,
#          Harmony or none) -> neighbors -> Louvain clusters -> UMAP and/or t-SNE (as ticked) -> cluster markers (wilcoxon,
#          only.pos, min.pct 0.10, logfc 0.25) and, with a control condition, treated-vs-control genes inside every cluster
#          (|log2FC| >= 0.25, adjusted p < 0.05; positive = higher in the treated condition; clusters need >= 10 cells in both).
#          One difference from the notebooks: only the variable features are scaled (marker tests use the normalized data anyway),
#          which keeps big datasets fast.
# Output : umap.png / tsne.png (clusters, sample, condition), qc_violin.png, elbow.png, cells_per_sample.csv, cluster_composition.csv,
#          cluster_markers.csv, condition_DEGs_by_cluster.csv, seurat_object.rds.
# (.find_10x, the sample labels inside barcodes and the other helpers shared with the design page are in analysis/singlecell.R)

# MatrixMarket file -> sparse matrix (genes x cells). With data.table installed the triplets are read in seconds and in a fraction of the memory of
# Matrix::readMM (72 million entries need about 3 GB instead of more than 7 GB); otherwise Matrix::readMM is used. Entries outside the matrix are dropped.
.read_mtx <- function(f) {
  gz <- grepl("\\.gz$", f); rd <- function() if (gz) gzfile(f, "rt") else file(f, "rt")
  slow <- function() { con <- rd(); on.exit(close(con)); methods::as(Matrix::readMM(con), "CsparseMatrix") }
  if (!requireNamespace("data.table", quietly = TRUE)) return(slow())
  con <- rd(); top <- tryCatch(readLines(con, n = 200), finally = close(con))
  hdr <- which(!startsWith(top, "%"))[1]
  if (!length(top) || is.na(hdr) || !grepl("coordinate", top[1]) || !grepl("general", top[1]) || !grepl("real|integer", top[1])) return(slow())
  dims <- suppressWarnings(as.numeric(strsplit(trimws(top[hdr]), "[[:space:]]+")[[1]]))
  if (length(dims) != 3 || anyNA(dims)) return(slow())
  src <- if (gz) list(cmd = paste("gzip -dc", shQuote(f))) else list(file = f)
  dt <- tryCatch(do.call(data.table::fread, c(src, list(skip = hdr, header = FALSE, sep = " ", colClasses = c("integer", "integer", "numeric"),
                                                         nThread = 1L, showProgress = FALSE))), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(dt) || ncol(dt) != 3) return(slow())
  i <- dt[[1]]; j <- dt[[2]]; x <- dt[[3]]; rm(dt)
  if (!is.numeric(x)) x <- suppressWarnings(as.numeric(x))
  ok <- !is.na(i) & !is.na(j) & i >= 1 & i <= dims[1] & j >= 1 & j <= dims[2]
  if (!all(ok)) { i <- i[ok]; j <- j[ok]; x <- x[ok] }
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = dims[1:2])
}

# Read one 10x matrix folder (matrix.mtx, features.tsv or genes.tsv, barcodes.tsv; each plain or .gz) into a genes x cells sparse matrix.
# Seurat's Read10X only accepts the all-gz layout; this accepts both. With several feature types only "Gene Expression" is kept.
.read_10x <- function(dir, log = function(...) invisible(NULL)) {
  pick <- function(...) { for (b in c(...)) for (f in c(b, paste0(b, ".gz"))) if (file.exists(file.path(dir, f))) return(file.path(dir, f)); NULL }
  mf <- pick("matrix.mtx"); ff <- pick("features.tsv", "genes.tsv"); bf <- pick("barcodes.tsv")
  if (is.null(mf) || is.null(ff) || is.null(bf)) stop("The 10x folder ", dir, " needs matrix.mtx, features.tsv (or genes.tsv) and barcodes.tsv (each may end in .gz).", call. = FALSE)
  rd <- function(f) if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt")
  m <- .read_mtx(mf)
  ft <- utils::read.delim(ff, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  con <- rd(bf); bc <- tryCatch(readLines(con), finally = close(con)); bc <- bc[nzchar(bc)]
  if (nrow(ft) != nrow(m) || length(bc) != ncol(m)) {
    if (nrow(ft) == ncol(m) && length(bc) == nrow(m)) {                      # some archives store cells x genes: turn it around
      m <- Matrix::t(m); log(sprintf("The matrix in %s was stored as cells x genes (%d x %d): it was turned around to genes x cells.", basename(dir), nrow(m), ncol(m)))
    } else stop("The 10x files in ", dir, " do not agree: matrix is ", nrow(m), " x ", ncol(m), ", features ", nrow(ft), ", barcodes ", length(bc), ".", call. = FALSE)
  }
  rownames(m) <- make.unique(as.character(if (ncol(ft) >= 2) ft[[2]] else ft[[1]])); colnames(m) <- bc
  if (ncol(ft) >= 3 && any(ft[[3]] == "Gene Expression")) m <- m[ft[[3]] == "Gene Expression", , drop = FALSE]
  m
}

# Counts that cannot be right - missing, infinite, negative or absurdly large (1e308 turns up in damaged files) - are set to 0 so that one bad entry
# cannot crash the whole analysis. Returns list(m, bad = entries changed, fractional = entries that are not whole numbers).
.clean_counts <- function(m, max_count = 1e7) {
  if (!inherits(m, "dgCMatrix")) m <- methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix")
  x <- m@x; bad <- !is.finite(x) | x < 0 | x > max_count
  n_bad <- sum(bad); n_frac <- sum(!bad & x != round(x))
  if (n_bad) { m@x[bad] <- 0; m <- Matrix::drop0(m) }
  list(m = m, bad = n_bad, fractional = n_frac, entries = length(x))
}

# The sample names listed in a sample sheet's "sample" column (NULL when there is no usable sheet).
.sheet_samples <- function(sample_sheet) {
  if (is.null(sample_sheet) || !nzchar(sample_sheet) || !file.exists(sample_sheet)) return(NULL)
  ss <- strym_read_sample_sheet(sample_sheet); if (inherits(ss, "error")) return(NULL)
  sc <- .find_col(ss, "^sample"); if (is.na(sc)) return(NULL)
  s <- trimws(as.character(ss[[sc]])); s[!is.na(s) & nzchar(s)]
}

# One matrix -> the sample matrices to analyze: counts cleaned; split by the sample label in the barcodes when the sample sheet lists those labels
# (labels the sheet does not list are left out; without a sheet all cells stay ONE sample); cells below the UMI filter dropped BEFORE any
# Seurat object is built (so empty samples cost nothing); a sample left with no cells is skipped, with the reason in the log.
# Returns list(mats = named list of matrices, read = cells read per sample, skipped = data.frame(sample, cells, reason)).
.sc_prepare <- function(m, name, sample_sheet, p, log) {
  nf <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  cl <- .clean_counts(m); m <- cl$m
  if (cl$bad) log(sprintf("%s: %s entries that were missing, infinite, negative or above 10 million were set to 0.", name, nf(cl$bad)))
  if (cl$fractional) log(sprintf("Note (%s): %s of %s entries are not whole numbers, so this does not look like raw counts (it may be normalized or corrected data). It is used as it is.", name, nf(cl$fractional), nf(cl$entries)))
  mats <- list(); labs <- strym_split_labels(colnames(m))
  if (length(labs)) {
    listed <- .sheet_samples(sample_sheet)
    use <- if (is.null(listed)) character() else labs[tolower(labs) %in% tolower(listed)]
    if (length(use)) {
      lab <- strym_barcode_labels(colnames(m))
      if (length(use) < length(labs)) log(sprintf("Left out (not in the sample sheet): %s.", paste(setdiff(labs, use), collapse = ", ")))
      for (L in use) {
        ms <- m[, which(lab == L), drop = FALSE]
        colnames(ms) <- strym_strip_label(colnames(ms), L)                                          # the label is added back when the samples are merged
        mats[[L]] <- ms
      }
      log(sprintf("%s holds %d samples; %s analyzed separately.", name, length(labs), if (length(use) == length(labs)) "all are" else sprintf("%d of them are", length(use))))
    } else {
      log(sprintf("Note: the cell barcodes in %s carry %d sample labels (%s%s), but %s, so all cells are analyzed as ONE sample. To analyze them separately, give StRym a sample sheet whose “sample” column lists these labels.",
                  name, length(labs), paste(head(labs, 4), collapse = ", "), if (length(labs) > 4) ", ..." else "",
                  if (is.null(listed)) "there is no sample sheet" else "the sample sheet does not list any of them"))
      mats[[name]] <- m
    }
  } else mats[[name]] <- m
  keep <- list(); read <- integer(); skipped <- data.frame(sample = character(), cells = integer(), reason = character(), stringsAsFactors = FALSE)
  for (s in names(mats)) {
    ms <- mats[[s]]; read[s] <- ncol(ms); u <- Matrix::colSums(ms)
    ok <- u > p$sc_min_counts & (p$sc_max_counts <= 0 | u < p$sc_max_counts)
    if (!any(ok)) {
      why <- sprintf("none of its %s cells has more than %g UMIs (the most any cell has is %s)", nf(ncol(ms)), p$sc_min_counts, nf(max(u, 0)))
      skipped[nrow(skipped) + 1, ] <- list(s, ncol(ms), why); log(sprintf("Sample %s skipped: %s.", s, why)); next
    }
    keep[[s]] <- ms[, ok, drop = FALSE]
  }
  list(mats = keep, read = read, skipped = skipped)
}

# A Seurat object for one prepared sample matrix, or NULL (with the reason in the log) when no cell passes the gene filter.
.sc_seurat <- function(ms, name, p, log) {
  o <- tryCatch(CreateSeuratObject(counts = ms, project = name, min.cells = p$sc_min_cells, min.features = p$sc_min_features), error = function(e) NULL)
  if (is.null(o) || ncol(o) == 0) { log(sprintf("Sample %s skipped: none of its cells has at least %g genes detected.", name, p$sc_min_features)); return(NULL) }
  o
}

# condition of every sample: from the sample sheet, else the name without its replicate number
.sc_conditions <- function(samples, sample_sheet) {
  cond <- sub("[-_.]?[0-9]+$", "", samples); cond[!nzchar(cond)] <- samples[!nzchar(cond)]
  if (nzchar(sample_sheet) && file.exists(sample_sheet)) {
    ss <- strym_read_sample_sheet(sample_sheet)
    if (!inherits(ss, "error")) {
      sc <- .find_col(ss, "^sample"); cc <- .find_col(ss, "^condition")
      if (!is.na(sc) && !is.na(cc)) { i <- match(tolower(samples), tolower(trimws(as.character(ss[[sc]])))); cond[!is.na(i)] <- trimws(as.character(ss[[cc]]))[i[!is.na(i)]] }
    }
  }
  cond
}

# Memory in use by this R process, for the log (Linux only). Big single cell data can use more memory than Ubuntu is given (7-8 GB by default);
# when a run is killed, the last of these lines shows how far it got and how much it was using.
.mem_note <- function(log, label) {
  st <- tryCatch(readLines("/proc/self/status", warn = FALSE), error = function(e) character()); mi <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character())
  kb <- function(x, key) { l <- grep(paste0("^", key, ":"), x, value = TRUE); if (length(l)) as.numeric(gsub("[^0-9]", "", l[1])) else NA_real_ }
  used <- kb(st, "VmRSS"); tot <- kb(mi, "MemTotal")
  if (!is.na(used) && !is.na(tot)) log(sprintf("Memory in use %s: %.1f GB of the %.1f GB Ubuntu has.", label, used / 1e6, tot / 1e6))
  invisible(NULL)
}
MARKER_CELLS <- 1000        # cell cap per group in the marker tests: clusters bigger than this are down-sampled (Seurat's usual advice for big data)

# Seurat's FindAllMarkers() one cluster at a time (same tests, same thresholds, same columns): the log shows the progress and the memory,
# memory is given back between clusters, and a cluster that fails is skipped instead of losing every cluster's result.
#
# Why not FindMarkers(<Seurat object>): it turns the WHOLE normalized matrix into a dense one before it picks the cells (20,000 genes x 22,000 cells = 3.7 GB for
# one test, whatever max.cells.per.ident says). So only the cells that are compared (at most MARKER_CELLS per side) are handed over, as a small matrix.
.de_wilcox <- function(dat, c1, c2, only_pos, min_pct) {
  if (length(c1) > MARKER_CELLS) c1 <- sample(c1, MARKER_CELLS)
  if (length(c2) > MARKER_CELLS) c2 <- sample(c2, MARKER_CELLS)
  sub <- dat[, c(c1, c2), drop = FALSE]                                    # (a plain matrix gives no result in FindMarkers: a small Seurat object it is)
  so <- CreateSeuratObject(counts = sub, project = "de", min.cells = 0, min.features = 0)
  so <- SeuratObject::SetAssayData(so, assay = "RNA", layer = "data", new.data = sub)
  SeuratObject::Idents(so) <- stats::setNames(factor(rep(c("a", "b"), c(length(c1), length(c2))), levels = c("a", "b")), colnames(so))
  suppressWarnings(FindMarkers(so, ident.1 = "a", ident.2 = "b", assay = "RNA", test.use = "wilcox", only.pos = only_pos, logfc.threshold = 0.25, min.pct = min_pct, verbose = FALSE))   # (ties in the counts make the Wilcoxon test warn "cannot compute exact p-value": it falls back to the normal approximation, which is fine for hundreds of cells)
}
.find_all_markers <- function(dat, clusters, log) {
  cl <- sort(unique(as.character(clusters))); res <- list(); set.seed(1)
  for (i in seq_along(cl)) {
    x <- tryCatch(.de_wilcox(dat, names(clusters)[clusters == cl[i]], names(clusters)[clusters != cl[i]], only_pos = TRUE, min_pct = 0.10), error = function(e) NULL)
    if (!is.null(x) && nrow(x)) { x <- x[x$p_val_adj < 0.05, , drop = FALSE]; if (nrow(x)) res[[length(res) + 1]] <- cbind(x, cluster = cl[i], gene = rownames(x), row.names = NULL) }
    rm(x); if (i %% 10 == 0 || i == length(cl)) { invisible(gc()); .mem_note(log, sprintf("after the marker tests of %d of %d clusters", i, length(cl))) }
  }
  if (length(res)) do.call(rbind, res) else data.frame()
}

# the number of PCs from the elbow rule of the notebooks
.elbow_pcs <- function(sdev) {
  pct <- sdev / sum(sdev) * 100; cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[-length(pct)] - pct[-1]) > 0.1), decreasing = TRUE)[1] + 1
  pcs <- suppressWarnings(min(c(co1, co2), na.rm = TRUE)); if (!is.finite(pcs)) pcs <- min(20, length(sdev))
  max(2, min(pcs, length(sdev)))
}

strym_step_seurat_clustering <- function(ctx) {
  for (pk in c("Seurat", "SeuratObject")) if (!requireNamespace(pk, quietly = TRUE))
    stop("The Seurat package is not installed in StRym's R environment. Run setup_modules.sh (see the README) to add it.", call. = FALSE)
  suppressPackageStartupMessages(library(Seurat))
  p <- strym_normalize_params(ctx$design$params); d <- ctx$design; out <- ctx$out_dir; log <- ctx$log
  want_umap <- "umap" %in% d$plots; want_tsne <- "tsne" %in% d$plots
  if (!want_umap && !want_tsne) want_umap <- TRUE                                      # run on its own: default to UMAP
  path <- ctx$data_path; if (!nzchar(path) || !file.exists(path)) stop("Choose your single cell / single nucleus data (a 10x folder or a Seurat .rds file).", call. = FALSE)

  # ---- 1. read the data ----
  files <- character()
  objs <- list(); before <- integer(); skipped <- data.frame(sample = character(), cells = integer(), reason = character(), stringsAsFactors = FALSE)
  add_matrix <- function(m, name) {                    # one matrix -> its samples: cleaned, split, UMI-filtered, then Seurat objects
    pr <- .sc_prepare(m, name, ctx$sample_sheet, p, log)
    skipped <<- rbind(skipped, pr$skipped); before <<- c(before, pr$read)
    for (s in names(pr$mats)) {
      o <- .sc_seurat(pr$mats[[s]], s, p, log)
      if (is.null(o)) skipped[nrow(skipped) + 1, ] <<- list(s, ncol(pr$mats[[s]]), "no cell has enough genes detected") else objs[[s]] <<- o
    }
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path); if (!inherits(obj, "Seurat")) stop("The .rds file is not a Seurat object.", call. = FALSE)
    obj <- tryCatch(suppressWarnings(SeuratObject::UpdateSeuratObject(obj)), error = function(e) obj)       # objects saved by an older Seurat version
    if (!"RNA" %in% SeuratObject::Assays(obj))
      stop(sprintf("The Seurat object has no “RNA” assay (it has: %s). StRym needs the raw counts in an assay called RNA.", paste(SeuratObject::Assays(obj), collapse = ", ")), call. = FALSE)
    SeuratObject::DefaultAssay(obj) <- "RNA"
    md <- colnames(obj@meta.data)                                      # (a missing column is an error in `obj$name`, so look at the names)
    if (!"sample" %in% md) obj$sample <- as.character(obj$orig.ident)  # a merged object usually has only orig.ident
    if (!"condition" %in% md) {
      cc <- intersect(c("source", "treatment", "group"), md)[1]
      if (!is.na(cc)) obj$condition <- as.character(obj@meta.data[[cc]])
      else { u <- unique(as.character(obj$sample)); obj$condition <- .sc_conditions(u, ctx$sample_sheet)[match(as.character(obj$sample), u)] }   # one condition per SAMPLE, then each cell gets its sample's
    }
    log(sprintf("Seurat object: %d cells, %d genes, %d samples.", ncol(obj), nrow(obj), length(unique(obj$sample))))
    before <- table(obj$sample)
  } else if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Reading a .h5 file needs the hdf5r package, which is not installed. Use the 10x matrix folder instead.", call. = FALSE)
    nm <- sub("\\.h5$", "", basename(path)); m <- Read10X_h5(path); if (is.list(m)) m <- m[[1]]
    add_matrix(m, nm); rm(m)
  } else if (!dir.exists(path)) {
    stop(sprintf("StRym cannot read “%s” as single cell data. It reads a 10x folder (matrix.mtx, barcodes.tsv and features.tsv, each plain or .gz), a 10x .h5 file, or a Seurat .rds object. %s", basename(path),
                 if (grepl("\\.(h5ad|loom)$", path, ignore.case = TRUE)) "For .h5ad or .loom files, export the counts as a 10x folder first (for example with scanpy's write_10x_mtx, or in R with SeuratDisk or zellkonverter)."
                 else if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", path, ignore.case = TRUE)) "A counts table with genes in rows and cells in columns can be turned into a 10x folder with Seurat's write10xCounts (DropletUtils) in R."
                 else ""), call. = FALSE)
  } else {
    smp <- .find_10x(path)
    if (!nrow(smp)) stop("No 10x matrix (matrix.mtx or matrix.mtx.gz with barcodes and features) found in that folder.", call. = FALSE)
    if (anyDuplicated(smp$sample)) smp$sample <- make.unique(smp$sample, sep = "_")
    log(sprintf("Reading %d 10x folder%s: %s", nrow(smp), if (nrow(smp) == 1) "" else "s", paste(smp$sample, collapse = ", ")))
    for (i in seq_len(nrow(smp))) { m <- .read_10x(smp$dir[i], log); add_matrix(m, smp$sample[i]); rm(m); invisible(gc()) }
  }
  if (!grepl("\\.rds$", path, ignore.case = TRUE)) {     # matrices: put the samples that survived together
    if (!length(objs)) stop(sprintf("Only 0 cells pass the UMI filter (%g to %s). Lower the minimum UMI count in Analysis settings (the log says what happened to each sample).",
                                    p$sc_min_counts, if (p$sc_max_counts > 0) p$sc_max_counts else "no maximum"), call. = FALSE)
    names(objs) <- make.unique(names(objs), sep = "_"); cond <- .sc_conditions(names(objs), ctx$sample_sheet)
    for (i in seq_along(objs)) { objs[[i]]$sample <- names(objs)[i]; objs[[i]]$condition <- cond[i] }
    log(sprintf("Samples: %s.", paste(sprintf("%s (%s)", names(objs), cond), collapse = ", ")))
    obj <- if (length(objs) == 1) objs[[1]] else merge(x = objs[[1]], y = objs[-1], add.cell.ids = names(objs))
    rm(objs); invisible(gc()); .mem_note(log, "after reading and merging the samples")
  }
  samples <- unique(as.character(obj$sample)); conds <- unique(as.character(obj$condition))
  if (length(conds) != d$n_conditions && length(conds) > 1) log(sprintf("Note: %d condition%s found in the data (%s) but %d set on the design page.", length(conds), if (length(conds) == 1) "" else "s", paste(conds, collapse = ", "), d$n_conditions))

  # ---- 2. quality control: UMI filter ----
  keep <- obj$nCount_RNA > p$sc_min_counts & (p$sc_max_counts <= 0 | obj$nCount_RNA < p$sc_max_counts)
  if (sum(keep) < 50) stop(sprintf("Only %d cells pass the UMI filter (%g to %s). Lower the minimum UMI count in Analysis settings.", sum(keep), p$sc_min_counts, if (p$sc_max_counts > 0) p$sc_max_counts else "no maximum"), call. = FALSE)
  obj <- subset(obj, cells = colnames(obj)[keep])
  qc <- do.call(rbind, lapply(samples, function(s) { x <- obj@meta.data[obj$sample == s, ]
    data.frame(sample = s, condition = as.character(x$condition[1]), cells_read = as.integer(before[s]), cells_kept = nrow(x),
               median_genes = median(x$nFeature_RNA), median_UMIs = median(x$nCount_RNA), stringsAsFactors = FALSE) }))
  if (nrow(skipped)) qc <- rbind(qc, data.frame(sample = skipped$sample, condition = .sc_conditions(skipped$sample, ctx$sample_sheet), cells_read = as.integer(skipped$cells),
                                                cells_kept = 0L, median_genes = NA_real_, median_UMIs = NA_real_, stringsAsFactors = FALSE))     # samples left out are listed too
  utils::write.csv(qc, file.path(out, "cells_per_sample.csv"), row.names = FALSE); files <- c(files, "cells_per_sample.csv")
  log(sprintf("QC: %d of %d cells kept (UMIs > %g%s)%s.", ncol(obj), sum(before), p$sc_min_counts, if (p$sc_max_counts > 0) sprintf(" and < %g", p$sc_max_counts) else "",
              if (nrow(skipped)) sprintf("; %d sample%s skipped (see cells_per_sample.csv)", nrow(skipped), if (nrow(skipped) == 1) "" else "s") else ""))
  ggplot2::ggsave(file.path(out, "qc_violin.png"), VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), group.by = "sample", pt.size = 0, ncol = 2),
                  width = 4 + 1.2 * length(samples), height = 4.5, dpi = 120); files <- c(files, "qc_violin.png")

  # ---- 3. normalize, variable genes, PCA, number of PCs ----
  multi <- length(samples) > 1 && !identical(p$sc_integration, "none")
  if (multi && length(Layers(obj[["RNA"]], search = "counts")) == 1) obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)     # one layer per sample for integration
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = min(2000, nrow(obj)), verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  npc <- min(50, ncol(obj) - 1, length(VariableFeatures(obj)) - 1)
  obj <- RunPCA(obj, npcs = npc, verbose = FALSE)
  pcs <- if (p$sc_npcs > 0) min(p$sc_npcs, npc) else .elbow_pcs(Stdev(obj, reduction = "pca"))
  log(sprintf("Using %d principal components (%s).", pcs, if (p$sc_npcs > 0) "set in Analysis settings" else "elbow rule")); .mem_note(log, "after PCA")
  ggplot2::ggsave(file.path(out, "elbow.png"), ElbowPlot(obj, ndims = npc) + ggplot2::geom_vline(xintercept = pcs, linetype = "dashed", colour = "#4156A1") +
                    ggplot2::labs(subtitle = sprintf("dashed line: %d PCs used", pcs)), width = 6, height = 4, dpi = 120); files <- c(files, "elbow.png")

  # ---- 4. integration, clusters, maps ----
  red <- "pca"
  if (multi) {
    method <- p$sc_integration
    min_n <- min(table(obj$sample))                                                       # CCA / RPCA weigh anchors by k.weight cells: it must fit the smallest sample
    if (min_n < 30) log(sprintf("Note: the smallest sample has only %d cells; integration may be unreliable.", min_n))
    dd <- max(2, min(30, npc, min_n - 1))                                                 # CCA / RPCA compare samples in `dims` dimensions: at most the cells of the smallest sample
    integrate <- function(method) do.call(IntegrateLayers, c(list(object = obj, method = switch(method, cca = CCAIntegration, rpca = RPCAIntegration, harmony = HarmonyIntegration),
                                                                   orig.reduction = "pca", new.reduction = paste0("integrated.", method), verbose = FALSE),
                                                              if (method %in% c("cca", "rpca")) list(k.weight = max(5, min(100, min_n - 1)), dims = 1:dd) else list()))
    res <- tryCatch(integrate(method), error = function(e) e)
    if (inherits(res, "error") && method != "harmony") {                                  # a method that cannot handle these samples must not stop the whole analysis
      log(sprintf("Note: %s integration did not work on these samples (%s); trying Harmony instead.", toupper(method), gsub("\n", " ", conditionMessage(res))))
      method <- "harmony"; res <- tryCatch(integrate(method), error = function(e) e)
    }
    if (inherits(res, "error")) log(sprintf("Integration did not work (%s): the samples are only merged.", gsub("\n", " ", conditionMessage(res))))
    else { obj <- res; red <- paste0("integrated.", method); log(sprintf("Integrated %d samples with %s.", length(samples), toupper(method))) }
    rm(res); invisible(gc()); .mem_note(log, "after integration")
  } else if (length(samples) > 1) log("Integration switched off: samples are only merged.")
  obj <- FindNeighbors(obj, dims = 1:pcs, reduction = red, verbose = FALSE)
  obj <- FindClusters(obj, resolution = p$sc_resolution, cluster.name = "clusters", verbose = FALSE)
  nclu <- length(unique(obj$clusters)); log(sprintf("%d clusters (resolution %g).", nclu, p$sc_resolution))
  maps <- character()
  if (want_umap) { obj <- RunUMAP(obj, dims = 1:pcs, reduction = red, reduction.name = "umap", verbose = FALSE); maps <- c(maps, "umap") }
  if (want_tsne) { obj <- RunTSNE(obj, dims = 1:pcs, reduction = red, reduction.name = "tsne", check_duplicates = FALSE,
                                  perplexity = max(2, min(30, floor((ncol(obj) - 1) / 3) - 1))); maps <- c(maps, "tsne") }
  for (r in maps) {
    pl <- list(DimPlot(obj, reduction = r, group.by = "clusters", label = TRUE, pt.size = 0.6) + ggplot2::ggtitle("Clusters"),
               DimPlot(obj, reduction = r, group.by = "sample", pt.size = 0.6) + ggplot2::ggtitle("Sample"))
    if (length(conds) > 1) pl[[3]] <- DimPlot(obj, reduction = r, group.by = "condition", pt.size = 0.6) + ggplot2::ggtitle("Condition")
    ggplot2::ggsave(file.path(out, paste0(r, ".png")), patchwork::wrap_plots(pl, nrow = 1), width = 5.2 * length(pl), height = 5, dpi = 120); files <- c(files, paste0(r, ".png"))
    if (length(conds) > 1) { ggplot2::ggsave(file.path(out, paste0(r, "_by_condition.png")),
                                             DimPlot(obj, reduction = r, group.by = "clusters", split.by = "condition", label = TRUE, pt.size = 0.6),
                                             width = 5 * length(conds), height = 5, dpi = 120); files <- c(files, paste0(r, "_by_condition.png")) }
  }
  comp <- as.data.frame.matrix(table(cluster = obj$clusters, sample = obj$sample)); comp <- cbind(cluster = rownames(comp), comp)
  utils::write.csv(comp, file.path(out, "cluster_composition.csv"), row.names = FALSE); files <- c(files, "cluster_composition.csv")

  # ---- 5. markers and condition effects ----
  # After integration the data sit in one layer per sample, and the marker tests need ONE matrix. JoinLayers() on the whole object (counts, data and scaled
  # data of every sample) needs several GB more than the object holds (it ran a 22,000-cell object out of 8 GB), so the tests get a small object made
  # from the normalized data alone. The saved object is joined too when there is memory for it.
  invisible(gc()); .mem_note(log, "after the maps")
  split_layers <- length(grep("^data", SeuratObject::Layers(obj[["RNA"]]))) > 1
  cell_of <- function(x) stats::setNames(as.character(x), colnames(obj))
  if (split_layers) {
    dl <- grep("^data", SeuratObject::Layers(obj[["RNA"]]), value = TRUE)
    mats <- lapply(dl, function(l) SeuratObject::LayerData(obj[["RNA"]], layer = l))
    common <- Reduce(intersect, lapply(mats, rownames))                    # a layer holds the genes found in ITS sample: join on the genes every sample has
    if (length(common) < max(vapply(mats, nrow, 1L))) log(sprintf("Note: the marker tests use the %s genes that every sample has.", format(length(common), big.mark = ",")))
    dat <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE])); rm(mats)
    mo <- CreateSeuratObject(counts = dat, project = "markers", min.cells = 0, min.features = 0)
    mo <- SeuratObject::SetAssayData(mo, assay = "RNA", layer = "data", new.data = dat); rm(dat)
    mo$clusters <- cell_of(obj$clusters)[colnames(mo)]; mo$condition <- cell_of(obj$condition)[colnames(mo)]
  } else mo <- obj
  Idents(mo) <- "clusters"
  marker_dat <- SeuratObject::LayerData(mo[["RNA"]], layer = "data"); marker_cl <- stats::setNames(as.character(mo$clusters), colnames(mo))
  if (length(conds) > 1) obj$cluster_condition <- paste(obj$clusters, obj$condition, sep = "__")
  # The Seurat object is SAVED NOW and then released: everything the tests still need (the normalized matrix, the cluster and condition of every cell)
  # is in `mo`, and the tests need the memory.
  n_cells <- ncol(obj); Idents(obj) <- "clusters"
  if (split_layers) {                                                      # join the layers of the saved object only when there is room for it
    need <- 2 * as.numeric(utils::object.size(obj[["RNA"]])) / 1e9
    mi <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character()); av <- grep("^MemAvailable:", mi, value = TRUE)
    avail <- if (length(av)) as.numeric(gsub("[^0-9]", "", av[1])) / 1e6 else Inf
    if (avail > need + 1) { obj <- JoinLayers(obj); .mem_note(log, "after joining the layers") }
    else log(sprintf("Note: the saved Seurat object keeps one layer per sample, because joining them needs about %.0f GB and %.0f GB are free. Run JoinLayers(obj) on a bigger computer before further differential expression.", need, avail))
  }
  saveRDS(obj, file.path(out, "seurat_object.rds")); files <- c(files, "seurat_object.rds"); rm(obj); invisible(gc()); .mem_note(log, "after saving the Seurat object")
  if (any(table(mo$clusters) > MARKER_CELLS)) log(sprintf("Note: marker tests use at most %s cells per group (larger clusters are sampled), which keeps big datasets fast and within memory.", format(MARKER_CELLS, big.mark = ",")))
  mk <- tryCatch(.find_all_markers(marker_dat, marker_cl, log), error = function(e) { log("Cluster markers could not be computed: ", conditionMessage(e)); NULL })
  if (!is.null(mk) && nrow(mk)) {
    mk$gene <- if ("gene" %in% names(mk)) mk$gene else rownames(mk)
    mk <- mk[order(as.character(mk$cluster), mk$p_val_adj, -mk$avg_log2FC), c("gene", "cluster", setdiff(names(mk), c("gene", "cluster")))]
    utils::write.csv(mk, file.path(out, "cluster_markers.csv"), row.names = FALSE); files <- c(files, "cluster_markers.csv")
  } else if (!is.null(mk)) log("No cluster markers passed the thresholds.")
  ndeg <- NA_integer_
  if (length(conds) > 1) {
    pairs <- tryCatch(strym_pairs(data.frame(condition = conds), p), error = function(e) { log("Condition effects skipped: ", conditionMessage(e)); NULL })
    if (!is.null(pairs)) {
      mo$cluster_condition <- paste(mo$clusters, mo$condition, sep = "__")
      res <- list(); n_cc <- table(mo$cluster_condition); cells_cc <- split(colnames(mo), mo$cluster_condition); set.seed(1); skipped <- character()
      for (i in seq_len(nrow(pairs))) for (k in sort(unique(as.character(mo$clusters)))) {
        i1 <- paste(k, pairs$treated[i], sep = "__"); i2 <- paste(k, pairs$control[i], sep = "__")
        if (is.na(n_cc[i1]) || is.na(n_cc[i2]) || n_cc[i1] < 10 || n_cc[i2] < 10) { skipped <- union(skipped, k); next }          # too few cells in a cluster to compare
        x <- tryCatch(.de_wilcox(marker_dat, cells_cc[[i1]], cells_cc[[i2]], only_pos = FALSE, min_pct = 0.01), error = function(e) NULL)
        if (!is.null(x) && nrow(x)) res[[length(res) + 1]] <- cbind(gene = rownames(x), cluster = k, comparison = paste0(pairs$treated[i], "_vs_", pairs$control[i]), x, row.names = NULL)
      }
      if (length(skipped)) log(sprintf("Clusters not compared between conditions (fewer than 10 cells of one of the conditions): %s. If the treatment itself separates cells into their own clusters, integrate the samples (CCA, RPCA or Harmony) to compare the conditions inside a cluster.", paste(skipped, collapse = ", ")))
      if (length(res)) { deg <- do.call(rbind, res); deg <- deg[deg$p_val_adj < 0.05, , drop = FALSE]                   # significant genes only
                         deg <- deg[order(deg$comparison, deg$cluster, deg$p_val_adj), ]; ndeg <- length(unique(deg$gene)); nrows_deg <- nrow(deg)
                         utils::write.csv(deg, file.path(out, "condition_DEGs_by_cluster.csv"), row.names = FALSE); files <- c(files, "condition_DEGs_by_cluster.csv")
                         log(sprintf("Condition effects: %d different genes (%d gene-by-cluster rows; a gene that changes in two clusters is one gene; positive log2FC = higher in the treated condition).", ndeg, nrows_deg))
      } else log("No genes differed between conditions inside the clusters (or clusters were too small to compare).")
    }
  }
  rm(mo, marker_dat); invisible(gc()); .mem_note(log, "at the end")
  list(status = "ok", files = files,
       message = sprintf("%d cells in %d samples (%d condition%s) -> %d clusters; %s%s.", n_cells, length(samples), length(conds), if (length(conds) == 1) "" else "s", nclu,
                         paste(toupper(maps), collapse = " + "), if (!is.na(ndeg)) sprintf("; %d condition-effect genes", ndeg) else ""))
}
