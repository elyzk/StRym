# pipeline.R -- StRym's automation engine: experimental design -> plan -> organized results folder.
#
#   design  = description of the experiment (strym_default_design())
#   plan    = which steps make sense for it (strym_plan(), rules in analysis/registry.R)
#   run     = strym_run(): runs the chosen steps, each writing into its own numbered folder of the project,
#             then writes run_summary.csv and index.html at the top of the project folder.

# (re)load the registry and the placeholder helper from <app_dir>/analysis
strym_load_registry <- function(app_dir) {
  source(file.path(app_dir, "analysis", "placeholder.R"), local = FALSE)
  source(file.path(app_dir, "analysis", "common.R"), local = FALSE)
  source(file.path(app_dir, "analysis", "derived_diff.R"), local = FALSE)
  source(file.path(app_dir, "analysis", "singlecell.R"), local = FALSE)
  source(file.path(app_dir, "analysis", "registry.R"), local = FALSE)
  invisible(TRUE)
}

# ---- the design -------------------------------------------------------------------------------------------
strym_default_design <- function() list(
  project_name = "", project_dir = "",
  data_type = "bulk",              # "bulk" (whole organ / tissue) or "single_cell"
  data_kind = "counts",            # bulk: "counts" | "fastq";  single cell: "10x" | "seurat"
  data_path = "", sample_sheet = "",
  expr_unit = "auto",              # bulk counts table: "auto" | "counts" (raw counts -> CPM) | "normalized" (FPKM / TPM: used as given)
  diff_path = character(),         # bulk counts table: Cuffdiff gene_exp.diff file(s) / folder / .zip made elsewhere (q-values for the entropy steps)
  fastq_files = character(),       # raw reads: the FASTQ files chosen (all in data_path); empty = every FASTQ file in that folder
  gtf_file = "", fasta_file = "",  # only for raw reads (quantification): annotation and genome
  read_type = "single",            # "single" | "paired"  (raw reads)
  adapter_file = "",               # raw short reads: a FASTA of adapter sequences to trim (blank = the list that comes with StRym: TruSeq, Nextera, MGI/BGI, small RNA ...)
  read_platform = "short",         # raw reads: "short" (Illumina, MGI ...: HISAT2 + Cufflinks) | "long" (Oxford Nanopore, PacBio, e.g. Plasmidsaurus: minimap2 + featureCounts)
  threads = 4L,
  n_conditions = 2L, condition_names = character(),
  n_timepoints = 1L, timepoint_names = character(),
  de_methods = character(),        # subset of deseq2, shannon_entropy, modified_shannon_entropy
  plots = character(),             # subset of pca, umap, tsne, trajectory
  grn_methods = character(),       # gene regulatory network inference: subset of genist, rtpstar
  params = list(),                 # analysis settings (see strym_default_params() in analysis/common.R)
  interactive_browser = FALSE)

# ---- analysis settings ----------------------------------------------------------------------------------------------
# Defaults are the ones of the baseline scripts. Everything can be changed in "Analysis settings" on the design page.
strym_default_params <- function() list(
  control = "",                  # name(s) of the control condition(s), comma separated; blank = detect from control_pattern
  control_pattern = "control|ctrl|(^|[^a-z])ctl([^a-z]|$)|mock|untreated|wild.?type|^wt$|^ck$|baseline|reference",
  pairs = "",                    # advanced: "Heat:Control, Heat2:Control2" (treated:control); overrides automatic pairing
  baseline_first = TRUE,         # the first time point is a baseline (before treatment): reported, never called
  min_fpkm = 2,                  # a gene must reach this mean FPKM in some sample
  min_log2fc = 1,                # a time point needs |log2((treated+1)/(control+1))| >= this ...
  q_cutoff = 0.05,               # ... and Cuffdiff q < this
  min_sig_tp = 2,                # "main set": significant at >= this many time points
  ctrl_log2fc = 0.585,           # MSE: a control outlier must differ >= 1.5-fold from the rest
  rescue_ratio = 2,              # keep a gene that varies in the control if its treatment effect is >= this x that drift
  s_max = 3,                     # MSE: maximum number of outlier time points per gene
  h_control = 0.99,              # Shannon: control entropy below this = changes over time in the control -> excluded
  h_heat = 0.99,                 # Shannon: treated entropy must be below this = changes over time
  drop_reps = "",                # samples to leave out of every step (comma separated sample names, e.g. PCA outliers)
  tux_fresh = FALSE,             # TRUE = ignore the finished parts of an earlier Tux run and start it from scratch
  pca_genes = 500, pca_min_mean = 1,
  browser_extra = "",            # interactive browser: more result folders to include (comma separated), e.g. a GENIST / RTP-STAR tab run saved outside the project
  pca_level = "replicate",       # PCA: "replicate" (every sample is a point) | "average" (replicates averaged per condition and time point) | "both"
  entropy_heatmaps = TRUE,       # Shannon / modified Shannon entropy: also draw a pheatmap of the called genes (log2 fold change per time point)
  heatmap_max_genes = 100,       # ... at most this many genes per comparison (the strongest first; the heatmap shrinks its labels as genes grow)
  genist_max_genes = 500,        # GENIST step: at most this many genes in the network (strongest entropy genes first; 0 = no limit)
  genist_tf_file = "",           # GENIST step: list of transcription-factor ids (first column); blank = every gene may regulate (slow for many genes)
  trajectory_genes = "", trajectory_max = 24,
  sc_min_features = 200, sc_min_cells = 3, sc_min_counts = 500, sc_max_counts = 0,   # single cell: 0 = no upper limit
  sc_resolution = 1, sc_integration = "cca", sc_npcs = 0)                             # sc_npcs 0 = choose automatically

strym_split_list <- function(x) { x <- unlist(strsplit(as.character(x), "[,;\n]")); x <- trimws(x); x[nzchar(x)] }

# Fill in defaults and coerce types. Safe on partial input (e.g. read back from JSON).
strym_normalize_params <- function(p) {
  base <- strym_default_params(); p <- as.list(p)
  for (k in names(base)) {
    v <- p[[k]]
    if (is.character(base[[k]]) && length(v) > 1) v <- paste(as.character(unlist(v)), collapse = ", ")     # a list given as a JSON array: "A_1", "B_1" -> "A_1, B_1"
    if (is.null(v) || length(v) != 1 || is.na(v)) { p[[k]] <- base[[k]]; next }
    p[[k]] <- if (is.logical(base[[k]])) isTRUE(as.logical(v)) || identical(trimws(as.character(v)), "1")
              else if (is.numeric(base[[k]])) { n <- suppressWarnings(as.numeric(v)); if (is.na(n)) base[[k]] else n }
              else trimws(as.character(v))
  }
  if (!p$sc_integration %in% c("none", "cca", "rpca", "harmony")) p$sc_integration <- "cca"
  p$pca_level <- tolower(p$pca_level); if (!p$pca_level %in% c("replicate", "average", "both")) p$pca_level <- "replicate"
  if (!is.finite(p$heatmap_max_genes) || p$heatmap_max_genes < 2) p$heatmap_max_genes <- base$heatmap_max_genes
  p[names(base)]
}


.split_names <- function(x) { x <- unlist(strsplit(as.character(x), "[,;\n]")); x <- trimws(x); x[nzchar(x)] }

# Fill in defaults, coerce types, tidy names. Safe to call on partial input (e.g. from JSON).
strym_normalize_design <- function(d) {
  base <- strym_default_design()
  for (k in names(base)) if (is.null(d[[k]])) d[[k]] <- base[[k]]
  d$n_conditions <- suppressWarnings(as.integer(d$n_conditions)); if (is.na(d$n_conditions)) d$n_conditions <- 0L
  d$n_timepoints <- suppressWarnings(as.integer(d$n_timepoints)); if (is.na(d$n_timepoints)) d$n_timepoints <- 0L
  d$condition_names <- .split_names(d$condition_names); d$timepoint_names <- .split_names(d$timepoint_names)
  d$de_methods <- as.character(unlist(d$de_methods)); d$plots <- as.character(unlist(d$plots))
  d$grn_methods <- as.character(unlist(d$grn_methods)); d$fastq_files <- as.character(unlist(d$fastq_files))
  d$params <- strym_normalize_params(d$params)
  d$diff_path <- trimws(as.character(unlist(d$diff_path))); d$diff_path <- d$diff_path[!is.na(d$diff_path) & nzchar(d$diff_path)]
  d$expr_unit <- tolower(trimws(as.character(d$expr_unit)[1])); if (is.na(d$expr_unit) || !d$expr_unit %in% c("auto", "counts", "normalized")) d$expr_unit <- "auto"
  d$interactive_browser <- isTRUE(as.logical(d$interactive_browser))
  for (k in c("project_name", "project_dir", "data_type", "data_kind", "data_path", "sample_sheet", "gtf_file", "fasta_file", "read_type", "adapter_file"))
    d[[k]] <- trimws(as.character(d[[k]]))
  if (!d$read_type %in% c("single", "paired")) d$read_type <- "single"
  d$read_platform <- tolower(trimws(as.character(d$read_platform)[1])); if (is.na(d$read_platform) || !d$read_platform %in% c("short", "long")) d$read_platform <- "short"
  d$threads <- suppressWarnings(as.integer(d$threads)); if (is.na(d$threads) || d$threads < 1) d$threads <- 1L
  d
}

# The result files of a project that StRym shows (the Results list of the design page and the interactive browser): only what the analysis steps
# wrote into their own numbered folders (01_..., 02_..., ...), never the other files of the output folder (a researcher may well pick Downloads
# or Documents as the output folder), nor the design record, the browser itself or the Tux internals.
strym_result_files <- function(project_dir, pattern) {
  if (!nzchar(project_dir) || !dir.exists(project_dir)) return(character())
  tops <- list.dirs(project_dir, recursive = FALSE, full.names = FALSE)
  tops <- tops[grepl("^[0-9]{2}_", tops) & !tops %in% c("00_design", "07_interactive_browser")]
  f <- unlist(lapply(tops, function(d) file.path(d, list.files(file.path(project_dir, d), pattern = pattern, recursive = TRUE, ignore.case = TRUE))))
  f <- f[!grepl("^01_quantification/tux/", f)]
  if (is.null(f)) character() else sort(f)
}
# a table file of network results written by GENIST / RTP-STAR (a .txt, not .csv / .tsv)
STRYM_NETWORK_TXT <- "(network|cytoscape|biograph).*\\.txt$"

# Where the quantification step leaves its counts table (used by every later step when starting from raw reads).
strym_counts_path <- function(design) file.path(design$project_dir, "01_quantification", "counts.csv")

# ---- the plan -----------------------------------------------------------------------------------------------
# One row per automatic step: id, label, group, desc, dir, applicable, default, note.
strym_plan <- function(design) {
  d <- strym_normalize_design(design)
  rows <- lapply(strym_steps(), function(s) {
    r <- s$rule(d)
    data.frame(id = s$id, label = s$label, group = s$group, desc = s$desc, dir = s$dir,
               applicable = isTRUE(r$applicable), default = isTRUE(r$default) && isTRUE(r$applicable),
               fit = isTRUE(r$applicable) && isTRUE(r$fit), why = if (is.null(r$why)) "" else r$why,
               note = r$note, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Steps this one builds on that already have results (an entropy step's *_main_genes.csv) in the project folder from an earlier run, so that
# "Run only this" works for GENIST or trajectories without running the entropy step again.
strym_existing_prereq <- function(d, s) {
  if (is.null(d$project_dir) || !nzchar(d$project_dir) || identical(s$needs_any, "*")) return(character())
  reg <- strym_steps(); ids <- vapply(reg, `[[`, "", "id"); dirs <- vapply(reg, `[[`, "", "dir")
  want <- intersect(s$needs_any, ids)
  want[vapply(want, function(i) { p <- file.path(d$project_dir, dirs[ids == i]); dir.exists(p) && length(list.files(p, pattern = "_main_genes\\.csv$")) > 0 }, TRUE)]
}

# For the chosen step ids: named character vector (names = ids) with the reason a step cannot run ("" = fine).
strym_blocked <- function(design, selected) {
  d <- strym_normalize_design(design); steps <- strym_steps(); lab <- setNames(vapply(steps, `[[`, "", "label"), vapply(steps, `[[`, "", "id"))
  out <- setNames(rep("", length(selected)), selected)
  have_counts <- !identical(d$data_kind, "fastq") || "quantify" %in% selected ||
    (nzchar(d$project_dir) && file.exists(strym_counts_path(d)))      # raw reads need counts first
  for (s in steps) if (s$id %in% selected) {
    r <- s$rule(d)
    if (!isTRUE(r$applicable)) { out[s$id] <- if (nzchar(r$note)) r$note else "Does not apply to this experiment."; next }
    if (s$id != "quantify" && !have_counts) { out[s$id] <- "Needs the counts from “Quantify reads”: tick that step (or run it first)."; next }
    if (length(s$needs_any) && !(is.function(s$needs_unless) && isTRUE(s$needs_unless(d)))) {
      pool <- if (identical(s$needs_any, "*")) setdiff(selected, s$id) else intersect(s$needs_any, selected)
      if (!length(pool) && !length(strym_existing_prereq(d, s))) out[s$id] <- if (identical(s$needs_any, "*")) "Needs at least one analysis step to be selected."
        else paste0("Needs one of: ", paste(lab[s$needs_any], collapse = ", "), " (tick one of them, or give a gene list in Analysis settings).")
    }
  }
  out
}

# Hand-off modules (their own tabs): id, label, tab, desc, applicable, note.
strym_module_plan <- function(design) {
  d <- strym_normalize_design(design)
  do.call(rbind, lapply(strym_modules(), function(m) {
    r <- m$rule(d)
    data.frame(id = m$id, label = m$label, tab = m$tab, desc = m$desc, applicable = isTRUE(r$applicable), note = r$note,
               fit = isTRUE(r$applicable) && isTRUE(r$fit), why = if (is.null(r$why)) "" else r$why, stringsAsFactors = FALSE)
  }))
}

# Network inference (GENIST / RTP-STAR) needs a list of genes: warn when no gene-selection method is ticked.
strym_check_networks <- function(design, selected) {
  d <- strym_normalize_design(design)
  if (!length(d$grn_methods)) return(character())
  labs <- c(genist = "GENIST", rtpstar = "RTP-STAR")
  if (length(intersect(c("deseq2", "shannon_entropy", "modified_shannon_entropy"), selected))) return(character())
  sprintf("%s needs a list of genes: choose a differential expression / entropy method above, or bring your own gene list to its tab.",
          paste(labs[d$grn_methods], collapse = " and "))
}

# Problems with the plan itself, given the user's final selection. Returns list(errors, warnings).
strym_check_selection <- function(design, selected) {
  errs <- character(); warns <- character()
  if (!length(selected)) return(list(errors = "Select at least one analysis step.", warnings = warns))
  b <- strym_blocked(design, selected); b <- b[nzchar(b)]
  lab <- setNames(vapply(strym_steps(), `[[`, "", "label"), vapply(strym_steps(), `[[`, "", "id"))
  for (i in names(b)) warns <- c(warns, sprintf("%s: %s It will be skipped.", lab[[i]], b[[i]]))
  if (length(b) == length(selected)) errs <- c(errs, "None of the selected steps can run for this experiment.")
  if ("genist_network" %in% selected && !"genist_network" %in% names(b) && !nzchar(strym_normalize_params(design$params)$genist_tf_file))
    warns <- c(warns, "GENIST without a transcription factor list (optional): every gene may regulate every other, so the network takes longer to infer, much longer beyond about 100 genes. Add a list in Analysis settings to speed it up.")
  list(errors = errs, warnings = warns)
}

# ---- checking the experiment -----------------------------------------------------------------------------------
# Quick look at a data file/folder for the checks panel. Returns a one-line description or an error string with attr "bad".
# Counting the rows of a big table takes seconds (the page asks again at every edit): each file is looked at once.
.peek_cache <- new.env()
strym_peek_data <- function(path) {
  info <- file.info(path); key <- paste(path, info$size, format(info$mtime, "%s"))
  if (!is.na(info$size) && !is.null(hit <- .peek_cache[[key]])) return(hit)
  res <- strym_peek_data_uncached(path); if (!is.na(info$size)) assign(key, res, envir = .peek_cache); res
}
strym_peek_data_uncached <- function(path) {
  bad <- function(msg) structure(msg, bad = TRUE)
  if (!file.exists(path)) return(bad("not found"))
  if (dir.exists(path)) return(sprintf("folder with %d items", length(list.files(path))))
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
  size <- format(structure(file.size(path), class = "object_size"), units = "auto")
  if (ext %in% c("csv", "tsv", "txt", "tab")) {
    first <- readLines(path, n = 1, warn = FALSE)
    sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (lengths(regmatches(first, gregexpr(";", first, fixed = TRUE))) > lengths(regmatches(first, gregexpr(",", first, fixed = TRUE)))) ";" else ","
    ncols <- length(strsplit(first, sep, fixed = TRUE)[[1]])
    if (grepl("\\.gz$", path, ignore.case = TRUE)) return(sprintf("compressed table with %d columns (%s)", ncols, size))      # counting rows would mean unpacking it
    nrows <- suppressWarnings(as.integer(sub("\\s.*", "", system2("wc", c("-l", shQuote(path)), stdout = TRUE)))) - 1L
    return(sprintf("table: %s rows × %d columns (%s)", format(nrows, big.mark = ","), ncols, size))
  }
  if (ext %in% c("xlsx", "xls")) return(sprintf("Excel workbook (%s)", size))
  sprintf("%s file (%s)", if (nzchar(ext)) toupper(ext) else "data", size)
}

# What the numbers of a counts table are, for the plan (which steps can run): the unit the researcher chose, or, for "work it out", what the first rows say.
# A table whose values have decimals (FPKM / TPM / CPM) is "normalized"; whole numbers are taken as "auto" (= raw counts). Looked at once per file.
# Excel files are not peeked at (reading them is slow): the run itself decides then.
.unit_cache <- new.env()
.guess_unit <- function(path) {
  if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path) || !file.exists(path) || dir.exists(path)) return("auto")
  info <- file.info(path); key <- paste(path, info$size, format(info$mtime, "%s")); hit <- .unit_cache[[key]]
  if (!is.null(hit)) return(hit)
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
  res <- "auto"
  if (ext %in% c("csv", "tsv", "txt", "tab")) res <- tryCatch({
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r"); l <- readLines(con, n = 401, warn = FALSE); close(con)
    tf <- tempfile(fileext = paste0(".", ext)); writeLines(l, tf); on.exit(unlink(tf), add = TRUE)
    tab <- read_table_any(tf); v <- unlist(lapply(tab[-1], .to_num), use.names = FALSE); v <- v[is.finite(v)]
    if (length(v) && mean(v != floor(v)) >= 0.001) "normalized" else "auto"
  }, error = function(e) "auto")
  assign(key, res, envir = .unit_cache); res
}
strym_unit_of <- function(d) if (identical(d$data_kind, "counts") && identical(d$expr_unit, "auto")) .guess_unit(d$data_path) else d$expr_unit

# Read a sample sheet (csv/tsv/xlsx) or return an error object.
strym_read_sample_sheet <- function(path) {
  tryCatch({
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx", "xls")) as.data.frame(readxl::read_excel(path, .name_repair = "minimal", col_types = "text"))
    else { first <- readLines(path, n = 1, warn = FALSE)
           sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (lengths(regmatches(first, gregexpr(";", first, fixed = TRUE))) > lengths(regmatches(first, gregexpr(",", first, fixed = TRUE)))) ";" else ","
           utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", colClasses = "character", na.strings = character()) }
  }, error = function(e) e)
}
.find_col <- function(df, pattern) { i <- grep(pattern, names(df), ignore.case = TRUE); if (length(i)) names(df)[i[1]] else NA_character_ }

# The output-folder rules that must hold even when a run is started without the page (command line, scheduler): never a file, never a drive root, never the input folder.
.project_dir_errors <- function(d) {
  errs <- character(); pd <- d$project_dir; if (!nzchar(pd)) return(errs)
  if (file.exists(pd) && !dir.exists(pd)) errs <- c(errs, "The project folder you chose is a file, not a folder.")
  if (grepl("^/+\\.?/*$", pd) || identical(normalizePath(pd, mustWork = FALSE), "/") || grepl("^[A-Za-z]:[\\\\/]*$", pd)) errs <- c(errs, "Choose a folder for the results, not the root of the drive.")
  if (nzchar(d$data_path) && dir.exists(d$data_path) && identical(normalizePath(pd, mustWork = FALSE), normalizePath(d$data_path, mustWork = FALSE)))
    errs <- c(errs, "The output folder is the folder that holds your input data. Choose a different output folder so your data and the results are kept apart.")
  errs
}

# list(errors, warnings, info) about the experiment description (not about which steps are chosen).
strym_validate_design <- function(design, check_paths = TRUE) {
  d <- strym_normalize_design(design); errs <- character(); warns <- character(); info <- character()
  if (!nzchar(d$project_name)) errs <- c(errs, "Give the project a name.")
  if (!nzchar(d$project_dir)) errs <- c(errs, "Choose a project folder (where all results will be saved).")
  else if (check_paths) {
    errs <- c(errs, .project_dir_errors(d))
    anc <- d$project_dir                             # missing folders are created, so look at the nearest existing one
    while (!dir.exists(anc) && dirname(anc) != anc) anc <- dirname(anc)
    if (!dir.exists(anc)) errs <- c(errs, "The project folder cannot be created there (is the drive mounted?).")
    else if (file.access(anc, 2) != 0) errs <- c(errs, sprintf("Cannot write to %s: choose another project folder.", anc))
    if (file.exists(file.path(d$project_dir, "run_summary.csv")))
      warns <- c(warns, "This folder already holds StRym results: the steps you run will replace their earlier results.")
  }
  if (!d$data_type %in% c("bulk", "single_cell")) errs <- c(errs, "Choose the type of RNA-seq data.")
  if (d$n_conditions < 1) errs <- c(errs, "The number of conditions must be at least 1.")
  if (d$n_timepoints < 1) errs <- c(errs, "The number of time points must be at least 1.")
  if (d$n_conditions > 20) warns <- c(warns, sprintf("%d conditions is a lot: please double-check.", d$n_conditions))
  if (d$n_conditions == 1) warns <- c(warns, "With a single condition there is nothing to compare: differential expression will be unavailable.")
  if (length(d$condition_names) && length(d$condition_names) != d$n_conditions)
    warns <- c(warns, sprintf("You listed %d condition names but set %d conditions.", length(d$condition_names), d$n_conditions))
  if (length(d$timepoint_names) && length(d$timepoint_names) != d$n_timepoints)
    warns <- c(warns, sprintf("You listed %d time-point names but set %d time points.", length(d$timepoint_names), d$n_timepoints))
  if (check_paths) {
    if (!nzchar(d$data_path)) errs <- c(errs, "Choose your expression data (counts table, single-cell object or reads folder).")
    else if (identical(d$data_kind, "fastq")) { r <- strym_validate_reads(d); errs <- c(errs, r$errors); warns <- c(warns, r$warnings); info <- c(info, r$info) }
    else { p <- strym_peek_data(d$data_path); if (isTRUE(attr(p, "bad"))) errs <- c(errs, paste("Data:", p)) else info <- c(info, paste("Data:", p)) }
    if (nzchar(d$sample_sheet)) {
      if (!file.exists(d$sample_sheet)) errs <- c(errs, "Sample sheet not found.")
      else {
        ss <- strym_read_sample_sheet(d$sample_sheet)
        if (inherits(ss, "error")) errs <- c(errs, paste("Sample sheet:", conditionMessage(ss)))
        else if (!nrow(ss)) errs <- c(errs, "The sample sheet has a header but no samples: add one row per sample.")
        else {
          cc <- .find_col(ss, "^condition"); tc <- .find_col(ss, "^time")
          if (is.na(cc)) errs <- c(errs, "The sample sheet needs a column named “condition”.")
          else { nc <- length(unique(ss[[cc]]))
                 info <- c(info, sprintf("Sample sheet: %d samples, %d conditions.", nrow(ss), nc))
                 if (nc != d$n_conditions) warns <- c(warns, sprintf("The sample sheet has %d conditions but you set %d.", nc, d$n_conditions))
                 typed <- setdiff(tolower(d$condition_names), tolower(trimws(as.character(ss[[cc]]))))                  # a typo in a condition name (draught for drought) must not pass unnoticed
                 if (length(typed)) warns <- c(warns, sprintf("You typed the condition name%s %s, but the sample sheet has: %s.", if (length(typed) == 1) "" else "s",
                                                              paste(d$condition_names[tolower(d$condition_names) %in% typed], collapse = ", "), paste(utils::head(unique(ss[[cc]]), 8), collapse = ", "))) }
          sc <- .find_col(ss, "^sample")
          if (is.na(sc) && identical(d$data_kind, "counts")) errs <- c(errs, sprintf("The sample sheet needs a column named 'sample' holding the column names of the table (it has: %s).", paste(utils::head(names(ss), 6), collapse = ", ")))
          if (!is.na(sc) && identical(d$data_kind, "counts")) {
            cs <- strym_counts_samples(d$data_path); nmv <- trimws(as.character(ss[[sc]]))
            if (length(cs) && !any(tolower(nmv) %in% tolower(cs)))
              errs <- c(errs, sprintf("None of the %d sample names in the sample sheet is a column of the counts table (the sheet has %s ...; the table has %s ...). The names must match.", nrow(ss), nmv[1], cs[1]))
            if (length(cs) && !all(nmv %in% cs))
              warns <- c(warns, sprintf("%d of the %d sample names in the sample sheet are not columns of the counts table (e.g. %s).",
                                        sum(!nmv %in% cs), nrow(ss), nmv[!nmv %in% cs][1]))
            back <- cs[!tolower(cs) %in% tolower(nmv)]
            if (length(cs) && length(back) && length(back) < length(cs))
              warns <- c(warns, sprintf("%d column(s) of the counts table are not in the sample sheet and will be left out of the analysis (e.g. %s).", length(back), back[1]))
            cd <- .find_col(ss, "^condition"); if (!is.na(cd) && any(is.na(ss[[cd]]) | !nzchar(trimws(as.character(ss[[cd]])))))
              errs <- c(errs, "The sample sheet has samples without a condition: fill in every row.")
          }
          if (!is.na(tc)) { nt <- length(unique(ss[[tc]]))
                 info <- c(info, sprintf("Sample sheet: %d time point%s.", nt, if (nt == 1) "" else "s"))
                 if (nt != d$n_timepoints) warns <- c(warns, sprintf("The sample sheet has %d time points but you set %d.", nt, d$n_timepoints))
                 typed <- setdiff(tolower(d$timepoint_names), tolower(trimws(as.character(ss[[tc]]))))
                 if (length(typed) && d$n_timepoints > 1) warns <- c(warns, sprintf("You typed the time-point name%s %s, but the sample sheet has: %s.", if (length(typed) == 1) "" else "s",
                                                              paste(d$timepoint_names[tolower(d$timepoint_names) %in% typed], collapse = ", "), paste(utils::head(unique(ss[[tc]]), 8), collapse = ", "))) }
          else if (d$n_timepoints > 1) warns <- c(warns, "The sample sheet has no “time_point” column although you have several time points.")
        }
      }
    } else if (d$n_conditions > 1 || d$n_timepoints > 1)
      warns <- c(warns, if (identical(d$data_kind, "fastq"))
        "No sample sheet: conditions are guessed from the file names (groups listed above). A sample sheet with “file”, “condition” and “time_point” columns makes this explicit: use the template below."
        else "No sample sheet: the analysis cannot tell which sample belongs to which condition / time point. Use the template below.")
    if (length(d$diff_path)) {                                                     # Cuffdiff results made elsewhere, next to the table
      miss <- d$diff_path[!file.exists(d$diff_path)]
      if (!identical(d$data_kind, "counts")) warns <- c(warns, "Differential-expression results are only used with a counts / expression table: they are ignored for this input.")
      else if (length(miss)) errs <- c(errs, paste("Differential-expression results not found:", miss[1]))
      else {
        n <- tryCatch(strym_count_diff_files(d$diff_path), error = function(e) 0L)
        if (n < 1) errs <- c(errs, "No Cuffdiff comparison files (gene_exp.diff) were found in the differential-expression results you gave.")
        else info <- c(info, sprintf("Differential-expression results: %d comparison file%s (q-values for the entropy steps).", n, if (n == 1) "" else "s"))
      }
    }
  }
  list(errors = unique(errs), warnings = unique(warns), info = info)
}

# ---- raw reads: samples, checks, counts --------------------------------------------------------------------------------
# The sample table the Tuxedo pipeline needs (group, file, mate) for the reads folder in the design.
# Uses the sample sheet if there is one (its "file" column names the FASTQ files, "condition" gives the group,
# and "time_point" is appended to the group when there are several time points); otherwise groups are guessed
# from the file names (as on the TUX tab). Extra columns: condition, time_point, replicate.
strym_quantify_samples <- function(design) {
  d <- strym_normalize_design(design); paired <- identical(d$read_type, "paired")
  if (nzchar(d$sample_sheet) && file.exists(d$sample_sheet)) {
    ss <- strym_read_sample_sheet(d$sample_sheet)
    if (inherits(ss, "error")) stop(conditionMessage(ss), call. = FALSE)
    fc <- .find_col(ss, "^file"); cc <- .find_col(ss, "^condition"); tc <- .find_col(ss, "^time"); rc <- .find_col(ss, "^replicate"); mc <- .find_col(ss, "^mate")
    sc <- .find_col(ss, "^sample")
    if (is.na(cc)) stop("The sample sheet needs a “condition” column.", call. = FALSE)
    if (is.na(fc) && is.na(sc)) stop("The sample sheet needs a “file” column naming each FASTQ file (or a “sample” column named after the files).", call. = FALSE)
    n <- nrow(ss)
    f <- if (!is.na(fc)) trimws(as.character(ss[[fc]])) else rep("", n); f[is.na(f)] <- ""
    sm <- if (!is.na(sc)) trimws(as.character(ss[[sc]])) else rep("", n); sm[is.na(sm)] <- ""
    avail <- sort(list.files(d$data_path, pattern = .FQ_EXT, full.names = TRUE))      # the reads actually in the folder
    if (length(d$fastq_files)) avail <- avail[basename(avail) %in% basename(d$fastq_files)]   # ... limited to the files chosen
    if (paired && length(avail)) { mi <- .mate_info(avail); avail <- avail[mi$first] }        # a sample name matches the first mate only (not R1 and R2)
    resolve <- function(file, sample) {
      if (nzchar(file)) {                                    # a file name was given (maybe without .gz, as Explorer shows it)
        for (cand in c(file, file.path(d$data_path, file), file.path(d$data_path, paste0(file, ".gz")),         # (a path relative to the data folder keeps same-named files in condition folders apart)
                       file.path(d$data_path, basename(file)), file.path(d$data_path, paste0(basename(file), ".gz"))))
          if (file.exists(cand) && !dir.exists(cand)) return(cand)
        hit <- avail[.stem(avail) == .stem(file)]
        return(if (length(hit) == 1) hit else file.path(d$data_path, basename(file)))  # not found: keep the name so it is reported
      }
      if (nzchar(sample)) {                                  # no file name: match the sample name to a file in the folder
        hit <- avail[.stem(avail) == sample]
        if (length(hit) != 1) hit <- avail[startsWith(.stem(avail), sample)]
        if (length(hit) == 1) return(hit)
      }
      NA_character_
    }
    full <- mapply(resolve, f, sm, USE.NAMES = FALSE)
    if (anyNA(full)) stop(sprintf(paste0("Sample sheet: no FASTQ file found for %s. Fill in the “file” column with each FASTQ name ",
                                         "(e.g. %s), or name each sample after its file (without .fastq.gz)."),
                                  paste(ifelse(nzchar(sm[is.na(full)]), sm[is.na(full)], paste("row", which(is.na(full)))), collapse = ", "),
                                  if (length(avail)) basename(avail[1]) else "sample1.fastq.gz"), call. = FALSE)
    cond <- trimws(as.character(ss[[cc]])); tp <- if (!is.na(tc)) trimws(as.character(ss[[tc]])) else rep("", nrow(ss))
    bad_c <- is.na(cond) | !nzchar(cond)
    if (any(bad_c)) stop(sprintf("Sample sheet: row%s %s %s no condition. Fill in every row (a row without a condition would be left out silently).", if (sum(bad_c) == 1) "" else "s", paste(which(bad_c), collapse = ", "), if (sum(bad_c) == 1) "has" else "have"), call. = FALSE)
    grp <- if (d$n_timepoints > 1 && !is.na(tc)) paste(cond, tp, sep = "_") else cond
    find_mate <- function(m) { m <- trimws(as.character(m)); for (cand in c(m, file.path(d$data_path, m), file.path(d$data_path, basename(m)))) if (!is.na(m) && nzchar(m) && file.exists(cand) && !dir.exists(cand)) return(cand); file.path(d$data_path, basename(m)) }
    mate <- if (!is.na(mc)) vapply(ss[[mc]], find_mate, "", USE.NAMES = FALSE)             # a mate given with its folder is found there
            else if (paired) .mate_info(full)$mate else rep(NA_character_, nrow(ss))
    out <- data.frame(group = grp, file = full, mate = mate, condition = cond, time_point = tp,
                      replicate = if (!is.na(rc)) ss[[rc]] else NA, stringsAsFactors = FALSE)
  } else {
    s <- strym_filter_fastq(infer_samples(d$data_path, paired), d$fastq_files, paired)
    out <- data.frame(group = s$group, file = s$file, mate = s$mate, condition = s$group, time_point = "", replicate = NA,
                      stringsAsFactors = FALSE)
  }
  out
}

# Checks for the "raw reads" case. list(errors, warnings, info). Needs tuxedo.R (infer_samples, estimate_space ...).
strym_validate_reads <- function(d) {
  errs <- character(); warns <- character(); info <- character()
  if (!dir.exists(d$data_path)) return(list(errors = "For raw reads, choose your FASTQ files (select several at once) or the FOLDER that holds them.", warnings = warns, info = info))
  gone <- d$fastq_files[!file.exists(d$fastq_files)]
  if (length(gone)) return(list(errors = sprintf("Chosen FASTQ file%s not found: %s (is the drive still connected?)", if (length(gone) == 1) "" else "s",
                                                 paste(basename(gone), collapse = ", ")), warnings = warns, info = info))
  for (x in list(list(d$gtf_file, "annotation (GTF/GFF) file"), list(d$fasta_file, "genome FASTA file")))
    if (!nzchar(x[[1]])) errs <- c(errs, paste("Choose the", x[[2]], "used to quantify the reads."))
    else if (!file.exists(x[[1]])) errs <- c(errs, paste("The", x[[2]], "was not found."))
  if (nzchar(d$fasta_file) && grepl("\\.(gz|zip)$", d$fasta_file, ignore.case = TRUE)) errs <- c(errs, "The genome FASTA must be uncompressed (unzip it first).")
  else if (nzchar(d$fasta_file) && file.exists(d$fasta_file)) { p <- genome_problem(d$fasta_file); if (nzchar(p)) errs <- c(errs, paste("Genome FASTA:", p)) }
  if (nzchar(d$gtf_file) && file.exists(d$gtf_file)) { p <- annotation_problem(d$gtf_file); if (nzchar(p)) errs <- c(errs, paste("Annotation:", p)) }
  long <- identical(d$read_platform, "long")
  if (long && nzchar(d$gtf_file) && file.exists(d$gtf_file)) {                     # featureCounts counts per gene_id: a GFF3 (ID= / Parent=) would fail after the alignment
    head_l <- readLines(d$gtf_file, n = 400, warn = FALSE); head_l <- head_l[!grepl("^#", head_l) & nzchar(head_l)]
    if (length(head_l) && !any(grepl("gene_id", head_l, fixed = TRUE)))
      errs <- c(errs, "For long reads the annotation must be a GTF with gene_id attributes (this one has none, so it looks like a GFF3). Convert it, for example with `gffread annotation.gff3 -T -o annotation.gtf`.")
  }
  if (!long && nzchar(d$adapter_file)) { if (!file.exists(d$adapter_file)) errs <- c(errs, "The adapter file was not found.")
    else if (!grepl("^>", readLines(d$adapter_file, n = 1, warn = FALSE)[1] %||% "")) errs <- c(errs, "The adapter file must be a FASTA file (it should start with a line beginning with >).") }
  tl <- if (long) check_long_tools() else check_tools(); if (!all(tl$found)) errs <- c(errs, paste("Missing tools:", paste(tl$tool[!tl$found], collapse = ", "), if (long) "(see \u201cLong reads\u201d in the README: setup_modules.sh installs them)." else "(see the TUX tab, Settings)."))
  s <- tryCatch(strym_quantify_samples(d), error = function(e) e)
  if (inherits(s, "error")) return(list(errors = c(errs, conditionMessage(s)), warnings = warns, info = info))
  if (!nrow(s)) return(list(errors = c(errs, "No FASTQ files (.fastq/.fq, optionally .gz) found in that folder."), warnings = warns, info = info))
  paired <- identical(d$read_type, "paired")
  ns <- normalize_samples(s[, c("group", "file", "mate")], paired)
  if (long && paired) errs <- c(errs, "Long reads are single reads: choose Single-end as the read type.")
  v <- validate_samples(ns, paired, platform = d$read_platform); errs <- c(errs, v$errors); warns <- c(warns, v$warnings)
  info <- c(info, sprintf("Reads: %d %s FASTQ file%s in %d group%s (%s).", nrow(ns), d$read_type, if (nrow(ns) == 1) "" else "s",
                          length(unique(ns$group)), if (length(unique(ns$group)) == 1) "" else "s",
                          paste(sprintf("%s ×%d", unique(ns$group), as.integer(table(ns$group)[unique(ns$group)])), collapse = ", ")))
  if (length(unique(ns$group)) != d$n_conditions * (if (d$n_timepoints > 1 && nzchar(d$sample_sheet)) d$n_timepoints else 1))
    warns <- c(warns, sprintf("The reads fall into %d groups, but %d condition%s × %d time point%s were expected.", length(unique(ns$group)),
                              d$n_conditions, if (d$n_conditions == 1) "" else "s", d$n_timepoints, if (d$n_timepoints == 1) "" else "s"))
  if (nzchar(d$project_dir) && !length(v$errors)) {
    est <- estimate_space(ns, paired); anc <- d$project_dir
    while (!dir.exists(anc) && dirname(anc) != anc) anc <- dirname(anc)
    free <- free_space_bytes(anc)
    line <- sprintf("Reads ≈ %s uncompressed · space needed ≈ %s · free on the project drive: %s", fmt_bytes(est$total), fmt_bytes(est$need), fmt_bytes(free))
    if (!is.na(free) && free < est$need) errs <- c(errs, paste("Not enough free space on the project drive.", line)) else info <- c(info, line)
  }
  list(errors = errs, warnings = warns, info = info)
}

# Genes x samples counts table from a Tuxedo run (Cuffdiff's per-replicate raw fragment counts).
# `tux_dir` is the output folder of run_tuxedo(). Gene names come from gene_exp.diff (ids if a name is missing/duplicated).
strym_counts_from_cuffdiff <- function(tux_dir) {
  used <- utils::read.csv(file.path(tux_dir, "samples_used.csv"), stringsAsFactors = FALSE)
  rg <- utils::read.delim(file.path(tux_dir, "diff_out", "genes.read_group_tracking"), stringsAsFactors = FALSE)
  de <- utils::read.delim(file.path(tux_dir, "diff_out", "gene_exp.diff"), stringsAsFactors = FALSE)
  key <- paste(used$group, used$k - 1)                                   # Cuffdiff numbers replicates 0, 1, ... within a group
  rg$sample <- used$sample[match(paste(rg$condition, rg$replicate), key)]
  if (anyNA(rg$sample)) stop("Could not match Cuffdiff replicates to the samples used.", call. = FALSE)
  wide <- stats::reshape(rg[, c("tracking_id", "sample", "raw_frags")], idvar = "tracking_id", timevar = "sample", direction = "wide")
  names(wide) <- sub("^raw_frags\\.", "", names(wide))
  wide <- wide[, c("tracking_id", used$sample), drop = FALSE]
  g <- de$gene[match(wide$tracking_id, de$test_id)]
  bad <- is.na(g) | g %in% c("-", "") | duplicated(g) | duplicated(g, fromLast = TRUE)
  data.frame(gene = ifelse(bad, wide$tracking_id, g), round(wide[, used$sample, drop = FALSE]), check.names = FALSE, row.names = NULL)
}

# Keep only the chosen FASTQ files (all of them when nothing was chosen). Mates count as chosen when either mate is. Renumbers replicates.
strym_filter_fastq <- function(s, chosen, paired) {
  if (!length(chosen) || !nrow(s)) return(s)
  keep <- basename(s$file) %in% basename(chosen) | (paired & !is.na(s$mate) & basename(s$mate) %in% basename(chosen))
  s <- s[keep, , drop = FALSE]
  if (!nrow(s)) return(s)
  normalize_samples(s[, c("group", "file", "mate")], paired)
}

# ---- the sample sheet StRym builds from the inputs --------------------------------------------------------------------
# Sample names in a counts table (the column names after the gene column), or NULL if it cannot be read.
.counts_samples_cache <- new.env()
strym_counts_samples <- function(path) {
  if (!nzchar(path) || !file.exists(path) || dir.exists(path)) return(NULL)
  info <- file.info(path); key <- paste(path, info$size, format(info$mtime, "%s"))        # the design page asks again and again: look at each file once
  if (!is.null(hit <- .counts_samples_cache[[key]])) return(if (identical(hit, NA)) NULL else hit)
  res <- .counts_samples_read(path); assign(key, if (is.null(res)) NA else res, envir = .counts_samples_cache); res
}
.counts_samples_read <- function(path) {
  suppressWarnings(tryCatch({
    ext <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
    nm <- if (ext %in% c("xlsx", "xls")) names(readxl::read_excel(path, n_max = 0, .name_repair = "minimal"))
          else { first <- readLines(path, n = 1, warn = FALSE)
                 strsplit(first, if (grepl("\t", first, fixed = TRUE)) "\t" else if (lengths(regmatches(first, gregexpr(";", first, fixed = TRUE))) > lengths(regmatches(first, gregexpr(",", first, fixed = TRUE)))) ";" else ",", fixed = TRUE)[[1]] }
    nm <- gsub('^"|"$', "", trimws(nm)); nm <- nm[-1]; nm <- nm[nzchar(nm)]
    tab <- if (file.size(path) < 30e6) tryCatch(read_table_any(path), error = function(e) NULL)    # small enough to look at the values too
    if (is.null(tab)) nm[!grepl(.ANNOT_NAMES, nm, ignore.case = TRUE)]                              # (otherwise the names have to do)
    else { k <- strym_count_columns(tab); nm[nm %in% trimws(names(tab)[k$keep])] }                 # text / annotation columns (gene_name, Chr, Start, End, Length ...) are not samples
  }, error = function(e) NULL))
}

# For each sample name: which of `tokens` appears in it as a separate word (NA when none).
# Time point written as the LAST part of every sample name (Pana_Control_T3, ctrl_5d, WT-day2, heat_0h_rep1): returns list(tp, cond) with the
# condition = the name without it (and without a trailing replicate such as _rep1), or NULL when the names do not all end that way
# or show fewer than 2 different time points. A replicate suffix (rep1 / r1) is never mistaken for a time point.
.infer_time_suffix <- function(smp) {
  rep_re <- "[_.-]?(rep|r)[0-9]+$"; time_re <- "^(.*?)[_.-]((t|tp|time|d|day|h|hr|hour|w|wk|week|min|m)?[0-9]+(\\.[0-9]+)?(h|hr|hrs|d|days?|w|wk|min|m|s)?)$"
  base <- sub(rep_re, "", smp, ignore.case = TRUE)
  # a bare replicate number after a time with a unit (Control_1h_1, Control_24h_2): the last number is the replicate, not the time
  unit_re <- "[_.-](t|tp|time|d|day|h|hr|hour|w|wk|week|min|m)?[0-9]+(\\.[0-9]+)?(h|hr|hrs|d|days?|w|wk|min|m|s)$|[_.-](t|tp|time|day|hour|week)[0-9]+(\\.[0-9]+)?$"
  b2 <- sub("[_.-][0-9]+$", "", base)
  if (all(grepl(unit_re, b2, ignore.case = TRUE)) && !all(grepl(unit_re, base, ignore.case = TRUE))) base <- b2
  m <- regmatches(base, regexec(time_re, base, ignore.case = TRUE))
  if (!length(m) || any(lengths(m) < 3)) return(NULL)
  cond <- vapply(m, `[`, "", 2); tp <- vapply(m, `[`, "", 3)
  if (any(!nzchar(cond)) || length(unique(tp)) < 2) return(NULL)
  list(tp = tp, cond = cond)
}

.find_token <- function(samples, tokens) {
  if (!length(tokens)) return(rep(NA_character_, length(samples)))
  esc <- gsub("([\\\\.|(){}^$*+?\\[\\]])", "\\\\\\1", tokens, perl = TRUE)
  vapply(samples, function(x) { m <- which(vapply(esc, function(t) grepl(paste0("(^|[^A-Za-z0-9])", t, "([^A-Za-z0-9]|$)"), x, ignore.case = TRUE), logical(1)))
                                if (length(m)) tokens[m[which.max(nchar(tokens[m]))]] else NA_character_ }, "", USE.NAMES = FALSE)   # the longest match: WT_Heat, not WT; 10, not 1
}

# Build the sheet (sample, condition, time_point, replicate[, file]) from what the researcher entered: the real FASTQ names, or the
# columns of the counts table. Conditions / time points come from the names the researcher typed (if they appear in the sample names),
# else from the groups guessed from the names. It is only a draft: the page shows it so the researcher can confirm or correct it.
# Returns list(sheet = data.frame or NULL, note = text).
strym_auto_sample_sheet <- function(design) {
  d <- strym_normalize_design(design); none <- function(msg) list(sheet = NULL, note = msg)
  if (!nzchar(d$data_path) || !file.exists(d$data_path)) return(none("Choose your input data above and the sample sheet appears here."))
  fastq <- identical(d$data_kind, "fastq")
  if (fastq) {
    s <- tryCatch(strym_quantify_samples(modifyList(d, list(sample_sheet = ""))), error = function(e) NULL)
    if (is.null(s) || !nrow(s)) return(none("No FASTQ files found yet: choose your FASTQ files or the folder that holds them."))
    smp <- .unique_sample_names(.stem(s$file), s$group); grp <- s$group           # same file name in two condition folders: the group is put in front
    files <- ifelse(startsWith(s$file, paste0(d$data_path, "/")), substring(s$file, nchar(d$data_path) + 2L), basename(s$file))   # relative to the data folder (a flat folder: just the name)
  } else if (identical(d$data_kind, "counts")) {
    smp <- strym_counts_samples(d$data_path)
    if (!length(smp)) return(none("The sample names could not be read from the counts table (is the first row a header?)."))
    grp <- trimws(sub("[_.-]?(rep|r)?[0-9]+$", "", smp, ignore.case = TRUE)); grp[!nzchar(grp)] <- smp[!nzchar(grp)]
    files <- NULL
  } else if (identical(d$data_kind, "10x") && dir.exists(d$data_path)) {       # a 10x folder: one sample per matrix folder, or per sample label inside the barcodes
    sm <- tryCatch(strym_10x_samples(d$data_path), error = function(e) NULL)
    if (is.null(sm) || !length(sm$samples)) return(none("Single cell/single nucleus data: no 10x matrix found yet. Choose the folder that holds matrix.mtx(.gz), barcodes.tsv(.gz) and features.tsv(.gz)."))
    smp <- sm$samples; files <- NULL
    grp <- trimws(sub("[_.-]?(rep|r)?[0-9]+$", "", smp, ignore.case = TRUE)); grp[!nzchar(grp)] <- smp[!nzchar(grp)]
    if (identical(sm$how, "barcodes")) {
      n <- length(smp); cn <- d$condition_names; tn <- d$timepoint_names
      cond <- .find_token(smp, cn)
      if (anyNA(cond)) cond <- if (length(cn) && length(unique(grp)) == length(cn)) cn[match(grp, unique(grp))] else grp
      tp <- if (d$n_timepoints <= 1) rep(if (length(tn)) tn[1] else "t1", n) else .find_token(smp, tn); tp[is.na(tp)] <- ""
      key <- paste(cond, tp); rep_n <- as.integer(ave(seq_len(n), key, FUN = seq_along))
      return(list(sheet = data.frame(sample = smp, condition = cond, time_point = tp, replicate = rep_n, stringsAsFactors = FALSE),
                  note = sprintf("%d samples were found inside the cell barcodes. StRym analyzes each one separately and then combines them: check that every sample has the right condition and correct any cell that is wrong.", n)))
    }
  } else return(none("Single cell/single nucleus data: StRym cannot list the samples inside this file. If you have several samples, supply a sample sheet below."))
  n <- length(smp); cn <- d$condition_names; tn <- d$timepoint_names
  cond <- .find_token(smp, cn)
  if (anyNA(cond)) cond <- if (length(cn) && length(unique(grp)) == length(cn)) cn[match(grp, unique(grp))] else grp
  tp <- if (d$n_timepoints <= 1) rep(if (length(tn)) tn[1] else "t1", n) else .find_token(smp, tn)
  tp[is.na(tp)] <- ""
  inferred <- FALSE
  if (!fastq && d$n_timepoints > 1 && !any(nzchar(tp))) {                         # no time-point names typed: look for one at the end of the sample names (..._T3, ..._5d, ..._day2)
    it <- .infer_time_suffix(smp)
    if (!is.null(it)) { tp <- it$tp; inferred <- TRUE; if (anyNA(.find_token(smp, cn))) cond <- it$cond }
  }
  key <- paste(cond, tp); seqn <- as.integer(ave(seq_len(n), key, FUN = seq_along))
  rep_n <- suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", smp, perl = TRUE)))
  rep_n <- if (!fastq && !inferred && !anyNA(rep_n) && !anyDuplicated(paste(key, rep_n))) rep_n else seqn
  out <- data.frame(sample = smp, condition = cond, time_point = tp, replicate = rep_n, stringsAsFactors = FALSE)
  if (fastq) out$file <- files
  note <- sprintf("%d samples found. Conditions and time points were worked out from the names%s: check them and correct any cell that is wrong.",
                  n, if (length(cn) || length(tn)) " you typed above" else "")
  list(sheet = out, note = note)
}

# ---- what a run will need: disk space and time (shown before it starts) --------------------------------------------------
# Calibrated on a real run on this computer: 4 single-end samples (6.1 GB gz = 79 million reads) took 3 h 17 min on 1 thread
# (alignment ~138 s per million reads, assembly ~165 s per sample, ~4 min for the index and the merge / Cuffdiff).
# Read counts are estimated from the file sizes (~77 bytes per read in a .gz, ~250 uncompressed). It is an estimate, not a promise.
strym_fmt_duration <- function(sec, up = FALSE) {                          # up = TRUE rounds the minutes up (for an upper bound)
  sec <- max(0, round(sec)); h <- sec %/% 3600; m <- (if (up) ceiling else round)((sec %% 3600) / 60)
  if (m == 60) { h <- h + 1; m <- 0 }
  if (h > 0) sprintf("%d h %02d min", h, m) else if (m >= 1) sprintf("%d min", m) else "under a minute"
}
strym_estimate_run <- function(design, steps) {
  d <- strym_normalize_design(design); out <- list(seconds = 0, reads_m = NA_real_, need = NA_real_, free = NA_real_, lines = character())
  if (nzchar(d$project_dir)) { anc <- d$project_dir; while (!dir.exists(anc) && dirname(anc) != anc) anc <- dirname(anc); out$free <- free_space_bytes(anc) }
  others <- setdiff(steps, "quantify")
  if ("quantify" %in% steps && identical(d$data_kind, "fastq")) {
    s <- tryCatch(strym_quantify_samples(d), error = function(e) NULL)
    if (!is.null(s) && nrow(s)) {
      paired <- identical(d$read_type, "paired"); f <- c(s$file, if (paired) s$mate); f <- f[!is.na(f) & file.exists(f)]
      long <- identical(d$read_platform, "long"); thr <- max(1, d$threads)
      bytes_per_read <- if (long) ifelse(grepl("\\.gz$", f), 2000, 4500) else ifelse(grepl("\\.gz$", f), 77, 250)   # a 1.5 kb read is ~2 kB compressed, ~4.5 kB as text
      reads_m <- sum(file.size(f) / bytes_per_read) / 1e6
      out$seconds <- if (long) 120 + reads_m * 600 / thr^0.85 + 60                                          # minimap2 + featureCounts: very roughly 10 min per million reads on one thread
                     else 120 + reads_m * 138 / thr^0.85 + nrow(s) * 165 / thr^0.5 + 120
      out$reads_m <- reads_m
      ns <- normalize_samples(s[, c("group", "file", "mate")], paired); out$need <- estimate_space(ns, paired)$need
      out$lines <- c(out$lines, sprintf(if (long) "Aligning about %.2f million long reads from %d sample%s with %d thread%s (minimap2), then counting them per gene (featureCounts)." else "Reading, cleaning and aligning about %.0f million reads from %d sample%s with %d thread%s, then Cufflinks, Cuffmerge and Cuffdiff.",
                                        reads_m, nrow(s), if (nrow(s) == 1) "" else "s", d$threads, if (d$threads == 1) "" else "s"))
    }
  }
  out$seconds <- out$seconds + 4 * length(others) + 40 * ("seurat_clustering" %in% others) + 20 * ("genist_network" %in% others)   # the other steps take seconds on a typical table (single cell and networks longer)
  out
}

# ---- sample sheet template ------------------------------------------------------------------------------------------
strym_sample_sheet_template <- function(design, replicates = 2L) {
  d <- strym_normalize_design(design)
  if (identical(d$data_kind, "fastq") && nzchar(d$data_path) && dir.exists(d$data_path)) {      # raw reads: list the real files
    s <- tryCatch(infer_samples(d$data_path, identical(d$read_type, "paired")), error = function(e) NULL)
    if (!is.null(s) && nrow(s)) {
      tn <- if (length(d$timepoint_names) >= 1) d$timepoint_names[1] else "t1"
      return(data.frame(sample = s$sample, condition = s$group, time_point = tn, replicate = s$k,
                        file = ifelse(startsWith(s$file, paste0(d$data_path, "/")), substring(s$file, nchar(d$data_path) + 2L), basename(s$file)),
                        stringsAsFactors = FALSE, row.names = NULL))                              # condition = guess from the names: edit it
    }
  }
  cn <- if (length(d$condition_names) == d$n_conditions) d$condition_names else paste0("condition", seq_len(max(1, d$n_conditions)))
  tn <- if (length(d$timepoint_names) == d$n_timepoints) d$timepoint_names else paste0("t", seq_len(max(1, d$n_timepoints)))
  g <- expand.grid(replicate = seq_len(replicates), time_point = tn, condition = cn, stringsAsFactors = FALSE)
  g$sample <- sprintf("%s_%s_rep%d", g$condition, g$time_point, g$replicate)
  g <- g[order(match(g$condition, cn), match(g$time_point, tn), g$replicate), ]
  data.frame(sample = g$sample, condition = g$condition, time_point = g$time_point, replicate = g$replicate,
             file = "", stringsAsFactors = FALSE, row.names = NULL)
}

# ---- running ---------------------------------------------------------------------------------------------------------------
.html <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", as.character(x), fixed = TRUE), fixed = TRUE), fixed = TRUE)

# The StRym logo as an <img> with the picture embedded in it (so a results page still shows it when it is moved or e-mailed).
# "" when there is no app folder or the logo file is missing, so the pages simply come without it.
strym_logo_img <- function(app_dir = NULL, height = 46) {
  if (is.null(app_dir) || length(app_dir) != 1 || !nzchar(app_dir)) return("")
  f <- file.path(app_dir, "www", "strym-logo-white.png")
  if (!file.exists(f)) return("")
  sprintf('<img alt="StRym" style="height:%dpx;width:auto;display:block;margin:0 0 10px" src="data:image/png;base64,%s">',
          as.integer(height), gsub("\n", "", jsonlite::base64_enc(readBin(f, "raw", file.size(f)))))
}

# One-page summary at the top of the project folder.
strym_write_index <- function(project_dir, design, summary, app_dir = NULL) {
  d <- strym_normalize_design(design)
  badge <- function(s) sprintf('<span class="b %s">%s</span>', s, s)
  rows <- vapply(seq_len(nrow(summary)), function(i) {
    r <- summary[i, ]
    files <- if (nzchar(r$files)) paste(sprintf('<a href="%s/%s">%s</a>', .html(r$folder), .html(strsplit(r$files, "; ")[[1]]), .html(strsplit(r$files, "; ")[[1]])), collapse = " · ") else ""
    sprintf('<tr><td><b>%s</b><br><a href="%s/">%s/</a></td><td>%s</td><td>%.1f s</td><td>%s<br><span class="m">%s</span></td></tr>',
            .html(r$label), .html(r$folder), .html(r$folder), badge(r$status), r$seconds, files, .html(r$message))
  }, character(1))
  facts <- c(Project = d$project_name, `Data type` = if (d$data_type == "bulk") "Whole organ/tissue RNA-seq" else "Single cell/single nucleus RNA-seq",
             Conditions = paste0(d$n_conditions, if (length(d$condition_names)) paste0(" (", paste(d$condition_names, collapse = ", "), ")") else ""),
             `Time points` = paste0(d$n_timepoints, if (length(d$timepoint_names)) paste0(" (", paste(d$timepoint_names, collapse = ", "), ")") else ""),
             Data = d$data_path, `Sample sheet` = if (nzchar(d$sample_sheet)) d$sample_sheet else "(none)",
             Created = format(Sys.time(), "%Y-%m-%d %H:%M"))
  # a quick look: what to look at first, and the key pictures
  look <- character()
  nf <- file.path(project_dir, "00_data_overview", "data_overview_notes.csv")
  if (file.exists(nf)) { nt <- tryCatch(utils::read.csv(nf, stringsAsFactors = FALSE), error = function(e) NULL)
                         if (!is.null(nt) && nrow(nt)) look <- c(look, nt$issue[nt$level == "warning"]) }
  bad <- summary[summary$status %in% c("failed", "skipped"), , drop = FALSE]
  if (nrow(bad)) look <- c(look, sprintf("%s (%s): %s", bad$label, bad$status, bad$message))
  glance <- if (length(look)) paste0('<h2>Things to look at</h2><div class="look"><ul>', paste0("<li>", .html(look), "</li>", collapse = ""), "</ul></div>")
            else '<div class="allgood">&#10004; Nothing needs attention: no problems were noticed in the data or the run.</div>'
  key <- c("00_data_overview/design_grid.png", "00_data_overview/sample_correlation.png", "02_differential_expression/deseq2/deseq2_volcano.png",
           "04_pca/pca_plot.png", "02_differential_expression/shannon_entropy/shannon_calls_per_timepoint.png", "05_trajectories/trajectories.png",
           "03_clustering/seurat_umap_tsne/umap_by_condition.png", "03_clustering/seurat_umap_tsne/tsne_by_condition.png")
  key <- key[file.exists(file.path(project_dir, key))]
  pics <- if (length(key)) paste0('<h2>At a glance</h2><div class="pics">', paste0(sprintf('<a href="%s"><img src="%s" alt="%s"></a>', .html(key), .html(key), .html(basename(key))), collapse = ""), "</div>") else ""
  html <- c('<!doctype html><html><head><meta charset="utf-8"><title>StRym results</title><style>',
    '.look{background:#FFF4DC;border-left:5px solid #E0A100;padding:4px 18px;margin:8px 0}.allgood{background:#E9F3D4;border-left:5px solid #586800;padding:10px 16px;margin:14px 0}',
    '.pics{display:flex;flex-wrap:wrap;gap:12px}.pics img{width:250px;max-height:380px;object-fit:contain;border:1px solid #D5DCEF;background:#fff;border-radius:4px}',

    'body{font-family:Arial,sans-serif;background:#F4F6FC;color:#1F2433;margin:0}',
    '.top{background:linear-gradient(120deg,#2F3A03,#6F7D1C);color:#fff;padding:26px 36px}.top h1{margin:0;font-size:26px}',
    '.top p{margin:6px 0 0;color:#E4E9BC}.wrap{padding:22px 36px;max-width:1100px}',
    'h2{color:#24347B;border-bottom:2px solid #BFCC46;padding-bottom:4px}table{border-collapse:collapse;width:100%;background:#fff}',
    'td,th{border:1px solid #E4EAF7;padding:8px 10px;text-align:left;vertical-align:top;font-size:14px}th{background:#E4EAF7;color:#24347B}',
    '.b{padding:2px 9px;border-radius:10px;font-size:12px;color:#fff}.ok{background:#586800}.placeholder{background:#427E93}',
    '.failed{background:#CC0000}.skipped{background:#808080}.m{color:#4D4D4D;font-size:12px}a{color:#344891}',
    '.foot{margin-top:28px;padding:14px 36px;background:#24347B;color:#E4EAF7;font-size:12px;text-align:center}</style></head><body>',
    sprintf('<div class="top">%s<h1>%s</h1><p>StRym — Streamlined transcriptomic Regulatory mapping</p></div><div class="wrap">', strym_logo_img(app_dir), .html(d$project_name)),
    '<h2>Your experiment</h2><table>', sprintf('<tr><th style="width:170px">%s</th><td>%s</td></tr>', .html(names(facts)), .html(facts)),
    '</table>', glance, pics, '<h2>Results</h2><table><tr><th>Step</th><th>Status</th><th>Time</th><th>Files</th></tr>', rows, '</table>',
    if (any(summary$status == "placeholder")) '<p class="m">Placeholder steps are stand-ins until the real analysis scripts are added (analysis/README.md).</p>' else '', '</div>',
    '<div class="foot">StRym was developed by Jordan Reed &amp; Stephanie M. Stewart. The platform is being maintained by the Sozzani lab at NC State University.</div></body></html>')
  writeLines(html, file.path(project_dir, "index.html"))
  file.path(project_dir, "index.html")
}

#' Run analysis steps for a design.
#' @param steps ids of the steps to run (default: the plan's defaults)
#' @return list(summary = data.frame, project_dir, index)
strym_run <- function(design, steps = NULL, app_dir = ".", progress = function(pct, label) {}, log = message) {
  d <- strym_normalize_design(design)
  v <- strym_validate_design(d, check_paths = FALSE); v$errors <- c(v$errors, .project_dir_errors(d))
  if (length(v$errors)) stop(paste(v$errors, collapse = " "), call. = FALSE)
  strym_load_registry(app_dir)
  reg <- strym_steps(); ids <- vapply(reg, `[[`, "", "id")
  if (is.null(steps)) { p <- strym_plan(d); steps <- p$id[p$default] }
  bad <- setdiff(steps, ids); if (length(bad)) stop("Unknown step(s): ", paste(bad, collapse = ", "), call. = FALSE)
  steps <- ids[ids %in% steps]                                         # registry order
  if (!length(steps)) stop("No analysis steps selected.", call. = FALSE)

  pd <- d$project_dir; dir.create(file.path(pd, "00_design"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(pd, "06_regulatory_networks"), showWarnings = FALSE)
  # A counts table without a sample sheet: the Experimental Design page drafts one from the column names (control_rep1, control_rep2 -> condition "control",
  # replicates 1 and 2) and uses it. A run started without the page (run_analysis_cli.R, the R functions) gets the same draft, so both give the same groups;
  # without it every column would be a condition of its own and DESeq2 would find "no replicates".
  if (identical(d$data_kind, "counts") && !nzchar(d$sample_sheet) && nzchar(d$data_path) && file.exists(d$data_path) && !dir.exists(d$data_path)) {
    a <- tryCatch(strym_auto_sample_sheet(d), error = function(e) NULL)
    if (!is.null(a) && is.data.frame(a$sheet) && nrow(a$sheet)) {
      f <- file.path(pd, "00_design", "sample_sheet_drafted_by_strym.csv"); utils::write.csv(a$sheet, f, row.names = FALSE); d$sample_sheet <- f
      log(sprintf("No sample sheet was given: StRym drafted one from the column names (%d samples, %d conditions). It is saved in %s.", nrow(a$sheet), length(unique(a$sheet$condition)), f))
    }
  }
  jsonlite::write_json(d, file.path(pd, "00_design", "design.json"), auto_unbox = TRUE, pretty = TRUE)
  plan <- strym_plan(d); plan$selected <- plan$id %in% steps
  utils::write.csv(plan[, c("id", "label", "group", "applicable", "default", "selected", "note")], file.path(pd, "00_design", "plan.csv"), row.names = FALSE)
  if (nzchar(d$sample_sheet) && file.exists(d$sample_sheet)) file.copy(d$sample_sheet, file.path(pd, "00_design", paste0("sample_sheet.", tools::file_ext(d$sample_sheet))), overwrite = TRUE)

  blocked <- strym_blocked(d, steps)
  done <- character(); results <- list()
  cur_data <- d$data_path; cur_sheet <- d$sample_sheet                 # what later steps receive (replaced by the quantification output)
  if (identical(d$data_kind, "fastq") && !"quantify" %in% steps && file.exists(strym_counts_path(d))) {
    cur_data <- strym_counts_path(d)                                   # quantification was run earlier: reuse its counts
    ss0 <- file.path(pd, "01_quantification", "sample_sheet.csv")
    if (!nzchar(cur_sheet) && file.exists(ss0)) cur_sheet <- ss0
  }
  wts <- ifelse(steps == "quantify", 10, 1); cum <- c(0, cumsum(wts)) / sum(wts)   # quantification is the long step
  for (k in seq_along(steps)) {
    s <- reg[[match(steps[k], ids)]]
    pct0 <- 5 + 90 * cum[k]; pct1 <- 5 + 90 * cum[k + 1]
    progress(pct0, sprintf("%s (%d of %d)", s$label, k, length(steps)))
    out_dir <- file.path(pd, s$dir); t0 <- Sys.time()
    res <- list(status = "skipped", files = character(), message = "")
    if (nzchar(blocked[[s$id]])) {
      res$message <- blocked[[s$id]]
    } else if (identical(d$data_kind, "fastq") && "quantify" %in% steps && s$id != "quantify" && !"quantify" %in% done && !file.exists(strym_counts_path(d))) {
      res$message <- "Skipped: the reads could not be quantified, so there are no counts to analyse."
    } else if (length(s$needs_any) && !identical(s$needs_any, "*") && !length(intersect(s$needs_any, done)) && !length(strym_existing_prereq(d, s)) &&
               !(is.function(s$needs_unless) && isTRUE(s$needs_unless(d)))) {
      res$message <- "Skipped: it uses the results of an earlier step of this run, and that step did not finish."
    } else {
      if (s$id != "quantify" && dir.exists(out_dir))                      # a re-run replaces the step's earlier results (never the long Tux run)
        unlink(list.files(out_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE), recursive = TRUE)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      pre <- if (identical(s$needs_any, "*")) done else intersect(s$needs_any, done)   # earlier steps this one builds on
      if (!length(pre) && !identical(s$needs_any, "*")) { pre <- strym_existing_prereq(d, s)                       # not in this run: use the results an earlier run left in the project folder
        if (length(pre)) log(sprintf("Using the results of %s from an earlier run in this project folder.", paste(pre, collapse = ", "))) }
      ctx <- list(design = d, data_path = cur_data, sample_sheet = cur_sheet, out_dir = out_dir, project_dir = pd, app_dir = app_dir,
                  prereq = setNames(lapply(pre, function(i) file.path(pd, reg[[match(i, ids)]]$dir)), pre),
                  log = log,
                  sub_progress = function(frac, label) progress(pct0 + (pct1 - pct0) * min(max(frac, 0), 1), sprintf("%s — %s", s$label, label)))
      env <- new.env(parent = globalenv())
      res <- tryCatch({
        source(file.path(app_dir, "analysis", s$script), local = env)
        r <- get(s$fn, envir = env)(ctx)
        # a step may hand a new data file / sample sheet on to the steps after it (quantification does)
        if (length(r$data_path) && nzchar(r$data_path)) cur_data <- r$data_path
        if (length(r$sample_sheet) && nzchar(r$sample_sheet) && !nzchar(d$sample_sheet)) cur_sheet <- r$sample_sheet
        list(status = r$status %||% "ok", files = as.character(r$files %||% character()), message = r$message %||% "")
      }, error = function(e) { log("ERROR in step ", s$label, ": ", conditionMessage(e)); list(status = "failed", files = character(), message = conditionMessage(e)) })
      if (res$status %in% c("ok", "placeholder")) done <- c(done, s$id)
    }
    results[[k]] <- data.frame(step = s$id, label = s$label, status = res$status,
                               seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), folder = s$dir,
                               files = paste(res$files, collapse = "; "), message = res$message, stringsAsFactors = FALSE)
    log(sprintf("[%s] %s: %s", res$status, s$label, res$message))
  }
  summary <- do.call(rbind, results)
  # keep the rows of earlier runs for steps that were not run this time (so a single-step run adds to the summary)
  all_rows <- summary; f <- file.path(pd, "run_summary.csv")
  old <- if (file.exists(f)) tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(old) && all(names(summary) %in% names(old))) {
    old <- old[!old$step %in% summary$step & old$step %in% ids, names(summary), drop = FALSE]
    old[is.na(old)] <- ""
    all_rows <- rbind(old, summary); all_rows <- all_rows[order(match(all_rows$step, ids)), , drop = FALSE]
  }
  utils::write.csv(all_rows, f, row.names = FALSE)
  index <- strym_write_index(pd, d, all_rows, app_dir)
  progress(100, "Done")
  invisible(list(summary = summary, all = all_rows, project_dir = pd, index = index))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- demo experiment ---------------------------------------------------------------------------------------------------------
# A small bulk RNA-seq experiment: 2 conditions (control / heat) x 3 time points x 2 replicates, 300 genes.
strym_make_demo_data <- function(dir, seed = 42) {
  set.seed(seed); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ss <- strym_sample_sheet_template(list(n_conditions = 2, condition_names = c("control", "heat"), n_timepoints = 3,
                                         timepoint_names = c("0h", "1h", "4h")), 2)
  ng <- 300; genes <- sprintf("GENE%03d", seq_len(ng)); base <- rgamma(ng, 2, 0.01)
  eff <- matrix(0, ng, nrow(ss))
  up <- sample(ng, 40); down <- sample(setdiff(seq_len(ng), up), 40)
  tfac <- c("0h" = 0, "1h" = .5, "4h" = 1)[ss$time_point]; heat <- ss$condition == "heat"
  eff[up, ] <- outer(rep(1, 40), ifelse(heat, 1.5 * tfac, 0)); eff[down, ] <- outer(rep(1, 40), ifelse(heat, -1.2 * tfac, 0))
  mu <- base * exp(eff)
  counts <- matrix(rnbinom(length(mu), mu = mu, size = 10), ng, nrow(ss), dimnames = list(genes, ss$sample))
  utils::write.csv(data.frame(gene = genes, counts, check.names = FALSE), file.path(dir, "demo_counts.csv"), row.names = FALSE)
  ss$file <- ss$sample
  utils::write.csv(ss, file.path(dir, "demo_sample_sheet.csv"), row.names = FALSE)
  invisible(list(counts = file.path(dir, "demo_counts.csv"), sample_sheet = file.path(dir, "demo_sample_sheet.csv")))
}
