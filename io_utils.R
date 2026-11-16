# io_utils.R -- shared input/output helpers for the GENIST and RTP-STAR modules.

# Read a table from .xlsx/.xls, .csv, or delimited .txt/.tsv. The first row is the header
# (as in MATLAB's readtable/dataset, which the original tools use). Column names are kept as-is.
read_table_any <- function(path, sheet = 1, header = TRUE) {
  if (!nzchar(path) || !file.exists(path)) stop("File not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
  df <- if (ext %in% c("xlsx", "xls", "xlsm")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required to read Excel files.")
    sheets <- readxl::excel_sheets(path)
    # the chosen sheet; when it is empty or has no numbers (a cover / notes sheet) the first sheet that does is used
    pick <- if (is.numeric(sheet)) sheet else match(sheet, sheets)
    has_num <- function(s) { x <- tryCatch(readxl::read_excel(path, sheet = s, col_names = FALSE, n_max = 30, .name_repair = "minimal"), error = function(e) NULL)
                             !is.null(x) && ncol(x) > 1 && any(vapply(x[-1], is.numeric, TRUE)) }
    if (!is.na(pick) && !has_num(pick)) { alt <- which(vapply(seq_along(sheets), has_num, TRUE)); if (length(alt)) { message("Note: sheet “", sheets[pick], "” has no numbers; using sheet “", sheets[alt[1]], "”."); pick <- alt[1] } }
    d <- as.data.frame(readxl::read_excel(path, sheet = pick, col_names = header, .name_repair = "minimal"), stringsAsFactors = FALSE)
    # a title or note above the real header (merged cells): blank header cells mean the real header is further down -- shift it up
    for (k in 1:8) {
      nm <- names(d); if (!header || !nrow(d) || mean(is.na(nm) | !nzchar(trimws(nm))) < 0.5) break           # a title row leaves most header cells blank; one blank cell is just an unnamed index column
      first <- as.character(unlist(d[1, ]))
      if (all(is.na(first) | !nzchar(trimws(first)))) { d <- d[-1, , drop = FALSE]; rownames(d) <- NULL; next }   # an empty row between the title and the real header
      if (sum(!is.na(first) & nzchar(first)) < max(2, ncol(d) / 2)) {              # another title / note line: skip it when a full row follows within a few lines
        nb <- vapply(seq_len(min(8, nrow(d))), function(i) sum(!is.na(d[i, ]) & nzchar(as.character(unlist(d[i, ])))), 0)
        if (!any(nb >= max(2, ncol(d) / 2))) break
        d <- d[-1, , drop = FALSE]; rownames(d) <- NULL; next }
      names(d) <- ifelse(is.na(first) | !nzchar(first), paste0("V", seq_along(first)), first); d <- d[-1, , drop = FALSE]; rownames(d) <- NULL
    }
    if (header && ncol(d) && (is.na(names(d)[1]) || !nzchar(trimws(names(d)[1])))) names(d)[1] <- "id"       # the unnamed index column of an export from R / pandas
    d
  } else {
    first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
    n_semi <- lengths(regmatches(first, gregexpr(";", first, fixed = TRUE))); n_comma <- lengths(regmatches(first, gregexpr(",", first, fixed = TRUE)))
    sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (n_semi > n_comma) ";" else if (n_comma > 0) "," else ""
    dec <- "."                                                       # European exports: semicolons and decimal commas (1,5)
    if (identical(sep, ";")) { body <- utils::head(readLines(path, n = 60, warn = FALSE, encoding = "UTF-8")[-1], 50); if (any(grepl("[0-9],[0-9]", body))) dec <- "," }
    d <- NULL
    if (file.size(path) > 5e6 && nzchar(sep) && requireNamespace("data.table", quietly = TRUE))             # big table: data.table's reader is ~50x faster (every step reads the table again)
      d <- tryCatch(data.table::fread(path, sep = sep, dec = dec, header = header, data.table = FALSE, check.names = FALSE, strip.white = TRUE,
                                      fill = TRUE, encoding = "UTF-8", quote = "\"", showProgress = FALSE), error = function(e) NULL, warning = function(w) NULL)
    if (is.null(d)) d <- utils::read.table(path, header = header, sep = sep, dec = dec, quote = "\"", comment.char = "", check.names = FALSE,
                                           stringsAsFactors = FALSE, fill = TRUE, strip.white = TRUE,
                                           fileEncoding = "UTF-8-BOM")
    # R's write.table(row.names = TRUE) writes a header ONE FIELD SHORTER than its rows, and read.table then turns the first column into row names:
    # the gene ids would silently disappear and every column shift. Put them back as the first column.
    rn <- attr(d, "row.names")
    if (header && is.character(rn)) { d <- cbind(data.frame(id = rn, stringsAsFactors = FALSE), d); rownames(d) <- NULL }
    d
  }
  if (!ncol(df)) stop("No columns found in ", basename(path), call. = FALSE)
  if (is.numeric(df[[1]])) df[[1]] <- formatC(df[[1]], format = "fg", digits = 15, flag = "")     # numeric ids: 100000, not 1e+05
  df[[1]] <- trimws(as.character(df[[1]]))
  df
}

# First column as a clean character vector (NA/empty removed, duplicates kept in order).
read_id_column <- function(path) {
  x <- read_table_any(path)[[1]]
  x[!is.na(x) & nzchar(x)]
}

# Numeric matrix from columns 2..end of a table (gene ids in column 1 become row names).
# Non-numeric cells become NA. Duplicate gene ids: first occurrence kept.
#
# A column that is not made of numbers (a gene symbol, a description, an annotation) is NOT data: it is left out, and its name is returned in
# attr(result, "dropped") so the caller can say so. (Otherwise it would count as an extra time point or sample made of missing values.)
table_to_matrix <- function(df, cols = seq_len(ncol(df))[-1]) {
  ids <- df[[1]]
  keep <- !is.na(ids) & nzchar(ids) & !duplicated(ids)
  is_num <- vapply(df[cols], function(v) {
    s <- trimws(as.character(v)); s <- s[!is.na(s) & nzchar(s)]
    length(s) > 0 && mean(!is.na(suppressWarnings(as.numeric(s)))) >= 0.5
  }, logical(1))
  dropped <- names(df)[cols][!is_num]; cols <- cols[is_num]
  m <- if (length(cols)) sapply(df[cols], function(v) suppressWarnings(as.numeric(v)), simplify = TRUE) else matrix(numeric(), nrow(df), 0)
  if (!is.matrix(m)) m <- matrix(m, nrow = nrow(df), dimnames = list(NULL, names(df)[cols]))
  rownames(m) <- ids
  m <- m[keep, , drop = FALSE]
  attr(m, "dropped") <- dropped
  m
}

# Say once (as a message) which non-numeric columns table_to_matrix() left out.
note_dropped_columns <- function(m, what) {
  d <- attr(m, "dropped")
  if (length(d)) message("Note: the ", what, " has non-numeric column(s) (", paste(head(d, 5), collapse = ", "), if (length(d) > 5) ", ..." else "", ") that are not data and were ignored.")
  invisible(m)
}

# Do two paths name the same file? (Windows drives, /mnt/<letter>, are not case sensitive.)
same_file <- function(a, b) {
  norm <- function(p) { p <- normalizePath(p, winslash = "/", mustWork = FALSE); if (grepl("^/mnt/[a-z](/|$)", p)) tolower(p) else p }
  length(a) == 1 && length(b) == 1 && !is.na(a) && !is.na(b) && nzchar(a) && nzchar(b) && identical(norm(a), norm(b))
}

# Problems when a file the program is about to WRITE is also one of the researcher's INPUT files (the rule: inputs are never changed).
# `inputs` = named character vector (label = path), `outputs` = the paths that will be written.
output_overwrites_input <- function(inputs, outputs) {
  p <- character()
  for (i in seq_along(inputs)) for (o in outputs)
    if (same_file(inputs[[i]], o)) { p <- c(p, sprintf("The results would overwrite your %s (%s). Choose another output folder or file name: StRym never changes your input files.", names(inputs)[i], basename(o))); break }
  p
}

# A name that is safe in a file name: letters, digits, . _ - only; runs of anything else become one "_".
strym_safe_name <- function(x, fallback = "project") {
  n <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(x)[1])))
  if (is.na(n) || !nzchar(n)) fallback else n
}

# Informative name of a network table for Cytoscape: <project>_<comparison>_<method>_network_for_Cytoscape.txt
#   strym_network_filename("Heat shock", "GENIST", "Heat_vs_Control")  ->  "Heat_shock_Heat_vs_Control_GENIST_network_for_Cytoscape.txt"
strym_network_filename <- function(project, method, what = "") {
  parts <- c(strym_safe_name(project, ""), if (nzchar(what)) strym_safe_name(what, ""), method, "network_for_Cytoscape")
  paste0(paste(parts[nzchar(parts)], collapse = "_"), ".txt")
}

# A quick picture of a network (a sanity check, not a figure for a paper): PNG next to the results, shown in the interactive browser.
# edges = data.frame(from, to, type = activates|inhibits|regulates, weight). Never stops the analysis: returns FALSE when it cannot draw.
strym_network_png <- function(edges, file, title = "", max_edges = 1500, size = 1500) {
  tryCatch({
    if (!requireNamespace("igraph", quietly = TRUE) || is.null(edges) || !nrow(edges)) return(FALSE)
    edges <- edges[!is.na(edges$from) & !is.na(edges$to), , drop = FALSE]
    if (!nrow(edges)) return(FALSE)
    w <- suppressWarnings(as.numeric(edges$weight)); w[!is.finite(w)] <- 1; edges$weight <- w
    total <- nrow(edges)
    if (total > max_edges) edges <- edges[order(-edges$weight, method = "radix")[seq_len(max_edges)], , drop = FALSE]    # a hairball of thousands of edges tells nothing: strongest only
    g <- igraph::graph_from_data_frame(data.frame(from = edges$from, to = edges$to), directed = TRUE)
    n <- igraph::vcount(g)
    set.seed(1)
    lay <- if (n <= 400) igraph::layout_with_fr(g, niter = 300) else if (n <= 3000) igraph::layout_with_drl(g) else igraph::layout_randomly(g)
    outdeg <- igraph::degree(g, mode = "out")
    vs <- 3 + 9 * sqrt(outdeg / max(1, max(outdeg)))
    cols <- c(activates = "#6F7D1C", inhibits = "#4156A1", regulates = "#9A9A9A")
    ec <- unname(cols[as.character(edges$type)]); ec[is.na(ec)] <- "#9A9A9A"
    lw <- 0.5 + 2.5 * edges$weight / max(edges$weight, 1e-12)
    lab_n <- if (n <= 60) n else 15                                                    # all names when few; only the main regulators when many
    lab <- rep(NA_character_, n); top <- order(-outdeg, method = "radix")[seq_len(min(lab_n, n))]; lab[top] <- igraph::V(g)$name[top]
    grDevices::png(file, width = size, height = size, res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(0.5, 0.5, 2.5, 0.5))
    plot(g, layout = lay, vertex.size = vs, vertex.color = "#BFCC46", vertex.frame.color = "#6F7D1C", vertex.label = lab, vertex.label.cex = if (n <= 60) 0.8 else 0.7,
         vertex.label.color = "#1F2433", vertex.label.dist = 0.4, edge.color = ec, edge.width = lw, edge.arrow.size = if (n <= 200) 0.4 else 0.2,
         edge.curved = if (n <= 60) 0.15 else 0, main = "")
    graphics::title(main = sprintf("%s%s", title, if (total > nrow(edges)) sprintf(" (strongest %d of %d edges shown)", nrow(edges), total) else ""), cex.main = 1)
    graphics::legend("bottomleft", legend = c("activates", "inhibits", "direction not determined"), col = unname(cols), lwd = 3, bty = "n", cex = 0.8)
    TRUE
  }, error = function(e) { message("Note: the network picture could not be drawn (", conditionMessage(e), ")."); FALSE })
}

# Case-sensitive, locale-independent ordering (matches MATLAB's sort on cellstr).
c_sort <- function(x) x[order(x, method = "radix")]
c_unique_sort <- function(x) c_sort(unique(x))

# Map gene ids to symbols using a 2-column table (id, symbol). Ids without a symbol keep their id.
# Multiple symbols for one id: the first is used.
map_symbols <- function(ids, symbol_df) {
  if (is.null(symbol_df) || !nrow(symbol_df) || ncol(symbol_df) < 2) return(ids)
  key <- as.character(symbol_df[[1]]); val <- as.character(symbol_df[[2]])
  ok <- !is.na(key) & !is.na(val) & nzchar(val)
  key <- key[ok]; val <- val[ok]
  pos <- match(ids, key)
  out <- ids
  hit <- !is.na(pos)
  out[hit] <- val[pos[hit]]
  out
}

# Progress / log markers understood by the GUI ("[PROGRESS 40] label").
emit_progress <- function(pct, label) if (!isTRUE(getOption("strym.quiet_progress"))) message(sprintf("[PROGRESS %d] %s", as.integer(round(pct)), label))   # quiet: a step that runs GENIST inside a bigger run keeps its own progress bar

# NA (like a MATLAB NaN comparison) counts as FALSE.
.isT <- function(x) !is.na(x) & x
