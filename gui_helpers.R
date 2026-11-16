# gui_helpers.R -- reusable Shiny pieces for the StRym modules (GENIST, RTP-STAR).
# Requires drives.R (clean_path, unmounted_drive_problem) and jobs.R.

# ---- a path box with a Browse button (file or folder) ----------------------------------------------
# The example in an empty path box: Windows paths on Windows (StRym runs inside WSL there), the usual home-folder path elsewhere.
path_placeholder <- function(example = "file.xlsx") switch(strym_platform(),
  wsl = sprintf("Paste a path (D:\\data\\%s or /mnt/d/...) or browse", example), mac = sprintf("Paste a path (/Users/me/data/%s) or browse", example),
  windows = sprintf("Paste a path (D:/data/%s) or browse", example), sprintf("Paste a path (/home/me/data/%s) or browse", example))

pathUI <- function(id, label, type = c("file", "dir", "any"), help = NULL, multiple = FALSE,
                   placeholder = path_placeholder()) {
  type <- match.arg(type); ns <- NS(id)
  drop_hint <- switch(type, dir = "Drag and drop a folder here, or browse.",
                      file = if (multiple) "Drag and drop files here, or browse." else "Drag and drop a file here, or browse.",
                      any = if (multiple) "Drag and drop files or a folder here, or browse." else "Drag and drop a file or folder here, or browse.")
  div(class = "dz", `data-ns` = id,                                            # drop target (see www/strym.js)
    textInput(ns("txt"), label, "", width = "100%", placeholder = placeholder),
    if (type == "dir") actionButton(ns("btn"), "Browse…", class = "btn-default btn-sm")
    else if (type == "file") actionButton(ns("btn"), if (multiple) "Browse (select several)…" else "Browse…", class = "btn-default btn-sm")
    else tagList(actionButton(ns("btn"), if (multiple) "Browse for files (select several)…" else "Browse for a file…", class = "btn-default btn-sm"),   # "any": file(s) OR a folder
                 actionButton(ns("btn_dir"), "Browse for a folder…", class = "btn-default btn-sm")),
    span(class = "dz-hint", drop_hint),
    uiOutput(ns("note")),
    if (!is.null(help)) div(class = "hint", help),
    div(style = "display:none", fileInput(ns("up"), NULL, multiple = TRUE)),    # receives small dropped files that cannot be located on disk
    style = "margin-bottom:10px")
}

# Returns list(path = reactive (cleaned), set = function(value)).
UPLOAD_DIR <- file.path(path.expand("~/.tuxnet_r"), "uploads")      # small dropped files that could not be located on disk
.UPLOAD_MAX_ONE <- 300 * 1024^2; .UPLOAD_MAX_ALL <- 600 * 1024^2     # bytes; bigger files must be found on disk or browsed

# `multiple = TRUE` lets the user pick several files at once (in one folder); the folder goes into the box and the chosen files are
# returned by $files() (used for FASTQ reads). Files and folders can also be dropped on the field.
# Returns list(path = reactive (cleaned), files = reactive (chosen files, may be empty), set, set_files).
pathServer <- function(id, roots, type = c("file", "dir", "any"),
                       filetypes = c("xlsx", "xls", "csv", "txt", "tsv", "tab"), multiple = FALSE) {
  type <- match.arg(type)
  moduleServer(id, function(input, output, session) {
    br <- strym_browser_attach(input, output, session, key = "bf")        # the Browse window (filebrowser.R); `roots` is no longer needed and only kept so callers do not change
    sel <- reactiveVal(character()); sel_dir <- reactiveVal("")

    set_txt <- function(v) updateTextInput(session, "txt", value = as.character(v))
    # a choice made by browsing or dropping: one path (file or folder) or several files from one folder
    apply_paths <- function(p) {
      p <- as.character(p); p <- p[nzchar(p)]; if (!length(p)) return(invisible())
      if (type != "dir" && multiple && (length(p) > 1 || (type == "any" && grepl("\\.(fastq|fq)(\\.gz)?$", p, ignore.case = TRUE) && !dir.exists(p)))) {
        d <- unique(dirname(p))
        if (length(d) > 1) return(showNotification("Please choose files that are all in the same folder.", type = "warning", duration = 8))
        sel(p); sel_dir(d); set_txt(d)
      } else { sel(character()); set_txt(p[1]) }
    }
    observeEvent(input$txt, { if (length(sel()) && !identical(clean_path(input$txt), sel_dir())) sel(character()) }, ignoreInit = TRUE)  # typing a new path drops the file choice
    observeEvent(input$btn, br$open(if (type == "dir") "dir" else "file", if (type == "dir") "Choose a folder" else if (multiple) "Choose files" else "Choose a file",
                                    exts = filetypes, multiple = multiple && type != "dir", start = isolate(input$txt), cb = function(p) apply_paths(p)))
    observeEvent(input$btn_dir, br$open("dir", "Choose a folder", start = isolate(input$txt), cb = function(p) { sel(character()); set_txt(p[1]) }))

    # ---- drag and drop: find the dropped file(s) on this computer; if they cannot be found and are small, take a copy ----
    observeEvent(input$dropped, {
      d <- input$dropped; nm <- as.character(unlist(d[["names"]])); if (!length(nm)) return()
      sz <- suppressWarnings(as.numeric(unlist(d[["sizes"]]))); isd <- as.logical(unlist(d[["dirs"]]))
      if (length(sz) != length(nm)) sz <- rep(NA_real_, length(nm)); if (length(isd) != length(nm)) isd <- rep(FALSE, length(nm))
      warn <- function(m) showNotification(m, type = "warning", duration = 10)
      if (type == "dir" && !all(isd)) return(warn("Please drop a FOLDER here."))
      if (type == "file" && any(isd)) return(warn("Please drop a file here, not a folder."))
      if (!multiple && length(nm) > 1) {
        warn(sprintf("This field takes one item: using “%s”.", nm[1])); nm <- nm[1]; sz <- sz[1]; isd <- isd[1]
      }
      p0 <- clean_path(isolate(input$txt))
      hints <- c(if (nzchar(p0) && dir.exists(p0)) p0, if (nzchar(p0)) dirname(p0), unlist(isolate(sel_dir())))
      nid <- showNotification("Looking for your file(s) on this computer…", duration = NULL, closeButton = FALSE, type = "message")
      found <- tryCatch(locate_dropped(nm, sz, isd, hints), error = function(e) setNames(rep(NA_character_, length(nm)), nm))
      removeNotification(nid)
      if (!anyNA(found)) {
        apply_paths(unname(found))
        showNotification(sprintf("Found: %s", paste(to_windows_path(unname(found)), collapse = ", ")), duration = 8)
      } else if (!any(isd) && all(is.finite(sz)) && all(sz <= .UPLOAD_MAX_ONE) && sum(sz) <= .UPLOAD_MAX_ALL) {
        session$sendCustomMessage("strym-upload", list(ns = d[["ns"]], names = as.list(nm)))   # small: upload a copy
        showNotification("Could not find the original, so a copy is being uploaded…", duration = 6)
      } else warn(paste0("StRym could not find ", paste0("“", nm[is.na(found)], "”", collapse = ", "),
                         " on this computer. A browser does not reveal where a dropped file is stored, and large files are not copied. ",
                         "Please use the Browse button for these."))
    })
    observeEvent(input$up, {
      u <- input$up; if (is.null(u) || !nrow(u)) return()
      dest <- file.path(UPLOAD_DIR, format(Sys.time(), "%Y%m%d_%H%M%S")); dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      to <- file.path(dest, basename(u$name)); file.copy(u$datapath, to, overwrite = TRUE)
      apply_paths(to)
      showNotification(sprintf("Copied into %s", dest), duration = 8)
    })

    path <- reactive(clean_path(input$txt))
    output$note <- renderUI({                     # small inline hint when the path is not usable
      p <- path(); if (!nzchar(p)) return(NULL)
      u <- unmounted_drive_problem(p)
      if (nzchar(u)) div(class = "bad", u)
      else if (type == "file" && !file.exists(p)) div(class = "bad", "File not found.")
      else if (type == "any" && !file.exists(p)) div(class = "bad", "File or folder not found.")
      else if (type == "dir" && !dir.exists(p)) { anc <- p; while (!dir.exists(anc) && dirname(anc) != anc) anc <- dirname(anc)
                                                 if (identical(anc, "/") || !dir.exists(anc) || identical(anc, p)) div(class = "bad", "No part of this path exists (is the drive connected?).")
                                                 else div(class = "hint", "This folder does not exist yet: StRym will create it.") }
      else if (length(sel())) div(class = "hint", sprintf("%d file%s selected: %s", length(sel()), if (length(sel()) == 1) "" else "s",
                                                           paste(basename(sel()), collapse = ", ")))
    })
    list(path = path, files = reactive(sel()),
         set = function(v) { sel(character()); updateTextInput(session, "txt", value = v) },
         set_files = function(v) { v <- as.character(v); sel(v); sel_dir(if (length(v)) unique(dirname(v)) else "")
                                   if (length(v)) updateTextInput(session, "txt", value = unique(dirname(v))[1]) })        # the box shows the folder the files are in (it was left empty)
  })
}

# ---- cached table reads (validation re-runs often; xlsx reads are slow) ------------------------------
.read_cache <- new.env()
read_table_cached <- function(path, header = TRUE) {
  info <- file.info(path)
  key <- paste(path, info$size, format(info$mtime, "%s"), header)
  if (exists(key, envir = .read_cache, inherits = FALSE)) return(get(key, envir = .read_cache))
  val <- tryCatch(read_table_any(path, header = header), error = function(e) e)
  if (!inherits(val, "error")) assign(key, val, envir = .read_cache)                   # an error is not remembered: a file Excel had locked is read again once it is free
  val
}

# ---- the Run tab: progress, log, buttons --------------------------------------------------------------
job_panel_ui <- function(ns, run_label = "Run") {
  div(class = "card",
      uiOutput(ns("gate")),
      actionButton(ns("run"), run_label, class = "btn-primary"),
      actionButton(ns("cancel"), "Cancel", class = "btn-danger"),
      br(), br(),
      uiOutput(ns("progress")),
      p(class = "hint", "A run keeps going if you close this page; reopen it to see progress."),
      tags$pre(class = "log", id = ns("logbox"), textOutput(ns("log"))))
}

# Threads typed by the user -> a whole number between 1 and the number of cores (blank / text / decimals are repaired).
whole_threads <- function(x) {
  n <- suppressWarnings(as.integer(round(as.numeric(x)[1])))
  mx <- suppressWarnings(parallel::detectCores()); if (is.na(mx) || mx < 1) mx <- 1L
  if (length(n) != 1 || is.na(n) || n < 1) n <- 1L
  min(n, as.integer(mx))
}

# Wires a job to its panel. `script` is the CLI runner file name. Returns reactives.
job_controller <- function(name, script, input, output, session, can_run, app_dir, on_finished = function(status) {}) {
  ns <- session$ns
  rv <- reactiveValues(state = job_load(name))
  tick <- reactiveTimer(1000)
  running <- reactive({ tick(); job_alive(rv$state, script) })
  status <- reactive({                       # idle | running | ok | failed | cancelled
    tick()
    if (is.null(rv$state)) "idle" else if (job_alive(rv$state, script)) "running" else job_outcome(rv$state$log)
  })
  log_text <- reactive({ tick(); job_log_tail(rv$state$log) })

  output$gate <- renderUI({
    if (running()) div(class = "running", "Running…")
    else if (can_run()) div(class = "ok", "Ready to run.")
    else div(class = "bad", "Not ready — see the checks on the Inputs tab.")
  })
  observe({
    session$sendCustomMessage("tux-toggle", list(id = ns("run"), enable = can_run() && !running()))
    session$sendCustomMessage("tux-toggle", list(id = ns("cancel"), enable = running()))
  })
  output$log <- renderText(log_text())
  output$progress <- renderUI({
    st <- status(); pr <- job_progress(log_text())
    pct <- switch(st, ok = 100, idle = 0, max(pr$pct, 4))
    cls <- switch(st, ok = "progress-bar-success", failed = "progress-bar-danger",
                  cancelled = "progress-bar-warning", "progress-bar-striped active")
    label <- switch(st, idle = "Not started", ok = "Finished — see the Results tab",
                    failed = "Failed — see the log below", cancelled = "Canceled", pr$label)
    div(class = "progress", div(class = paste("progress-bar", cls), style = paste0("width:", pct, "%"), label))
  })

  prev <- reactiveVal("idle")
  observeEvent(status(), {
    st <- status()
    if (prev() == "running" && st %in% c("ok", "failed", "cancelled")) {
      if (!is.null(rv$state$out) && !is.na(rv$state$out) && dir.exists(rv$state$out))
        file.copy(rv$state$log, file.path(rv$state$out, paste0(name, "_run.log")), overwrite = TRUE)
      on_finished(st)
    }
    prev(st)
  })
  observeEvent(input$cancel, {
    if (job_alive(rv$state, script)) { kill_tree(rv$state$pid); showNotification("Run canceled.", type = "warning") }
  })
  list(state = reactive(rv$state), status = status, running = running, reset = function() { rv$state <- NULL },
       start = function(cfg, out, result) { rv$state <- tryCatch(job_start(name, script, cfg, app_dir, out, result),
                                                         error = function(e) { showNotification(conditionMessage(e), type = "error", duration = 12); rv$state }) })   # an error here would end the session (grey page)
}

# ---- network picture ---------------------------------------------------------------------------------------
# edges: data.frame(from, to, type = activates|inhibits|regulates, weight). extra_nodes: isolated nodes.
network_widget <- function(edges, extra_nodes = character(), height = "560px") {
  ids <- unique(c(edges$from, edges$to, extra_nodes))
  if (!length(ids)) return(NULL)
  outdeg <- as.numeric(table(factor(edges$from, levels = ids)))
  nodes <- data.frame(id = ids, label = ids, value = 8 + 4 * outdeg,
                      title = sprintf("%s — %d outgoing edge%s", ids, outdeg, ifelse(outdeg == 1, "", "s")),
                      stringsAsFactors = FALSE)
  cols <- c(activates = "#6F7D1C", inhibits = "#4156A1", regulates = "#808080")   # Genomic Green / Bio-Indigo / gray
  w <- edges$weight; w[!is.finite(w)] <- 1
  e <- data.frame(from = edges$from, to = edges$to, color = unname(cols[edges$type]),
                  width = 1 + 4 * w / max(w, 1e-12), arrows = "to",
                  title = sprintf("%s %s %s", edges$from, edges$type, edges$to), stringsAsFactors = FALSE)
  visNetwork::visNetwork(nodes, e, height = height, width = "100%") |>
    visNetwork::visIgraphLayout(layout = "layout_with_fr", randomSeed = 1) |>
    visNetwork::visNodes(shape = "dot", color = list(background = "#BFCC46", border = "#586800")) |>
    visNetwork::visEdges(smooth = FALSE) |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1), nodesIdSelection = TRUE) |>
    visNetwork::visInteraction(navigationButtons = TRUE)
}

network_legend <- function() {
  sw <- function(col, txt) span(style = "margin-right:16px", span(style = sprintf("display:inline-block;width:22px;height:4px;background:%s;vertical-align:middle;margin-right:5px", col)), txt)
  div(class = "hint", sw("#6F7D1C", "activates"), sw("#4156A1", "inhibits / represses"), sw("#808080", "direction of effect not determined"),
      "· larger node = more outgoing edges · thicker edge = stronger")
}
