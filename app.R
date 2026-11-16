# StRym (Streamlined transcriptomic Regulatory mapping): the Shiny app.
#
# Pages (navigation bar): About, Experimental Design (mod_design.R), Analysis modules (TUX tab written in this file, GENIST in mod_genist.R, RTP-STAR in
# mod_rtpstar.R) and Drives (drives.R). Long jobs run as separate R processes (jobs.R) so a run keeps going when the page is closed.
# Order of this file: packages and helper sources, TUX run-state helpers, the Drives page text (it depends on the platform), the page layout (ui), the server.
# Launch:  start_gui.bat (Windows)  or  bash start_gui.sh (inside WSL/Linux; both start Rscript and serve http://localhost:3838).
suppressPackageStartupMessages({
  library(shiny); library(DT); library(processx); library(jsonlite)
})
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a          # base R from 4.4
options(shiny.maxRequestSize = 600 * 1024^2)                 # small dropped files that cannot be located are uploaded (never FASTQ-sized ones)
app_dir <- normalizePath(if (exists("app_dir_override")) app_dir_override else getwd())
for (f in c("tuxedo.R", "drives.R", "io_utils.R", "jobs.R", "genist.R", "rtpstar.R", "pipeline.R", "filebrowser.R", "gui_helpers.R",
            "about.R", "mod_genist.R", "mod_rtpstar.R", "mod_design.R")) source(file.path(app_dir, f))
strym_load_registry(app_dir)                                 # analysis steps + rules (analysis/registry.R)
invisible(try(drives_nowait(), silent = TRUE))                # ask Windows about the drives in the BACKGROUND (a few seconds): starting the app and opening a Browse window never wait for it
invisible(try({ strym_places(); strym_list_dir(path.expand("~")) }, silent = TRUE))        # first use of the pieces behind the Browse window now, so the first click is as quick as the rest

STATE_DIR  <- path.expand("~/.tuxnet_r")                    # lives on the Linux disk
STATE_FILE <- file.path(STATE_DIR, "run_state.json")
dir.create(STATE_DIR, showWarnings = FALSE)
default_threads <- max(1L, min(4L, parallel::detectCores() - 1L, na.rm = TRUE))
STEP_LABELS <- c("Cleaning and aligning reads", "Assembling transcripts",
                 "Merging assemblies", "Differential expression")

# ---- run-state helpers (a run keeps going if the browser tab is closed) ----
pid_alive <- function(pid) pid_alive_pattern(pid, "run_cli")       # works with /proc (Linux, WSL) and with ps (macOS); the pattern guards against pid reuse
load_state <- function() if (file.exists(STATE_FILE)) tryCatch(fromJSON(STATE_FILE), error = function(e) NULL)
run_outcome <- function(logfile) {                            # "ok" | "failed" | "cancelled"
  txt <- if (file.exists(logfile)) paste(readLines(logfile, warn = FALSE), collapse = "\n") else ""
  if (grepl("Done. Results in", txt, fixed = TRUE)) "ok"
  else if (grepl("ERROR:", txt, fixed = TRUE)) "failed" else "cancelled"
}

# ---- the Drives page: what is shown depends on the platform (see strym_platform() in drives.R) ----------------------------------------------
# Windows (StRym inside WSL): drives need to be mounted into Ubuntu, so there are Mount / Release buttons. Mac and Linux: the system mounts the
# drives; the page only lists them. Only the headings are bold; the longer instructions belong in the protocol.
platform <- strym_platform()
drives_intro <- switch(platform,
  wsl = "StRym can use every drive that Windows can see. A drive plugged in after StRym started has to be mounted once.",
  mac = "The drives your Mac has mounted are listed here (they live under /Volumes).",
  windows = "The drives that have a letter are listed here.",
  "The drives your system has mounted are listed here (usually under /media, /run/media or /mnt).")
drives_card <- div(class = "card", h4("Drives"), p(class = "hint", drives_intro),
  tableOutput("drives"),
  if (platform == "wsl") fluidRow(
    column(5, selectInput("mount_letter", "Drive", choices = "…", width = "100%")),
    column(7, div(style = "margin-top:25px",
                  actionButton("mount", "Mount drive"), actionButton("unmount", "Release drive"), actionButton("drives_refresh", "Refresh"))))
  else actionButton("drives_refresh", "Refresh"),
  uiOutput("mount_msg"))

drives_help <- if (platform == "wsl") { div(class = "card", h4("Using an external drive"),
  tags$ol(
    tags$li("Plug the drive in and wait until Windows File Explorer shows it."),
    tags$li("Press Refresh. If the status says not mounted or needs reconnecting, select the drive and press Mount drive."),
    tags$li(HTML("When the status says ready, type a path such as <code>D:\\data</code> in any box, or use Browse."))),
  h4("Before you unplug it"),
  tags$ol(
    tags$li("Make sure no StRym run is using the drive."),
    tags$li("Select it here and press Release drive, then eject it in Windows.")),
  h4("If it does not work"),
  tags$ul(
    tags$li("Check that File Explorer shows the drive and that it has a drive letter."),
    tags$li("Press Refresh and try Mount drive again."),
    tags$li(HTML("As a last resort close this page, run <code>wsl --shutdown</code> in PowerShell and start StRym again. A run in progress stops; saved results are kept."))),
  p(class = "hint", "The protocol document has the full details, with screenshots."))
} else { div(class = "card", h4("Using an external drive"),
  tags$ol(
    tags$li("Plug the drive in and wait until it appears in your file manager."),
    tags$li("Press Refresh. The drive is listed with its folder."),
    tags$li(HTML(sprintf("Type a path such as <code>%s</code> in any box, or use Browse.",
                         switch(platform, mac = "/Volumes/MyDrive/data", windows = "D:/data", "/media/me/MyDrive/data"))))),
  h4("Before you unplug it"),
  tags$ol(
    tags$li("Make sure no StRym run is using the drive."),
    tags$li(switch(platform, mac = "Eject it in Finder.", windows = "Eject it with Safely Remove Hardware.", "Unmount or eject it in your file manager."))),
  h4("If it does not appear"),
  tags$ul(
    tags$li("Check that your file manager shows the drive."),
    if (platform == "linux") tags$li(HTML("If it is not mounted, mount it in the file manager or with <code>mount</code>, then press Refresh.")),
    tags$li("Press Refresh again after reconnecting the drive.")),
  p(class = "hint", "The protocol document has the full details."))
}

ui <- navbarPage(
  title = tags$img(src = "strym-logo-white.png", alt = "StRym"), id = "module", windowTitle = "StRym",
  selected = "About", collapsible = TRUE,
  header = tags$head(
    tags$link(rel = "icon", type = "image/png", href = "favicon.png?v=9b"), tags$link(rel = "apple-touch-icon", href = "favicon-192.png?v=9b"),   # ?v=9b: browsers cache tab icons hard
    tags$link(rel = "stylesheet", href = "strym.css?v=11"), tags$script(src = "strym.js?v=11"),                  # ?v=11: a browser keeps old copies of these files; a new number makes it fetch the new ones
    tags$script(HTML("Shiny.addCustomMessageHandler('tux-toggle',function(m){
      var e=document.getElementById(m.id); if(e){e.disabled=!m.enable;}});
      setInterval(function(){document.querySelectorAll('pre.log').forEach(function(e){e.scrollTop=e.scrollHeight;});},800);"))),
  tabPanel("About", about_ui()),
  tabPanel("Experimental Design", br(), design_ui("design")),
  navbarMenu("Analysis modules",
  tabPanel("TUX · RNA-seq", tabsetPanel(id = "tabs",
    tabPanel("1 · Inputs", br(),
      fluidRow(
        column(6,
          div(class = "card", h4("Data"),
            textInput("dir_txt", "Folder with FASTQ files (read-only)", "", width = "100%",
                      placeholder = path_placeholder("")),
            actionButton("dir_in", "Browse…", class = "btn-default btn-sm"),
            p(class = "hint", "Either all FASTQ files in one folder (groups are guessed from file names ",
              "and can be edited), or one sub-folder per condition."),
            radioButtons("read_type", "Read type",
                         c("Single-end" = "single",
                           "Paired-end (mates named A_/B_, _R1/_R2 or _1/_2)" = "paired"))),
          div(class = "card", h4("Reference"),
            textInput("gtf_txt", "Annotation (GTF/GFF)", "", width = "100%"),
            actionButton("gtf", "Browse…", class = "btn-default btn-sm"),
            textInput("fasta_txt", "Genome (FASTA, uncompressed)", "", width = "100%"),
            actionButton("fasta", "Browse…", class = "btn-default btn-sm")),
          div(class = "card", h4("Output"),
            textInput("out_txt", "Output folder (can be on any drive with enough space)", "", width = "100%"),
            actionButton("out_in", "Browse…", class = "btn-default btn-sm"),
            div(class = "path", textOutput("out_win", inline = TRUE)))),
        column(6,
          div(class = "card", h4("Checks"), uiOutput("validation")),
          div(class = "card", h4("Samples"),
            p(class = "hint", "Double-click a Group to edit it. Leave Group blank to skip a sample."),
            DTOutput("samples"),
            actionButton("rescan", "Re-detect samples"))))),
    tabPanel("2 · Settings", br(),
      fluidRow(
        column(6, div(class = "card",
          h4("Read cleaning"),
          selectInput("cleaner", "Cleaning program", c("Automatic" = "auto", "fastq-mcf (TuxNet)" = "fastq-mcf", "fastp (optional)" = "fastp"), "auto"),
          numericInput("q", "Quality threshold", 30, 0, 41, 1),
          numericInput("l", "Minimum sequence length", 30, 1, 500, 1),
          numericInput("w", "Window size", 5, 1, 100, 1),
          p(class = "hint", "TuxNet defaults: 30 / 30 / 5. Automatic uses fastq-mcf when it is installed and fastp otherwise (macOS has no fastq-mcf). The two programs trim slightly differently, so results can differ a little between them."))),
        column(6, div(class = "card",
          h4("Compute"),
          numericInput("threads", "Threads", default_threads, min = 1, max = max(1L, parallel::detectCores()), step = 1),
          h4("Advanced"),
          checkboxInput("keep_cleaned", "Keep cleaned FASTQ files (uses much more disk space)", FALSE),
          checkboxInput("skip_space", "Skip the free-space check", FALSE),
          textInput("tool_dir", "Tool folder (blank = use PATH)", ""),
          uiOutput("tools_status"))))),
    tabPanel("3 · Run", br(),
      div(class = "card",
        uiOutput("run_gate"),
        actionButton("run", "Run pipeline", class = "btn-primary"),
        actionButton("cancel", "Cancel", class = "btn-danger"),
        actionButton("clear", "Delete previous outputs…"),
        br(), br(),
        uiOutput("progress"),
        p(class = "hint", "Long runs keep going if you close this page; reopen it to see progress."),
        tags$pre(class = "log", id = "logbox", textOutput("log")))),
    tabPanel("4 · Results", br(),
      div(class = "card",
        uiOutput("results_header"),
        fluidRow(column(5, selectInput("cmp", "Comparison", "All", width = "100%")),
                 column(5, br(), checkboxInput("sig_only", "Only significant genes (q < 0.05)", FALSE))),
        DTOutput("results"),
        downloadButton("dl", "Download table (CSV)")),
      div(class = "card", h4("Volcano plot"), plotOutput("volcano", height = 380)),
      div(class = "card", h4("Inputs for RTP-STAR"), uiOutput("rtp_inputs")))
  )),
  tabPanel("GENIST", genist_ui("genist")),
  tabPanel("RTP-STAR", rtpstar_ui("rtp"))),
  tabPanel("Drives", br(), fluidRow(column(8, drives_card, drives_help))),
  footer = div(class = "strym-footer",
    "StRym version 11. StRym was developed by Jordan Reed & Stephanie M. Stewart. ",
    "The platform is being maintained by the Sozzani lab at NC State University.")
)

server <- function(input, output, session) {
  roots <- NULL   # the modules still take this argument, but the Browse window (filebrowser.R) builds its own list of places
  genist_server("genist", roots, app_dir)
  tux_to_rtp <- reactiveVal(NULL)                                               # set by "Use them in the RTP-STAR tab" on the TUX results page
  rtpstar_server("rtp", roots, app_dir, from_tux = tux_to_rtp)
  design_server("design", roots, app_dir)
  observeEvent(input$go_design, updateNavbarPage(session, "module", "Experimental Design"))   # About page button
  # reopened page: open the module whose run is still going (TUX handles its own tab below)
  observeEvent(TRUE, {
    if (job_alive(job_load("analysis"), "run_analysis_cli.R")) updateNavbarPage(session, "module", "Experimental Design")
    else if (job_alive(job_load("rtpstar"), "run_rtpstar_cli.R")) updateNavbarPage(session, "module", "RTP-STAR")
    else if (job_alive(job_load("genist"), "run_genist_cli.R")) updateNavbarPage(session, "module", "GENIST")
  }, once = TRUE)
  br <- strym_browser_attach(input, output, session, key = "tuxbf")       # the Browse window (filebrowser.R), shared by the four pickers of the TUX tab

  # ---- drives ----
  drives_tick <- reactiveVal(0)
  drives_seen <- reactiveVal(0)
  drives <- reactive({                                                          # never waits for Windows: the table appears at once and fills in when the background answer arrives
    drives_tick(); invalidateLater(if (drives_pending()) 1500 else 15000)
    drives_for_table()
  })
  drive_status <- function(d) if (platform != "wsl") ifelse(d$state == "ready", "ready", "not responding: plug it in again, then press Refresh") else ifelse(d$state == "ready", paste("ready:", d$mount),
                              ifelse(d$state == "stale", "needs reconnecting: press Mount drive",
                              ifelse(d$state == "removed", "removed from Windows: press Release drive", "not mounted: press Mount drive")))
  output$drives <- renderTable({
    d <- drives()
    shiny::validate(shiny::need(nrow(d) > 0, if (platform == "wsl") "No Windows drives found. Press Refresh; if the drive is plugged in, see “If it does not work” below." else "No external drives found. Plug the drive in and press Refresh."))
    data.frame(Drive = if (platform == "wsl") paste0(d$letter, ":") else d$mount, Name = d$label,
               Size = ifelse(is.na(d$size), "", vapply(d$size, fmt_bytes, "")), Free = ifelse(is.na(d$free), "", vapply(d$free, fmt_bytes, "")),
               Status = drive_status(d))
  })
  observe({                                                                    # every drive except C: (the Windows system drive, always connected); only on WSL, where drives are mounted by hand
    req(platform == "wsl")
    d <- drives(); d <- d[d$letter != "C", , drop = FALSE]
    if (!nrow(d)) return(updateSelectInput(session, "mount_letter", choices = c("(no other drives found)" = "")))
    nm <- paste0(d$letter, ":", ifelse(nzchar(d$label), paste0(" ", d$label), ""), " (", c(ready = "ready", `not mounted` = "not mounted", stale = "needs reconnecting", removed = "removed")[d$state], ")")
    sel <- d$letter[d$state != "ready"][1]; if (is.na(sel)) sel <- d$letter[1]
    updateSelectInput(session, "mount_letter", choices = setNames(d$letter, nm), selected = sel)
  })
  mount_msg <- reactiveVal(NULL)
  observeEvent(input$mount, {
    if (!nzchar(input$mount_letter)) return(mount_msg(list(ok = FALSE, message = "There is no drive to mount. Plug the drive in, then press Refresh.")))
    nid <- showNotification("Connecting the drive…", duration = NULL, closeButton = FALSE, type = "message"); r <- mount_drive(input$mount_letter); removeNotification(nid); mount_msg(r); drives_tick(drives_tick() + 1)
  })
  observeEvent(input$unmount, {
    if (!nzchar(input$mount_letter)) return(mount_msg(list(ok = FALSE, message = "Pick the drive to release in the Drive box.")))
    busy <- job_alive(job_load("analysis"), "run_analysis_cli.R") || job_alive(job_load("rtpstar"), "run_rtpstar_cli.R") ||
            job_alive(job_load("genist"), "run_genist_cli.R") || (!is.null(load_state()$pid) && pid_alive(load_state()$pid))
    if (busy) return(mount_msg(list(ok = FALSE, message = "A StRym run is in progress, and it may be using this drive. Wait for it to finish (or cancel it) before releasing a drive.")))
    nid <- showNotification("Releasing the drive…", duration = NULL, closeButton = FALSE, type = "message"); r <- unmount_drive(input$mount_letter); removeNotification(nid); mount_msg(r); drives_tick(drives_tick() + 1)
  })
  observeEvent(input$drives_refresh, { refresh_drives(); drives_tick(drives_tick() + 1); mount_msg(NULL) })
  output$mount_msg <- renderUI({
    r <- mount_msg(); if (is.null(r)) return(NULL)
    div(class = if (r$ok) "ok" else "bad", r$message)
  })

  # ---- paths: Browse buttons fill the text boxes, which are the source of truth ----
  observeEvent(input$dir_in, br$open("dir", "Choose the folder with the FASTQ files", start = isolate(input$dir_txt), cb = function(p) updateTextInput(session, "dir_txt", value = p[1])))
  observeEvent(input$out_in, br$open("dir", "Choose the output folder", start = isolate(input$out_txt), cb = function(p) updateTextInput(session, "out_txt", value = p[1])))
  observeEvent(input$gtf, br$open("file", "Choose the annotation (GTF / GFF)", exts = c("gtf", "gff", "gff3"), start = isolate(input$gtf_txt), cb = function(p) updateTextInput(session, "gtf_txt", value = p[1])))
  observeEvent(input$fasta, br$open("file", "Choose the genome (FASTA)", exts = c("fa", "fasta", "fna"), start = isolate(input$fasta_txt), cb = function(p) updateTextInput(session, "fasta_txt", value = p[1])))

  inputs_folder <- reactive(clean_path(input$dir_txt))
  out_folder    <- reactive(clean_path(input$out_txt))
  gtf_path      <- reactive(clean_path(input$gtf_txt))
  fasta_path    <- reactive(clean_path(input$fasta_txt))
  paired <- reactive(input$read_type == "paired")
  tool_dir <- reactive(if (nzchar(trimws(input$tool_dir))) trimws(input$tool_dir) else NULL)

  # Suggest an output folder next to the data until the user types their own.
  auto_out <- reactiveVal("")
  observeEvent(list(inputs_folder(), drives_tick()), {
    f <- inputs_folder()
    if (dir.exists(f) && (!nzchar(input$out_txt) || identical(input$out_txt, auto_out()))) {
      auto_out(file.path(f, "tuxedo_results")); updateTextInput(session, "out_txt", value = auto_out())
    }
  })
  output$out_win <- renderText(if (nzchar(out_folder()) && strym_is_wsl()) paste("Windows path:", to_windows_path(out_folder())) else "")

  tools <- reactive(check_tools(tool_dir(), input$cleaner %||% "auto"))
  output$tools_status <- renderUI({
    t <- tools()
    if (all(t$found)) div(class = "ok", "✔ All 7 required tools found.")
    else div(class = "bad", "✘ Missing: ", paste(t$tool[!t$found], collapse = ", "))
  })

  # ---- samples: inferred from the folder, groups editable ----
  sheet <- reactiveVal(NULL)       # data.frame(sample, group, file, mate)
  rescan <- reactiveVal(0)
  observeEvent(input$rescan, rescan(rescan() + 1))
  observeEvent(list(inputs_folder(), paired(), rescan(), drives_tick()), {
    f <- inputs_folder()
    if (!nzchar(f) || !dir.exists(f)) return(sheet(NULL))
    s <- tryCatch(infer_samples(f, paired()), error = function(e) NULL)
    sheet(if (is.null(s) || !nrow(s)) NULL else
      data.frame(sample = s$sample, group = s$group, file = s$file, mate = s$mate, stringsAsFactors = FALSE))
  })
  output$samples <- renderDT({
    s <- sheet()
    shiny::validate(shiny::need(!is.null(s), "No FASTQ files found in this folder yet."))
    show <- data.frame(Sample = s$sample, Group = s$group,
                       Size = vapply(file.size(s$file), fmt_bytes, ""), stringsAsFactors = FALSE)
    datatable(show, rownames = FALSE, selection = "none", editable = list(target = "cell", disable = list(columns = c(0, 2))),
              options = list(pageLength = 25, dom = "tp"))
  })
  observeEvent(input$samples_cell_edit, {
    i <- input$samples_cell_edit; s <- sheet()
    if (!is.null(s) && i$col == 1) { s$group[i$row] <- as.character(i$value); sheet(s) }
  })
  eff_samples <- reactive({          # what will actually run: blanks dropped, groups normalized
    s <- sheet(); if (is.null(s)) return(NULL)
    normalize_samples(s[, c("group", "file", "mate")], paired())
  })

  # ---- validation ----
  rv <- reactiveValues(state = load_state(), done_tick = 0)
  validation <- reactive({
    rv$done_tick; drives_tick()
    tryCatch({
    errs <- character(); warns <- character(); info <- character()
    chk <- function(path, what, dir = FALSE) {
      if (!nzchar(path)) return(paste0("Choose ", what, "."))
      u <- unmounted_drive_problem(path); if (nzchar(u)) return(u)
      if (dir) { if (!dir.exists(path)) return(paste(what, "not found.")) }
      else if (!file.exists(path)) return(paste(what, "not found."))
      ""
    }
    add <- function(m) if (nzchar(m)) errs <<- c(errs, m)
    add(chk(inputs_folder(), "a folder with FASTQ files", TRUE))
    add(chk(gtf_path(), "an annotation (GTF/GFF) file"))
    add(chk(fasta_path(), "a genome FASTA file"))
    if (nzchar(fasta_path()) && grepl("\\.(gz|zip)$", fasta_path(), ignore.case = TRUE))
      errs <- c(errs, "The genome FASTA must be uncompressed (unzip it first).")
    if (!nzchar(out_folder())) errs <- c(errs, "Choose an output folder.")
    else {
      u <- unmounted_drive_problem(out_folder())
      if (nzchar(u)) errs <- c(errs, u)
      else { oc <- check_output_folder(out_folder()); errs <- c(errs, oc$errors); warns <- c(warns, oc$warnings) }
    }
    s <- eff_samples()
    in_ok <- nzchar(inputs_folder()) && dir.exists(inputs_folder()) && !nzchar(unmounted_drive_problem(inputs_folder()))
    if (!is.null(s) && in_ok) {
      v <- validate_samples(s, paired()); errs <- c(errs, v$errors); warns <- c(warns, v$warnings)
      est <- estimate_space(s, paired())
      if (nzchar(out_folder()) && !nzchar(unmounted_drive_problem(out_folder()))) {
        free <- free_space_bytes(out_folder())
        line <- sprintf("%d samples in %d groups · input ≈ %s · space needed ≈ %s · free on output drive: %s",
                        nrow(s), length(unique(s$group)), fmt_bytes(est$total), fmt_bytes(est$need), fmt_bytes(free))
        if (!is.na(free) && free < est$need && !isTRUE(input$skip_space))
          errs <- c(errs, paste("Not enough free space on the output drive.", line))
        else info <- c(info, line)
      }
    } else if (is.null(s) && in_ok) errs <- c(errs, paste0(
      "No FASTQ files found in that folder (files must end in .fastq/.fq, optionally .gz).",
      if (paired()) " Paired-end mates must be named A_x/B_x, x_R1/x_R2 or x_1/x_2." else ""))
    if (!all(tools()$found)) errs <- c(errs, "Required tools are missing (see Settings → Advanced).")
    list(errors = unique(errs), warnings = unique(warns), info = info)
    }, error = function(e)         # never leave the page hanging on an unexpected problem
      list(errors = paste("Unexpected problem while checking inputs:", conditionMessage(e)),
           warnings = character(), info = character()))
  })
  output$validation <- renderUI({
    v <- validation()
    if (!length(v$errors) && !length(v$warnings))
      return(tagList(div(class = "ok", "✔ Everything looks good."), lapply(v$info, function(i) div(class = "hint", i))))
    tagList(lapply(v$errors, function(e) div(class = "bad", "✘ ", e)),
            lapply(v$warnings, function(w) div(class = "warn", "⚠ ", w)),
            # notes (warnings) are advice, not blockers: the researcher can read them and run anyway. Errors cannot be bypassed.
            if (!length(v$errors) && length(v$warnings)) checkboxInput("ack_warn", "Run anyway: I have read the notes above", value = isTRUE(isolate(input$ack_warn))),
            lapply(v$info, function(i) div(class = "hint", i)))
  })

  # ---- run management (state lives in a file so a closed tab doesn't lose the run) ----
  status_tick <- reactiveTimer(1000)
  state_alive <- function(st) !is.null(st) && !is.null(st$pid) && pid_alive(st$pid)
  running <- reactive({ status_tick(); state_alive(rv$state) })
  status <- reactive({            # "idle" | "running" | "ok" | "failed" | "cancelled"
    status_tick()
    if (is.null(rv$state)) return("idle")
    if (state_alive(rv$state)) "running" else run_outcome(rv$state$log)
  })

  can_run <- reactive(length(validation()$errors) == 0 && (!length(validation()$warnings) || isTRUE(input$ack_warn)) && !running())
  output$run_gate <- renderUI({
    if (running()) div(class = "running", "Running…")
    else if (can_run()) div(class = "ok", "Ready to run.")
    else if (length(validation()$errors) == 0 && length(validation()$warnings)) div(class = "warn", "Read the notes on tab 1, then tick “Run anyway” if they are fine.")
    else div(class = "bad", "Not ready — see the checks on tab 1.")
  })
  observe({
    session$sendCustomMessage("tux-toggle", list(id = "run", enable = can_run()))
    session$sendCustomMessage("tux-toggle", list(id = "cancel", enable = running()))
  })

  observeEvent(input$run, {
    req(can_run())
    s <- eff_samples()
    cfg <- list(inputs_folder = inputs_folder(), output_folder = out_folder(),
                gtf_file = gtf_path(), fasta_file = fasta_path(), read_type = input$read_type,
                samples = s[, c("group", "file", "mate")],
                quality_thresh = input$q, min_seq_length = input$l, window_size = input$w,
                cleaner = input$cleaner %||% "auto", threads = whole_threads(input$threads), adapter_fasta = file.path(app_dir, "IlluminaAdaptorSeq.fasta"),
                tool_dir = if (is.null(tool_dir())) "" else tool_dir(),
                keep_cleaned = isTRUE(input$keep_cleaned), skip_space_check = isTRUE(input$skip_space))
    cfgfile <- file.path(STATE_DIR, "run_config.json")
    write_json(cfg, cfgfile, auto_unbox = TRUE, na = "null")
    logfile <- file.path(STATE_DIR, "last_run.log")
    extra <- if (!is.null(tool_dir())) paste(file.path(tool_dir(), c("ea-utils/clipper", ".", "hisat2", "samtools", "cufflinks")), collapse = ":")
    p <- process$new("Rscript", c(file.path(app_dir, "run_cli.R"), cfgfile),
                     stdout = logfile, stderr = "2>&1", cleanup = FALSE, cleanup_tree = FALSE,
                     env = c(Sys.getenv(), PATH = paste(c(Sys.getenv("PATH"), extra), collapse = ":")))
    st <- list(pid = p$get_pid(), log = logfile, out = out_folder(), started = format(Sys.time()))
    write_json(st, STATE_FILE, auto_unbox = TRUE); rv$state <- st
    updateTabsetPanel(session, "tabs", "3 · Run")
  })

  observeEvent(input$cancel, {
    if (state_alive(rv$state)) { kill_tree(rv$state$pid); showNotification("Run canceled.", type = "warning") }
  })

  observeEvent(input$clear, {
    req(dir.exists(out_folder()))
    showModal(modalDialog(title = "Delete previous outputs?",
      "This permanently deletes Cleaned, Aligned, Cufflinks, merged_asm and diff_out inside the OUTPUT folder. Your FASTQ files are never touched.",
      footer = tagList(modalButton("Keep"), actionButton("clear_ok", "Delete", class = "btn-danger"))))
  })
  observeEvent(input$clear_ok, {
    removeModal(); clear_outputs(out_folder()); rv$done_tick <- rv$done_tick + 1
    showNotification("Previous outputs deleted.")
  })

  log_text <- reactive({
    status_tick()
    if (is.null(rv$state) || !file.exists(rv$state$log)) return("")
    l <- readLines(rv$state$log, warn = FALSE)
    paste(tail(l, 400), collapse = "\n")
  })
  output$log <- renderText(log_text())

  output$progress <- renderUI({
    st <- status(); txt <- log_text()
    last <- function(re) { m <- regmatches(txt, gregexpr(re, txt))[[1]]; if (length(m)) tail(m, 1) else NA }
    stp <- last("\\[STEP [0-9]/4\\]"); smp <- last("\\[SAMPLE [0-9]+/[0-9]+\\]")
    n <- if (is.na(stp)) 0L else as.integer(sub("\\[STEP ([0-9])/4\\]", "\\1", stp))
    frac <- 0
    if (n %in% 1:2 && !is.na(smp)) {
      ij <- as.integer(strsplit(gsub("[^0-9/]", "", smp), "/")[[1]]); frac <- (ij[1] - 1) / ij[2]
    }
    pct <- switch(st, ok = 100, idle = 0, max(0, ((n - 1) + frac) / 4 * 100))
    cls <- switch(st, ok = "progress-bar-success", failed = "progress-bar-danger",
                  cancelled = "progress-bar-warning", "progress-bar-striped active")
    label <- switch(st, idle = "Not started", ok = "Finished — see the Results tab",
                    failed = "Failed — see the log below", cancelled = "Canceled",
                    if (n > 0) paste0("Step ", n, "/4: ", STEP_LABELS[n],
                                      if (n %in% 1:2 && !is.na(smp)) paste0(" · sample ", sub("\\[SAMPLE (.*)\\]", "\\1", smp)) else "")
                    else "Starting…")
    div(class = "progress", div(class = paste("progress-bar", cls), style = paste0("width:", max(pct, 4), "%"), label))
  })
  prev_status <- reactiveVal("idle")
  observeEvent(status(), {
    st <- status()
    if (prev_status() == "running" && st %in% c("ok", "failed", "cancelled")) {
      rv$done_tick <- rv$done_tick + 1
      if (!is.null(rv$state$out) && dir.exists(rv$state$out)) file.copy(rv$state$log, file.path(rv$state$out, "run.log"), overwrite = TRUE)
      if (st == "ok") updateTabsetPanel(session, "tabs", "4 · Results")
    }
    prev_status(st)
  })

  # ---- Results (read from the output folder, so past runs can be viewed too) ----
  diff_file <- reactive({
    rv$done_tick; drives_tick()
    f <- file.path(out_folder(), "diff_out", "gene_exp.diff")
    shiny::validate(shiny::need(file.exists(f), "No results in this output folder yet. Run the pipeline first."))
    f
  })
  diff_table <- reactive({
    d <- read.delim(diff_file(), check.names = FALSE, stringsAsFactors = FALSE)
    d$log2FC <- suppressWarnings(as.numeric(d[["log2(fold_change)"]]))
    d$comparison <- paste(d$sample_1, "vs", d$sample_2)
    d
  })
  observe({
    d <- tryCatch(diff_table(), error = function(e) NULL)
    updateSelectInput(session, "cmp", choices = c("All", if (!is.null(d)) unique(d$comparison)))
  })
  shown <- reactive({
    d <- diff_table()
    if (!is.null(input$cmp) && input$cmp != "All") d <- d[d$comparison == input$cmp, , drop = FALSE]
    if (isTRUE(input$sig_only)) d <- d[d$q_value < 0.05, , drop = FALSE]
    d[, c("gene_id", "locus", "sample_1", "sample_2", "value_1", "value_2", "log2FC", "p_value", "q_value", "significant")]
  })
  output$results_header <- renderUI({
    d <- diff_table()
    h4(sprintf("%d tests · %d significant (q < 0.05) · %s", nrow(d), sum(d$q_value < 0.05, na.rm = TRUE),
               to_windows_path(file.path(out_folder(), "diff_out"))))
  })
  output$results <- renderDT(datatable(shown(), rownames = FALSE, options = list(pageLength = 15)) |>
                               formatSignif(c("value_1", "value_2", "log2FC", "p_value", "q_value"), 3))
  # The genes file and the expression file RTP-STAR needs, written at the end of the Tux run (tux_write_rtpstar_inputs in tuxedo.R)
  rtp_dir <- reactive({ rv$done_tick; drives_tick(); file.path(out_folder(), "rtpstar_inputs") })
  rtp_files <- reactive({ d <- rtp_dir(); g <- file.path(d, "rtpstar_genes_from_tux.csv"); e <- file.path(d, "rtpstar_expression_from_tux.csv")
    if (file.exists(g) && file.exists(e)) list(genes = g, expr = e, readme = file.path(d, "README_rtpstar_inputs.txt")) })
  output$rtp_inputs <- renderUI({
    f <- rtp_files(); d <- rtp_dir(); nm <- file.path(d, "rtpstar_inputs_NOT_made.txt")
    if (!is.null(f)) {
      sm <- if (file.exists(f$readme)) { l <- readLines(f$readme, warn = FALSE); l[grepl("differentially expressed genes", l)][1] } else NA_character_
      return(tagList(
        p(class = "hint", "Made from these results for the RTP-STAR tab: a genes file (the differentially expressed genes) and a gene expression file (one column per sample)."),
        if (!is.na(sm)) div(class = "ok", style = "font-weight:400", sm),
        div(class = "path", to_windows_path(d)),
        downloadButton("dl_rtp_genes", "Download the genes file", class = "btn-default btn-sm"), downloadButton("dl_rtp_expr", "Download the expression file", class = "btn-default btn-sm"), " ",
        actionButton("rtp_use", "Use them in the RTP-STAR tab", class = "btn-primary btn-sm"),
        p(class = "hint", style = "margin-top:6px", "Other cutoffs? Open the RTP-STAR tab and build them again under “Start from your Tux results”.")))
    }
    if (file.exists(nm)) return(div(class = "warn", "⚠ ", paste(readLines(nm, warn = FALSE), collapse = " ")))
    div(class = "hint", "They are made at the end of a Tux run (genes file + gene expression file for RTP-STAR). None yet in this output folder.")
  })
  output$dl_rtp_genes <- downloadHandler("rtpstar_genes_from_tux.csv", function(file) file.copy(rtp_files()$genes, file, overwrite = TRUE))
  output$dl_rtp_expr <- downloadHandler("rtpstar_expression_from_tux.csv", function(file) file.copy(rtp_files()$expr, file, overwrite = TRUE))
  observeEvent(input$rtp_use, { f <- rtp_files(); req(f); tux_to_rtp(list(genes = f$genes, expr = f$expr)); updateNavbarPage(session, "module", "RTP-STAR") })
  output$dl <- downloadHandler("tuxedo_differential_expression.csv",
                               function(file) write.csv(shown(), file, row.names = FALSE))
  output$volcano <- renderPlot({
    d <- shown(); d <- d[is.finite(d$log2FC) & is.finite(d$p_value), ]
    sig <- d$q_value < 0.05
    plot(d$log2FC, -log10(pmax(d$p_value, 1e-300)), pch = 19, col = ifelse(sig, "#4156A1", "#B2B2B2"),
         xlab = "log2 fold change", ylab = "-log10(p)", main = "Differential expression")
    abline(v = c(-1, 1), lty = 3, col = "grey60")
    legend("topright", c("q < 0.05", "n.s."), pch = 19, col = c("#4156A1", "#B2B2B2"), bty = "n")
  })

  # Reopened page: remember the last run's output folder (so its results can be viewed), and
  # if that run is still going, jump to the Run tab.
  observeEvent(TRUE, {
    st <- rv$state
    if (!is.null(st$out) && !nzchar(isolate(input$out_txt))) {
      auto_out(st$out); updateTextInput(session, "out_txt", value = st$out)
    }
    if (state_alive(st)) updateTabsetPanel(session, "tabs", "3 · Run")
  }, once = TRUE)
}

shinyApp(ui, server)
