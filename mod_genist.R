# mod_genist.R -- the GENIST tab of StRym (gene regulatory networks from time-course data).
#
# genist_ui(id)      four tabs: 1 Inputs (start from Tux results, required files, optional files, output), 2 Parameters, 3 Run, 4 Results.
# genist_server(...) fills the boxes (example data, Tux results, the Arabidopsis TF list and gene names), checks the files while you type (check(): errors block the
#                    Run button, "notes" such as a missing TF list are only advice), starts run_genist_cli.R in the background, and shows the network (visNetwork),
#                    the edge table and the clusters. The network table is named after the genes file unless you type a name (strym_network_filename()).
# The algorithm is in genist.R; this file is only the screen.

cytoscape_note <- function(what) div(class = "card cyto-note", h4(paste("Gene Regulatory Networks with", what)),
  p(style = "margin:0", HTML(paste0("The outputs of ", what, " are designed to be explored in <strong>Cytoscape</strong>. Download the network table and open it in Cytoscape ",
    "(File → Import → Network from File) to lay out, style and analyze the network."))))

genist_ui <- function(id) {
  ns <- NS(id)
  tagList(cytoscape_note("GENIST"), tabsetPanel(id = ns("tabs"),
    tabPanel("1 · Inputs", br(),
      fluidRow(
        column(6,
          div(class = "card", h4("Start from your Tux results"),
            pathUI(ns("tux"), "Tux results", type = "dir",
                   help = paste0("Your StRym project folder (the one that holds 01_quantification), or the Tux “diff_out” folder. ",
                                 "StRym then fills in the genes file and the time course file below from it. Nothing in that folder is changed.")),
            uiOutput(ns("tux_opts")), uiOutput(ns("tux_result"))),
          div(class = "card", h4("Required"),
            pathUI(ns("genes"), "Genes file",
                   help = "Excel/CSV file with the genes to include in the network in the first column (first row = header). Made for you when you start from Tux results above."),
            pathUI(ns("tc"), "Time course file",
                   help = "Expression of each gene across time or developmental points: gene names in column 1, one column per time point, earliest first (at least 3). It may contain many more genes than the genes file. Made for you when you start from Tux results above.")),
          div(class = "card", h4("Optional"),
            pathUI(ns("clu"), "Clustering data file",
                   help = "Expression across spatial regions / conditions / treatments (gene names in column 1). When given, genes are first grouped into clusters and a network is built per cluster, then joined through cluster hubs."),
            pathUI(ns("tf"), "Transcription factors file",
                   help = "Gene names in the first column. When given, only these genes can act as regulators. Optional: without it every gene may regulate every other gene, and the run takes longer."),
            pathUI(ns("sym"), "Gene name file",
                   help = "Gene id in column 1, a readable name (symbol) in column 2; used in the result tables and plot."),
            actionLink(ns("use_tf"), "Use the Arabidopsis TF list"), HTML("&nbsp;·&nbsp;"),
            actionLink(ns("no_tf"), "No TF list (optional; slower)"), HTML("&nbsp;·&nbsp;"),
            actionLink(ns("use_sym"), "Use the Arabidopsis gene names"))),
        column(6,
          div(class = "card", h4("Output"),
            pathUI(ns("out"), "Outputs folder", type = "dir"),
            textInput(ns("outfile"), "Output file name", "cytoscape_table.txt", width = "100%"),
            uiOutput(ns("overwrite_ui")),
            p(class = "hint", "Results: the network table (Cytoscape format: source, target, width, sign) and clusters.xlsx.")),
          div(class = "card", h4("Checks"), uiOutput(ns("validation"))),
          div(class = "card", h4("Try it"),
            actionButton(ns("example"), "Load the Arabidopsis example data"),
            p(class = "hint", "The example data shipped with TuxNet (genes, time course, clustering).")))))
    ,
    tabPanel("2 · Parameters", br(),
      fluidRow(
        column(6, div(class = "card", h4("Network inference"),
          selectInput(ns("lapse"), "Time lapse",
                      c("Next time point (1)" = "1", "Same time point (0)" = "0", "Either (0 and 1)" = "both"), "1"),
          p(class = "hint", "Which regulator activity can explain a gene's change: 1 = a regulator that changed in the PREVIOUS time point; 0 = in the SAME time point."),
          numericInput(ns("fc"), "Regulation fold-change threshold", 1.1, min = 1, step = 0.05),
          p(class = "hint", "Minimum fold change between consecutive time points for a gene to count as changing (default 1.1)."),
          numericInput(ns("pct"), "Regulation time percent", 0.1, min = 0, max = 1, step = 0.05),
          p(class = "hint", "Minimum fraction of the time points in which a regulator and a target must show matching changes (default 0.1)."))),
        column(6, div(class = "card", h4("Discretization and filtering"),
          numericInput(ns("levels"), "Discretization levels", 3, min = 2, step = 1),
          p(class = "hint", "Number of expression levels used to compute the Bayesian scores; integer, recommended 2 or 3 (default 3)."),
          numericInput(ns("low"), "Bottom fraction of regulations to drop", 0.05, min = 0, max = 1, step = 0.05),
          p(class = "hint", "Regulations weaker than this fraction of the average weight are discarded (default 0.05)."),
          h4("Compute"),
          numericInput(ns("cores"), "Threads", max(1L, min(4L, parallel::detectCores() - 1L)),
                       min = 1, max = max(1L, parallel::detectCores()), step = 1)))),
      p(class = "hint", "Defaults are those shown on the original TuxNet GENIST tab.")),
    tabPanel("3 · Run", br(), job_panel_ui(ns, "Run GENIST")),
    tabPanel("4 · Results", br(),
      div(class = "card", uiOutput(ns("res_header")), network_legend(), visNetwork::visNetworkOutput(ns("net"), height = "560px")),
      div(class = "card", h4("Network table"), DTOutput(ns("edges")),
          downloadButton(ns("dl_edges"), "Download network table (for Cytoscape)")),
      div(class = "card", h4("Clusters"), DTOutput(ns("clusters"))))))
}

genist_server <- function(id, roots, app_dir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    genes <- pathServer("genes", roots); tc <- pathServer("tc", roots); clu <- pathServer("clu", roots)
    tf <- pathServer("tf", roots); sym <- pathServer("sym", roots); out <- pathServer("out", roots, "dir")

    ref <- function(f) file.path(app_dir, "reference", f)
    observeEvent(TRUE, {                              # like the original tab, the Arabidopsis TF list and gene names are filled in at the start
      if (file.exists(ref("Arabidopsis_TFs.xlsx"))) tf$set(ref("Arabidopsis_TFs.xlsx"))
      if (file.exists(ref("GeneNames.xlsx"))) sym$set(ref("GeneNames.xlsx"))
    }, once = TRUE)
    observeEvent(input$use_tf,  tf$set(ref("Arabidopsis_TFs.xlsx")))
    observeEvent(input$no_tf,   tf$set(""))
    observeEvent(input$use_sym, sym$set(ref("GeneNames.xlsx")))
    ex_t0 <- reactiveVal(NULL)                                                   # after the example is loaded, tick "Replace the existing results" once that box has been redrawn
    observe({ t0 <- ex_t0(); req(!is.null(t0)); el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (el > 10) { ex_t0(NULL); return() }
      invalidateLater(400)
      if (el >= 1.5 && !is.null(input$overwrite)) { updateCheckboxInput(session, "overwrite", value = TRUE); ex_t0(NULL) } })
    observeEvent(input$example, {
      ex <- function(f) file.path(app_dir, "example_data", f)
      genes$set(ex("GENIST_genesfile.xlsx")); tc$set(ex("GENIST_timecoursefile.xlsx"))
      clu$set(ex("GENIST_clusteringfile.xlsx")); tf$set(ref("Arabidopsis_TFs.xlsx")); sym$set(ref("GeneNames.xlsx"))
      out$set(file.path(app_dir, "example_output", "GENIST"))
      ex_t0(Sys.time())                                                                # earlier example results in that folder may be replaced
      showNotification("Example data loaded. Check the Parameters, then press Run on the Run tab.", duration = 6)
    })
    # suggest an output folder next to the genes file until the user sets one
    observeEvent(genes$path(), {
      p <- genes$path()
      if (nzchar(p) && file.exists(p) && !nzchar(isolate(out$path()))) out$set(file.path(dirname(p), "GENIST_results"))
    })
    # an informative network table name from the genes file ("<genes file>_GENIST_network_for_Cytoscape.txt") instead of cytoscape_table.txt,
    # unless the researcher typed their own name or results already exist under the current name (their Results tab keeps working)
    last_auto <- reactiveVal("")
    observeEvent(genes$path(), {
      p <- genes$path(); if (!nzchar(p)) return()
      new <- strym_network_filename(tools::file_path_sans_ext(basename(p)), "GENIST")
      cur <- isolate(input$outfile); o <- isolate(out$path())
      if (!is.null(cur) && cur %in% c("cytoscape_table.txt", isolate(last_auto())) && !(nzchar(o) && file.exists(file.path(o, cur))))
        updateTextInput(session, "outfile", value = new)
      last_auto(new)
    })

    # ---- start from Tux results: build the genes file and the time course file ----
    tux <- pathServer("tux", roots, "dir")
    tux_msg <- reactiveVal(NULL)
    tux_found <- reactive({
      p <- tux$path(); if (!nzchar(p) || !dir.exists(p)) return(NULL)
      tryCatch({ f <- tux_find_results(p); list(found = f, groups = tux_groups(f)) }, error = function(e) e)
    })
    observeEvent(tux$path(), tux_msg(NULL), ignoreInit = TRUE)
    output$tux_opts <- renderUI({
      x <- tux_found(); if (is.null(x)) return(NULL)
      if (inherits(x, "error")) return(div(class = "bad", "✘ ", conditionMessage(x)))
      g <- x$groups; s <- g$default
      tagList(
        div(class = "hint", sprintf("Found the Tux results in %s.", to_windows_path(x$found$diff_dir))),
        if (g$guessed) div(class = "warn", "⚠ No StRym design record was found next to these results, so the time order below is Tux's own order. ",
                           "Check it and edit it if the time points are in another order."),
        if (length(g$series) > 1) selectInput(ns("tux_series"), "Condition to build the time course from", names(g$series), s, width = "100%"),
        textInput(ns("tux_order"), "Time points in order (earliest first)", paste(names(g$series[[s]]), collapse = ", "), width = "100%"),
        div(class = "hint", "Separate the names with commas. At least 3 time points are needed."),
        fluidRow(
          column(6, numericInput(ns("tux_q"), "Significance cutoff (q-value)", 0.05, min = 0.0001, max = 1, step = 0.01)),
          column(6, numericInput(ns("tux_lfc"), "Smallest change (log2 fold change)", 1, min = 0, step = 0.25))),
        fluidRow(
          column(6, numericInput(ns("tux_fpkm"), "Smallest expression (FPKM) at some time point", 1, min = 0, step = 0.5)),
          column(6, numericInput(ns("tux_max"), "Most genes to keep (0 = no limit)", 500, min = 0, step = 100))),
        div(class = "hint", "A gene is chosen when Tux finds it significantly different between two of the time points by at least this much. ",
            "If more genes pass than the limit, the ones that change most are kept. We recommend limiting the number of genes to 500."),
        actionButton(ns("tux_build"), "Build the genes file and the time course file", class = "btn-primary"))
    })
    observeEvent(input$tux_series, {
      x <- tux_found()
      if (is.list(x) && !inherits(x, "error") && input$tux_series %in% names(x$groups$series))
        updateTextInput(session, "tux_order", value = paste(names(x$groups$series[[input$tux_series]]), collapse = ", "))
    }, ignoreInit = TRUE)
    output$tux_result <- renderUI({ m <- tux_msg(); if (!is.null(m) && !m$ok) div(class = "bad", style = "font-weight:400", "✘ ", m$text) })   # only problems are shown here: the boxes below fill in on success
    observeEvent(input$tux_build, {
      x <- tux_found(); req(is.list(x), !inherits(x, "error"))
      num <- function(v, d) if (is.null(v) || is.na(v)) d else v
      odir <- out$path()
      if (!nzchar(odir)) {                                                    # default: next to the other network results of the project
        odir <- if (!is.null(x$found$project_dir)) file.path(x$found$project_dir, "06_regulatory_networks", "GENIST") else file.path(dirname(x$found$diff_dir), "GENIST")
      }
      r <- tryCatch(genist_inputs_from_tux(tux$path(), file.path(odir, "inputs"), order = strsplit(if (is.null(input$tux_order)) "" else input$tux_order, ",")[[1]],
                                           series = input$tux_series, q_cutoff = num(input$tux_q, 0.05), min_log2fc = num(input$tux_lfc, 1),
                                           min_fpkm = num(input$tux_fpkm, 1), max_genes = num(input$tux_max, 500)),
                    error = function(e) e)
      if (inherits(r, "error")) { tux_msg(list(ok = FALSE, text = conditionMessage(r))); return(invisible()) }
      if (!nzchar(out$path())) out$set(odir)
      genes$set(r$genes_file); tc$set(r$time_course_file)
      if (!nzchar(tf$path()) && file.exists(ref("Arabidopsis_TFs.xlsx"))) tf$set(ref("Arabidopsis_TFs.xlsx"))
      if (!nzchar(sym$path()) && file.exists(ref("GeneNames.xlsx"))) sym$set(ref("GeneNames.xlsx"))
      tux_msg(list(ok = TRUE, text = paste0(r$summary, " Saved in ", to_windows_path(dirname(r$genes_file)),
                                           ". The genes file and time course file below are filled in; check the Parameters, then press Run.")))
    })

    results_tick <- reactiveVal(0)          # bumped whenever a run finishes (refreshes checks and results)
    result_path <- reactive(file.path(out$path(), input$outfile))
    exists_prev <- reactive({ results_tick(); nzchar(out$path()) && file.exists(result_path()) })
    output$overwrite_ui <- renderUI(if (exists_prev())
      checkboxInput(ns("overwrite"), "Replace the existing results in this folder", FALSE))

    # ---- checks, with a short preview of what was found ----
    check <- reactive({
      errs <- character(); warns <- character(); info <- character()
      p <- genist_check_inputs(genes$path(), tc$path(), clu$path(), tf$path(), sym$path(), out$path(), input$outfile)
      errs <- c(errs, p)
      if (!length(p)) {
        g <- read_table_cached(genes$path()); t <- read_table_cached(tc$path())
        for (x in list(g, t)) if (inherits(x, "error")) errs <- c(errs, conditionMessage(x))
        if (!length(errs)) {
          gi <- unique(g[[1]][!is.na(g[[1]]) & nzchar(g[[1]])]); ti <- t[[1]]
          tm <- table_to_matrix(t)                                                    # only the numeric columns are time points (a symbol column is not one)
          in_t <- sum(gi %in% ti); info <- c(info, sprintf("%d genes listed; %d found in the time course (%d time points).", length(gi), in_t, ncol(tm)))
          if (length(attr(tm, "dropped"))) info <- c(info, sprintf("Ignored as not data: %s.", paste(head(attr(tm, "dropped"), 5), collapse = ", ")))
          if (ncol(tm) < 3) errs <- c(errs, "The time course needs at least 3 time points (numeric columns).")
          if (in_t < 2) errs <- c(errs, "Fewer than 2 of the listed genes are in the time course file (do the gene names match?).")
          if (nzchar(clu$path())) {
            cl <- read_table_cached(clu$path())
            if (inherits(cl, "error")) errs <- c(errs, conditionMessage(cl))
            else { info <- c(info, sprintf("%d of the listed genes are in the clustering file (%d conditions).", sum(gi %in% cl[[1]]), ncol(cl) - 1))
                   if (sum(gi %in% cl[[1]]) < 2) errs <- c(errs, "Fewer than 2 of the listed genes are in the clustering file.") }
          }
          if (nzchar(tf$path())) { tfs <- read_table_cached(tf$path())
            if (!inherits(tfs, "error")) info <- c(info, sprintf("%d of the listed genes are transcription factors.", sum(gi %in% tfs[[1]]))) }
        }
      }
      for (pth in c(genes$path(), tc$path(), clu$path(), tf$path(), sym$path(), out$path()))
        if (nzchar(pth)) { u <- unmounted_drive_problem(pth); if (nzchar(u)) errs <- c(errs, u) }
      num_ok <- function(x, lo, hi = Inf) is.numeric(x) && length(x) == 1 && !is.na(x) && x >= lo && x <= hi                  # a blank or out-of-range box must be a message, never a crash in the run
      if (!num_ok(input$fc, 1)) errs <- c(errs, "Fold-change threshold must be a number of at least 1.")
      if (!num_ok(input$pct, 0, 1)) errs <- c(errs, "Regulation time percent must be a number between 0 and 1.")
      if (!num_ok(input$low, 0, 1)) errs <- c(errs, "The fraction of weak regulations to drop must be a number between 0 and 1.")
      if (!num_ok(input$cores, 1)) errs <- c(errs, "Cores must be a whole number of at least 1.")
      if (is.na(input$levels) || input$levels < 2 || input$levels != round(input$levels)) errs <- c(errs, "Discretization levels must be an integer of at least 2.")
      if (exists_prev() && !isTRUE(input$overwrite)) warns <- c(warns, "Results already exist in this folder: tick “Replace the existing results”, or choose another folder.")
      notes <- if (!nzchar(tf$path())) "No transcription factor list: every gene can act as a regulator, so the network takes longer to infer (much longer for more than about 100 genes). You can still run it." else character()   # advice only, never blocks
      list(errors = unique(errs), warnings = warns, info = info, notes = notes)
    })
    output$validation <- renderUI({
      v <- check(); nts <- lapply(v$notes, function(w) div(class = "warn", "⚠ ", w))
      if (!length(v$errors) && !length(v$warnings))
        return(tagList(div(class = "ok", "✔ Everything looks good."), nts, lapply(v$info, function(i) div(class = "hint", i))))
      tagList(lapply(v$errors, function(e) div(class = "bad", "✘ ", e)), lapply(v$warnings, function(w) div(class = "warn", "⚠ ", w)), nts,
              lapply(v$info, function(i) div(class = "hint", i)))
    })
    can_run <- reactive({ v <- check(); !length(v$errors) && !length(v$warnings) })

    # ---- run ----
    job <- job_controller("genist", "run_genist_cli.R", input, output, session, can_run, app_dir,
                          on_finished = function(st) { results_tick(results_tick() + 1)
                            if (st == "ok") updateTabsetPanel(session, "tabs", "4 · Results") })
    observeEvent(input$run, {
      req(can_run())
      cfg <- list(genes_file = genes$path(), time_course_file = tc$path(), clustering_file = clu$path(),
                  tf_file = tf$path(), symbol_file = sym$path(), output_folder = out$path(), output_file = input$outfile,
                  time_lapse = if (input$lapse == "both") "both" else as.numeric(input$lapse),
                  fc_thresh = input$fc, time_percent = input$pct, n_levels = input$levels, low_conn = input$low,
                  cores = whole_threads(input$cores))
      job$start(cfg, out$path(), result_path())
      updateTabsetPanel(session, "tabs", "3 · Run")
    })
    # reopened page: remember the last run's output location; jump to the Run tab if it is still going
    observeEvent(TRUE, {
      st <- job$state()
      cf <- file.path(JOB_DIR, "genist_config.json")
      cfg <- if (!is.null(st) && file.exists(cf)) tryCatch(jsonlite::fromJSON(cf), error = function(e) NULL)
      if (!is.null(cfg)) {                              # restore the last run's inputs and settings
        s <- function(x) if (is.null(x) || is.na(x)) "" else as.character(x)
        genes$set(s(cfg$genes_file)); tc$set(s(cfg$time_course_file)); clu$set(s(cfg$clustering_file))
        tf$set(s(cfg$tf_file)); sym$set(s(cfg$symbol_file)); out$set(s(cfg$output_folder))
        updateTextInput(session, "outfile", value = s(cfg$output_file))
        updateSelectInput(session, "lapse", selected = s(cfg$time_lapse))
        updateNumericInput(session, "fc", value = cfg$fc_thresh); updateNumericInput(session, "pct", value = cfg$time_percent)
        updateNumericInput(session, "levels", value = cfg$n_levels); updateNumericInput(session, "low", value = cfg$low_conn)
        updateNumericInput(session, "cores", value = cfg$cores)
      }
      if (job_alive(st, "run_genist_cli.R")) updateTabsetPanel(session, "tabs", "3 · Run")
    }, once = TRUE)

    # ---- results ----
    net <- reactive({
      results_tick(); req(nzchar(out$path()))
      f <- result_path()
      shiny::validate(shiny::need(file.exists(f), "No results in this folder yet. Run GENIST first."))
      d <- read.delim(f, stringsAsFactors = FALSE)
      d$type <- ifelse(d$sign > 0, "activates", ifelse(d$sign < 0, "inhibits", "regulates"))
      d
    })
    output$res_header <- renderUI({
      d <- net(); e <- d[d$width > 0, ]
      h4(sprintf("%d genes · %d regulations (%d activating, %d inhibiting, %d undetermined)",
                 length(unique(c(d$source, d$target))), nrow(e), sum(e$sign > 0), sum(e$sign < 0), sum(e$sign == 0)))
    })
    output$net <- visNetwork::renderVisNetwork({
      d <- net(); e <- d[d$width > 0, ]
      iso <- setdiff(unique(c(d$source, d$target)), c(e$source, e$target))
      network_widget(data.frame(from = e$source, to = e$target, type = e$type, weight = e$width, stringsAsFactors = FALSE), iso)
    })
    output$edges <- renderDT({
      d <- net(); d <- d[d$width > 0, c("source", "target", "width", "sign", "type")]
      datatable(d, rownames = FALSE, options = list(pageLength = 10)) |> formatSignif("width", 3)
    })
    output$dl_edges <- downloadHandler(function() input$outfile, function(file) file.copy(result_path(), file))
    output$clusters <- renderDT({
      results_tick(); f <- file.path(out$path(), "clusters.xlsx")
      shiny::validate(shiny::need(file.exists(f), "No cluster file yet."))
      datatable(as.data.frame(readxl::read_excel(f, sheet = 2, .name_repair = "minimal")), rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE))
    })
  })
}
