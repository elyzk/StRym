# mod_rtpstar.R -- the RTP-STAR tab of StRym (GENIE3-based regulatory networks with direction).
#
# rtpstar_ui(id)      four tabs: 1 Inputs (start from Tux results, required files, optional files, output), 2 Parameters, 3 Run, 4 Results.
# rtpstar_server(...) like the GENIST tab, plus "Start from your Tux results": it finds the genes file and the expression file a Tux run wrote (rtpstar_inputs/) or
#                     builds them with your cutoffs (rtpstar_inputs_from_tux() in rtpstar.R); the TUX results page can hand them over (from_tux).
# The algorithm is in rtpstar.R; this file is only the screen.

# Where the RTP-STAR results of inputs made from a Tux run go by default: the StRym project' 06_regulatory_networks/RTP-STAR when the files belong to a project,
# else a RTP-STAR_results folder beside the rtpstar_inputs folder (never inside it).
rtp_default_out <- function(genes_file) {
  d <- dirname(normalizePath(genes_file, winslash = "/", mustWork = FALSE)); up <- d
  for (i in 1:5) { up <- dirname(up); if (file.exists(file.path(up, "00_design", "design.json"))) return(file.path(up, "06_regulatory_networks", "RTP-STAR")) }
  file.path(dirname(d), "RTP-STAR_results")
}

rtpstar_ui <- function(id) {
  ns <- NS(id)
  tagList(cytoscape_note("RTP-STAR"), tabsetPanel(id = ns("tabs"),
    tabPanel("1 · Inputs", br(),
      fluidRow(
        column(6,
          div(class = "card", h4("Start from your Tux results"),
            pathUI(ns("tux"), "Tux results", type = "dir",
                   help = paste0("Your StRym project folder (the one that holds 01_quantification), or the Tux “diff_out” folder. StRym then fills in the genes file and the gene expression file below ",
                                 "from it (the differentially expressed genes, and the expression of every replicate). Nothing in that folder is changed.")),
            uiOutput(ns("tux_opts")), uiOutput(ns("tux_result"))),
          div(class = "card", h4("Required"),
            pathUI(ns("genes"), "Genes file",
                   help = "Excel/CSV file with the (differentially expressed) genes to use, in the first column (first row = header). Made for you when you start from Tux results above."),
            pathUI(ns("expr"), "Gene expression file",
                   help = "Expression used for the inference: gene names in column 1, then one column per sample (replicates in separate columns work best). May hold more genes than the genes file. A column named TF (1/0) marks the regulators; otherwise provide a TF file. Made for you when you start from Tux results above.")),
          div(class = "card", h4("Optional"),
            pathUI(ns("tf"), "Transcription factors file",
                   help = "Gene names in the first column; only these genes can regulate others. Optional: without it (and without a TF column) every gene is a potential regulator and the run takes longer."),
            pathUI(ns("sym"), "Gene name file", help = "Gene id in column 1, readable name in column 2."),
            pathUI(ns("clu"), "Clustering data file",
                   help = "Expression across spatial regions / conditions (gene names in column 1). Genes are clustered first and a network is inferred per cluster."),
            pathUI(ns("tc"), "Time course file",
                   help = "Expression over time (gene names in column 1, at least 3 time points). Used to label each edge as activating or inhibiting."),
            actionLink(ns("use_tf"), "Use the Arabidopsis TF list"), HTML("&nbsp;·&nbsp;"),
            actionLink(ns("no_tf"), "No TF list (optional; slower)"), HTML("&nbsp;·&nbsp;"),
            actionLink(ns("use_sym"), "Use the Arabidopsis gene names"))),
        column(6,
          div(class = "card", h4("Output"),
            pathUI(ns("out"), "Outputs folder", type = "dir"),
            textInput(ns("outfile"), "Results file name", "biograph.txt", width = "100%"),
            uiOutput(ns("overwrite_ui")),
            p(class = "hint", "Results: one line per edge (“A activates B”), plus the cluster files.")),
          div(class = "card", h4("Checks"), uiOutput(ns("validation"))),
          div(class = "card", h4("Try it"),
            actionButton(ns("example"), "Load the Arabidopsis example data"),
            p(class = "hint", "The example data shipped with TuxNet (226 genes; spatial clustering and a 5-point time course).")))))
    ,
    tabPanel("2 · Parameters", br(),
      fluidRow(
        column(6, div(class = "card", h4("Repeats"),
          numericInput(ns("iters"), "Number of iterations", 1, min = 1, step = 1),
          p(class = "hint", "How many times to run RTP-STAR and combine the results. Spatial clustering is random, so up to 100 iterations are recommended with it; with temporal clustering fewer than 10."),
          conditionalPanel("input.iters > 1", ns = ns,
            numericInput(ns("maxprop"), "Edge proportion to keep", 0.33, min = 0, max = 1, step = 0.01),
            p(class = "hint", "Edges that appear in at least this fraction of the iterations are kept (default 0.33)."))),
          div(class = "card", h4("Compute"),
            numericInput(ns("cores"), "Threads", max(1L, min(4L, parallel::detectCores() - 1L)),
                         min = 1, max = max(1L, parallel::detectCores()), step = 1),
            numericInput(ns("trees"), "Trees per random forest", 10000, min = 100, step = 1000),
            p(class = "hint", "TuxNet uses 10000. Fewer trees are faster but noisier."))),
        column(6,
          conditionalPanel("(input['clu-txt'] != '') || (input['preset-txt'] != '')", ns = ns,
            div(class = "card", h4("Clustering"),
              conditionalPanel("input['preset-txt'] == ''", ns = ns,
                selectInput(ns("ctype"), "Clustering type", c("Spatial (k-means)" = "S", "Temporal (hierarchical)" = "T"), "S"),
                numericInput(ns("seed"), "Clustering seed", 0, min = 0, step = 1),
                p(class = "hint", "0 = a new random start every run; any other number reproduces the same clusters.")),
              checkboxInput(ns("hubs"), "Connect the clusters through their hubs", TRUE),
              p(class = "hint", "The gene with the most outgoing edges in each cluster is a hub; hubs are linked to each other."))),
          conditionalPanel("(input['clu-txt'] == '') && (input['preset-txt'] == '')", ns = ns,
            div(class = "card", p(class = "hint", "Clustering options appear once a clustering data file is chosen on the Inputs tab."))),
          div(class = "card", h4("Advanced"),
            pathUI(ns("preset"), "Use my own clusters file instead of clustering",
                   help = "Gene in column 1, anything in column 2, cluster number in column 3 (as written by a previous run).")))),
      p(class = "hint", "Defaults are those of TuxNet's RTP-STAR tab.")),
    tabPanel("3 · Run", br(), job_panel_ui(ns, "Run RTP-STAR")),
    tabPanel("4 · Results", br(),
      div(class = "card", uiOutput(ns("res_header")), network_legend(), visNetwork::visNetworkOutput(ns("net"), height = "560px")),
      div(class = "card", h4("Edge table"), DTOutput(ns("edges")),
          downloadButton(ns("dl_edges"), "Download results file")),
      div(class = "card", h4("Clusters"), DTOutput(ns("clusters"))))))
}

# from_tux: optional reactiveVal set by the TUX tab ("Use them in the RTP-STAR tab") to list(genes = <genes file>, expr = <expression file>).
rtpstar_server <- function(id, roots, app_dir, from_tux = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    genes <- pathServer("genes", roots); expr <- pathServer("expr", roots); tf <- pathServer("tf", roots)
    sym <- pathServer("sym", roots); clu <- pathServer("clu", roots); tc <- pathServer("tc", roots)
    out <- pathServer("out", roots, "dir")
    preset <- pathServer("preset", roots)
    preset_path <- reactive(preset$path())

    ref <- function(f) file.path(app_dir, "reference", f)
    if (!is.null(from_tux)) observeEvent(from_tux(), { v <- from_tux(); req(is.list(v), file.exists(v$genes), file.exists(v$expr))
      if (!nzchar(out$path())) out$set(rtp_default_out(v$genes))
      genes$set(v$genes); expr$set(v$expr); showNotification("The genes file and the expression file made from the Tux results are filled in. Check the Parameters, then press Run.", duration = 8) }, ignoreNULL = TRUE)
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
      genes$set(ex("RTPSTAR_genesfile.xlsx")); expr$set(ex("RTPSTAR_expressionfile.xlsx"))
      clu$set(ex("RTPSTAR_clusteringfile.xlsx")); tc$set(ex("RTPSTAR_timecoursefile.xlsx"))
      tf$set(ref("Arabidopsis_TFs.xlsx")); sym$set(ref("GeneNames.xlsx"))
      out$set(file.path(app_dir, "example_output", "RTP-STAR"))
      ex_t0(Sys.time())                                                                # earlier example results in that folder may be replaced
      showNotification("Example data loaded. Check the Parameters, then press Run on the Run tab.", duration = 6)
    })
    observeEvent(genes$path(), {
      p <- genes$path()
      if (nzchar(p) && file.exists(p) && !nzchar(isolate(out$path()))) out$set(file.path(dirname(p), "RTP-STAR_results"))
    })
    last_auto <- reactiveVal("")                                                       # informative result name from the genes file (see mod_genist.R)
    observeEvent(genes$path(), {
      p <- genes$path(); if (!nzchar(p)) return()
      new <- strym_network_filename(tools::file_path_sans_ext(basename(p)), "RTPSTAR")
      cur <- isolate(input$outfile); o <- isolate(out$path())
      if (!is.null(cur) && cur %in% c("biograph.txt", isolate(last_auto())) && !(nzchar(o) && file.exists(file.path(o, cur))))
        updateTextInput(session, "outfile", value = new)
      last_auto(new)
    })

    # ---- start from Tux results: build the genes file and the expression file ----
    tux <- pathServer("tux", roots, "dir")
    tux_msg <- reactiveVal(NULL)
    tux_found <- reactive({
      p <- tux$path(); if (!nzchar(p) || !dir.exists(p)) return(NULL)
      tryCatch(tux_find_results(p), error = function(e) e)
    })
    prepared <- reactive({                                              # the files the Tux run itself wrote (rtpstar_inputs/), when they are there
      x <- tux_found(); if (is.null(x) || inherits(x, "error")) return(NULL)
      d <- file.path(dirname(x$diff_dir), "rtpstar_inputs"); g <- file.path(d, "rtpstar_genes_from_tux.csv"); e <- file.path(d, "rtpstar_expression_from_tux.csv")
      if (file.exists(g) && file.exists(e)) list(genes = g, expr = e, dir = d) else NULL
    })
    observeEvent(tux$path(), tux_msg(NULL), ignoreInit = TRUE)
    output$tux_opts <- renderUI({
      x <- tux_found(); if (is.null(x)) return(NULL)
      if (inherits(x, "error")) return(div(class = "bad", "✘ ", conditionMessage(x)))
      pr <- prepared()
      tagList(
        div(class = "hint", sprintf("Found the Tux results in %s.", to_windows_path(x$diff_dir))),
        if (!is.null(pr)) div(class = "ok", style = "font-weight:400", "✔ The Tux run already made these files for RTP-STAR (", span(class = "path", to_windows_path(pr$dir)), "). ",
                              actionButton(ns("tux_use"), "Use them", class = "btn-primary btn-sm")),
        div(class = "hint", if (!is.null(pr)) "Or build them again with other cutoffs:" else "Build the genes file and the expression file:"),
        fluidRow(
          column(6, numericInput(ns("tux_q"), "Significance cutoff (q-value)", 0.05, min = 0.0001, max = 1, step = 0.01)),
          column(6, numericInput(ns("tux_lfc"), "Smallest change (log2 fold change)", 1, min = 0, step = 0.25))),
        fluidRow(
          column(6, numericInput(ns("tux_fpkm"), "Smallest expression (FPKM) in some sample", 1, min = 0, step = 0.5)),
          column(6, numericInput(ns("tux_max"), "Most genes to keep (0 = no limit)", 500, min = 0, step = 100))),
        div(class = "hint", "A gene is chosen when Tux finds it significantly different between two of the groups by at least this much. If more genes pass than the limit, the ones that change most are kept."),
        actionButton(ns("tux_build"), "Build the genes file and the expression file", class = "btn-primary"))
    })
    output$tux_result <- renderUI({ m <- tux_msg(); if (!is.null(m) && !m$ok) div(class = "bad", style = "font-weight:400", "✘ ", m$text) })
    observeEvent(input$tux_use, { pr <- prepared(); req(pr); if (!nzchar(out$path())) out$set(rtp_default_out(pr$genes)); genes$set(pr$genes); expr$set(pr$expr); tux_msg(NULL)
      showNotification("Filled in from the Tux run. Check the Parameters, then press Run.", duration = 6) })
    observeEvent(input$tux_build, {
      x <- tux_found(); req(is.list(x), !inherits(x, "error"))
      num <- function(v, d) if (is.null(v) || is.na(v)) d else v
      odir <- out$path()
      if (!nzchar(odir)) odir <- if (!is.null(x$project_dir)) file.path(x$project_dir, "06_regulatory_networks", "RTP-STAR") else file.path(dirname(x$diff_dir), "RTP-STAR")
      r <- tryCatch(rtpstar_inputs_from_tux(tux$path(), file.path(odir, "inputs"), q_cutoff = num(input$tux_q, 0.05), min_log2fc = num(input$tux_lfc, 1),
                                            min_fpkm = num(input$tux_fpkm, 1), max_genes = num(input$tux_max, 500)), error = function(e) e)
      if (inherits(r, "error")) { tux_msg(list(ok = FALSE, text = conditionMessage(r))); return(invisible()) }
      if (!nzchar(out$path())) out$set(odir)
      genes$set(r$genes_file); expr$set(r$expression_file)
      tux_msg(list(ok = TRUE, text = r$summary))
      showNotification(paste0(r$summary, " Saved in ", to_windows_path(dirname(r$genes_file)), ". Check the Parameters, then press Run."), duration = 12)
    })

    results_tick <- reactiveVal(0)          # bumped whenever a run finishes (refreshes checks and results)
    result_path <- reactive(file.path(out$path(), input$outfile))
    exists_prev <- reactive({ results_tick(); nzchar(out$path()) && file.exists(result_path()) })
    output$overwrite_ui <- renderUI(if (exists_prev())
      checkboxInput(ns("overwrite"), "Replace the existing results in this folder", FALSE))

    check <- reactive({
      errs <- character(); warns <- character(); info <- character()
      p <- rtpstar_check_inputs(genes$path(), expr$path(), clu$path(), tc$path(), sym$path(), tf$path(), out$path(),
                                input$outfile, preset_path())
      errs <- c(errs, p)
      if (!length(p)) {
        g <- read_table_cached(genes$path()); e <- read_table_cached(expr$path())
        for (x in list(g, e)) if (inherits(x, "error")) errs <- c(errs, conditionMessage(x))
        if (!length(errs)) {
          gi <- unique(g[[1]][!is.na(g[[1]]) & nzchar(g[[1]])])
          hit <- gi[gi %in% e[[1]]]
          info <- c(info, sprintf("%d genes listed; %d found in the expression file (%d data columns).", length(gi), length(hit), ncol(e) - 1))
          if (length(hit) < 2) errs <- c(errs, "Fewer than 2 of the listed genes are in the expression file (do the gene names match?).")
          if ("TF" %in% names(e)) info <- c(info, "The expression file has a TF column: it defines the regulators.")
          else if (nzchar(tf$path())) { t <- read_table_cached(tf$path())
            if (!inherits(t, "error")) info <- c(info, sprintf("%d of the found genes are transcription factors.", sum(hit %in% t[[1]]))) }
          else warns <- c(warns, "No TF file or TF column (optional): every gene will be treated as a potential regulator, so the random forests are larger and the run takes longer. You can still run it.")
          if (nzchar(tc$path())) { t <- read_table_cached(tc$path())
            if (inherits(t, "error")) errs <- c(errs, conditionMessage(t))
            else { tm <- table_to_matrix(t)                                      # only numeric columns are time points
                   info <- c(info, sprintf("Time course: %d of the found genes, %d time points%s.", sum(hit %in% t[[1]]), ncol(tm),
                                           if (length(attr(tm, "dropped"))) sprintf(" (ignored as not data: %s)", paste(head(attr(tm, "dropped"), 5), collapse = ", ")) else ""))
                   if (ncol(tm) < 3) errs <- c(errs, "The time course needs at least 3 time points (numeric columns).") } }
          if (nzchar(clu$path())) { cl <- read_table_cached(clu$path())
            if (inherits(cl, "error")) errs <- c(errs, conditionMessage(cl))
            else { n <- sum(hit %in% cl[[1]]); info <- c(info, sprintf("Clustering: %d of the found genes are in the clustering file.", n))
                   if (n < 3) errs <- c(errs, "Fewer than 3 of the listed genes are in the clustering file.") } }
        }
      }
      for (pth in c(genes$path(), expr$path(), clu$path(), tc$path(), tf$path(), sym$path(), out$path(), preset_path()))
        if (nzchar(pth)) { u <- unmounted_drive_problem(pth); if (nzchar(u)) errs <- c(errs, u) }
      if (is.null(input$iters) || is.na(input$iters) || input$iters < 1 || input$iters != round(input$iters)) errs <- c(errs, "Number of iterations must be an integer of at least 1.")
      num_ok <- function(x, lo, hi = Inf) is.numeric(x) && length(x) == 1 && !is.na(x) && x >= lo && x <= hi                  # a blank or out-of-range box must be a message, never a crash in the run
      if (!num_ok(input$maxprop, 0, 1)) errs <- c(errs, "The edge proportion to keep must be a number between 0 and 1.")
      if (!num_ok(input$trees, 100)) errs <- c(errs, "Trees per random forest must be a number of at least 100.")
      if (!num_ok(input$cores, 1)) errs <- c(errs, "Threads must be a whole number of at least 1.")
      if (!num_ok(input$seed, 0)) errs <- c(errs, "The clustering seed must be a number of at least 0.")
      if (exists_prev() && !isTRUE(input$overwrite)) warns <- c(warns, "Results already exist in this folder: tick “Replace the existing results”, or choose another folder.")
      list(errors = unique(errs), warnings = unique(warns), info = info)
    })
    output$validation <- renderUI({
      v <- check()
      if (!length(v$errors) && !length(v$warnings))
        return(tagList(div(class = "ok", "✔ Everything looks good."), lapply(v$info, function(i) div(class = "hint", i))))
      tagList(lapply(v$errors, function(e) div(class = "bad", "✘ ", e)), lapply(v$warnings, function(w) div(class = "warn", "⚠ ", w)),
              lapply(v$info, function(i) div(class = "hint", i)))
    })
    # the "no TF" note is only advice; results-exist is a blocker unless confirmed
    can_run <- reactive({
      v <- check(); blockers <- v$warnings[grepl("already exist", v$warnings)]
      !length(v$errors) && !length(blockers)
    })

    job <- job_controller("rtpstar", "run_rtpstar_cli.R", input, output, session, can_run, app_dir,
                          on_finished = function(st) { results_tick(results_tick() + 1)
                            if (st == "ok") updateTabsetPanel(session, "tabs", "4 · Results") })
    observeEvent(input$run, {
      req(can_run())
      use_clu <- nzchar(clu$path()) || nzchar(preset_path())
      cfg <- list(genes_file = genes$path(), expression_file = expr$path(), output_folder = out$path(),
                  output_file = input$outfile, clustering_file = clu$path(), time_course_file = tc$path(),
                  symbol_file = sym$path(), tf_file = tf$path(), preset_clusters_file = preset_path(),
                  numiters = input$iters, maxprop = if (is.na(input$maxprop)) 0.33 else input$maxprop,
                  connect_hubs = if (use_clu) isTRUE(input$hubs) else FALSE,
                  clustering_seed = if (use_clu && !is.na(input$seed) && input$seed > 0) input$seed else NULL,
                  clustering_type = if (use_clu) input$ctype else "S",
                  ntrees = input$trees, cores = whole_threads(input$cores))
      job$start(cfg, out$path(), result_path())
      updateTabsetPanel(session, "tabs", "3 · Run")
    })
    observeEvent(TRUE, {
      st <- job$state()
      cf <- file.path(JOB_DIR, "rtpstar_config.json")
      cfg <- if (!is.null(st) && file.exists(cf)) tryCatch(jsonlite::fromJSON(cf), error = function(e) NULL)
      if (!is.null(cfg)) {                              # restore the last run's inputs and settings
        s <- function(x) if (is.null(x) || is.na(x)) "" else as.character(x)
        genes$set(s(cfg$genes_file)); expr$set(s(cfg$expression_file)); tf$set(s(cfg$tf_file)); sym$set(s(cfg$symbol_file))
        clu$set(s(cfg$clustering_file)); tc$set(s(cfg$time_course_file)); preset$set(s(cfg$preset_clusters_file))
        out$set(s(cfg$output_folder)); updateTextInput(session, "outfile", value = s(cfg$output_file))
        updateNumericInput(session, "iters", value = cfg$numiters); updateNumericInput(session, "maxprop", value = cfg$maxprop)
        updateNumericInput(session, "trees", value = cfg$ntrees); updateNumericInput(session, "cores", value = cfg$cores)
        if (!is.null(cfg$clustering_type)) updateSelectInput(session, "ctype", selected = cfg$clustering_type)
        if (!is.null(cfg$connect_hubs)) updateCheckboxInput(session, "hubs", value = isTRUE(cfg$connect_hubs))
        updateNumericInput(session, "seed", value = if (is.null(cfg$clustering_seed)) 0 else cfg$clustering_seed)
      }
      if (job_alive(st, "run_rtpstar_cli.R")) updateTabsetPanel(session, "tabs", "3 · Run")
    }, once = TRUE)

    # ---- results: lines "A activates B" or "A activates B 0.67" ----
    net <- reactive({
      results_tick(); req(nzchar(out$path()))
      f <- result_path()
      shiny::validate(shiny::need(file.exists(f), "No results in this folder yet. Run RTP-STAR first."))
      l <- trimws(readLines(f, warn = FALSE)); l <- l[nzchar(l)]
      shiny::validate(shiny::need(length(l) > 0, "The run finished but found no edges."))
      p <- strsplit(l, " ", fixed = TRUE)
      data.frame(from = vapply(p, `[`, "", 1), type = vapply(p, `[`, "", 2), to = vapply(p, `[`, "", 3),
                 weight = vapply(p, function(x) if (length(x) >= 4) suppressWarnings(as.numeric(x[4])) else NA_real_, numeric(1)),
                 stringsAsFactors = FALSE)
    })
    output$res_header <- renderUI({
      d <- net()
      h4(sprintf("%d genes · %d edges (%d activating, %d inhibiting, %d undetermined)", length(unique(c(d$from, d$to))), nrow(d),
                 sum(d$type == "activates"), sum(d$type == "inhibits"), sum(d$type == "regulates")))
    })
    output$net <- visNetwork::renderVisNetwork({ d <- net(); network_widget(d) })
    output$edges <- renderDT({
      d <- net(); names(d) <- c("Regulator", "Effect", "Target", "Proportion of runs")
      if (all(is.na(d[[4]]))) d <- d[-4]
      datatable(d, rownames = FALSE, options = list(pageLength = 10))
    })
    output$dl_edges <- downloadHandler(function() input$outfile, function(file) file.copy(result_path(), file))
    output$clusters <- renderDT({
      results_tick(); f <- file.path(out$path(), "clusters.csv")
      shiny::validate(shiny::need(file.exists(f), "No cluster file (clustering was not used, or several iterations were run: see clusters1.csv, …)."))
      datatable(read.csv(f, check.names = FALSE), rownames = FALSE, options = list(pageLength = 10))
    })
  })
}
