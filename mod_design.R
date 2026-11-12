# mod_design.R -- the "Experimental Design" page: describe the experiment, review the plan, run, browse results.

design_ui <- function(id) {
  ns <- NS(id)
  steps <- strym_steps(); groups <- unique(vapply(steps, `[[`, "", "group"))
  plan_rows <- lapply(groups, function(g) tagList(
    div(class = "plan-group", g),
    lapply(Filter(function(s) s$group == g, steps), function(s) div(class = "plan-row",
      div(class = "plan-main",
          checkboxInput(ns(paste0("sel_", s$id)), label = span(tags$b(s$label), uiOutput(ns(paste0("badge_", s$id)), inline = TRUE)), value = FALSE, width = "100%"),
          div(class = "hint", s$desc), uiOutput(ns(paste0("note_", s$id)))),
      actionButton(ns(paste0("one_", s$id)), "Run only this", class = "btn-default btn-sm")))))
  # one choice in "What would you like?": a checkbox + (rendered by the server) why it is or is not available for this design
  opt_row <- function(id, label, value = FALSE) div(class = "opt-row",
    checkboxInput(ns(paste0("opt_", id)), tags$b(label), value, width = "100%"), uiOutput(ns(paste0("optnote_", id))))

  fluidRow(
    column(6,
      div(style = "margin-bottom:10px", actionButton(ns("clear_all"), "Clear all fields / start a new project", class = "btn-default btn-sm"),
          span(class = "hint", style = "margin-left:8px", "Empties the form. Files on disk are not touched.")),
      div(class = "card", h4("1 · Your project and output folder"),
        textInput(ns("name"), "Project name", "", width = "100%", placeholder = "e.g. Heat shock in Arabidopsis roots"),
        pathUI(ns("dir"), "Output folder (where StRym saves all results)", type = "dir",
               help = "This is where your RESULTS go: StRym creates numbered sub-folders here. You choose your input data in the next section. Your input files are never changed.")),
      div(class = "card", h4("2 · Your experiment"),
        radioButtons(ns("type"), "Type of RNA-seq data",
                     c("Whole organ/tissue RNA-seq" = "bulk", "Single cell/single nucleus RNA-seq" = "single_cell")),
        conditionalPanel("input.type == 'bulk'", ns = ns,
          selectInput(ns("kind_bulk"), "What do you have?",
                      c("A counts table (genes × samples)" = "counts", "Raw reads (FASTQ files) — quantify them with the TUX module first" = "fastq"), width = "100%")),
        conditionalPanel("input.type == 'single_cell'", ns = ns,
          selectInput(ns("kind_sc"), "What do you have?",
                      c("A 10x folder or .h5 file" = "10x", "A Seurat object (.rds)" = "seurat"), width = "100%")),
        pathUI(ns("data"), "Input data", type = "any", multiple = TRUE,
               help = "A counts table (CSV / TSV / Excel, genes in rows), a single-cell file or 10x folder. For raw reads, select all your FASTQ files at once (Ctrl- or Shift-click in the file window, or drag and drop them) — or choose the folder that holds them."),
        conditionalPanel("input.type == 'bulk' && input.kind_bulk == 'counts'", ns = ns,
          selectInput(ns("unit"), "What are the numbers in the table?",
                      c("Work it out (whole numbers = raw counts, decimals = already normalized)" = "auto", "Raw counts" = "counts",
                        "Already normalized (FPKM / TPM): use as they are" = "normalized"), width = "100%"),
          pathUI(ns("diff"), "Differential-expression results from another run (optional)", type = "any", multiple = TRUE,
                 help = "Only needed for the Shannon entropy steps when you did not start from raw reads. Cuffdiff comparison files (gene_exp.diff): one file, several, a folder of them or a .zip. They can be one per genotype or time point; StRym matches them to your conditions and time points by name and takes their q-values.")),
        conditionalPanel("input.type == 'bulk' && input.kind_bulk == 'fastq'", ns = ns,
          div(class = "refbox",
            tags$b("References"),
            p(class = "hint", "StRym runs the Tuxedo pipeline for you as the first step and turns your reads into a counts table."),
            pathUI(ns("gtf"), "Annotation (GTF / GFF)", type = "file"),
            pathUI(ns("fasta"), "Genome (FASTA, uncompressed)", type = "file"),
            actionLink(ns("use_ref"), "Use the Arabidopsis (TAIR10) reference that comes with StRym"),
            pathUI(ns("adapter"), "Adapter sequences to trim (optional)", type = "file",
                   help = "A FASTA file of the adapters of your library kit. Blank uses the list that comes with StRym (Illumina TruSeq, Nextera, MGI / BGI, small RNA and others). Only for short reads."),
            selectInput(ns("platform"), "Sequencing platform", c("Short reads (Illumina, MGI / BGI ...)" = "short", "Long reads (Oxford Nanopore, PacBio; e.g. Plasmidsaurus)" = "long"), width = "100%"),
            p(class = "hint", "Short reads are cleaned, aligned with HISAT2 and compared with Cufflinks/Cuffdiff. Long reads are aligned with minimap2 and counted with featureCounts; the q-values for the entropy steps are then calculated from the counts. If your provider already gives you a counts table (Plasmidsaurus does), use \u201cA counts table\u201d instead: it is quicker.")
            ,
            radioButtons(ns("read_type"), "Read type", c("Single-end" = "single", "Paired-end (A_/B_, _R1/_R2 or _1/_2)" = "paired"), inline = TRUE),
            numericInput(ns("threads"), "Threads", max(1L, min(4L, parallel::detectCores() - 1L)),
                         min = 1, max = max(1L, parallel::detectCores()), step = 1))),
        fluidRow(
          column(6, numericInput(ns("ncond"), "Number of conditions", 2, min = 1, step = 1)),
          column(6, numericInput(ns("ntime"), "Number of time points", 1, min = 1, step = 1))),
        textInput(ns("cnames"), "Condition names (optional)", "", width = "100%", placeholder = "control, heat"),
        textInput(ns("tnames"), "Time-point names (optional)", "", width = "100%", placeholder = "0h, 1h, 4h")),
      div(class = "card sheet-card", h4("3 · Sample sheet"),
        p(class = "hint", "StRym builds this from the details above. Check that every sample has the right condition and time point."),
        uiOutput(ns("sheet_note")),
        DTOutput(ns("sheet_tbl")),
        div(style = "margin:10px 0 4px",
            actionButton(ns("save_sheet"), "Save to the project folder", class = "btn-primary btn-sm"), " ",
            downloadButton(ns("dl_sheet"), "Download a copy", class = "btn-default btn-sm")),
        uiOutput(ns("sheet_status")),
        hr(),
        pathUI(ns("ss"), "Already have a sample sheet? Use yours instead (optional)", type = "file",
               help = "One row per sample, with a “condition” column and, for time courses, a “time_point” column. For raw reads add a “file” column with each FASTQ file name. Clear the box to let StRym build one again.")),
      div(class = "card", h4("4 · What would you like?"),
        div(class = "opt-group", "Differential expression and gene selection"),
        opt_row("deseq2", "DESeq2"),
        opt_row("shannon_entropy", "Shannon entropy"),
        opt_row("modified_shannon_entropy", "Modified Shannon entropy"),
        div(class = "opt-group", "Plots and clustering"),
        opt_row("qc", "Data overview: samples per condition, depth, sample similarity", TRUE),
        opt_row("pca", "PCA"),
        div(class = "opt-row",
          tags$b("Seurat clustering maps"),
          div(class = "opt-inline", checkboxInput(ns("opt_umap"), "UMAP", FALSE), checkboxInput(ns("opt_tsne"), "t-SNE", FALSE)),
          span(class = "hint", "Choose UMAP, t-SNE, or both."), uiOutput(ns("optnote_seurat"))),
        opt_row("trajectory", "Trajectories"),
        div(class = "opt-group", "Gene regulatory network (GRN) inference"),
        opt_row("genist", "GENIST: network from a time course"),
        opt_row("rtpstar", "RTP-STAR: directed network (random forests)"),
        div(class = "opt-group", "Share"),
        checkboxInput(ns("browser"), "Build an interactive web browser to explore the results", FALSE),
        div(class = "hint", "Saved as a single .html file in the project folder that you can open, export and share.")),
      div(class = "card", h4("Analysis settings (optional)"),
        tags$details(
          tags$summary("Control condition, cutoffs, samples to leave out, gene list, single cell options"),
          p(class = "hint", style = "margin-top:8px", "The defaults are those of the original analysis scripts. Change them only if you need to."),
          conditionalPanel("input.type == 'bulk'", ns = ns,
            div(class = "opt-group", "Time-course analyses (entropy, PCA, trajectories)"),
            textInput(ns("p_control"), "Control condition(s)", "", width = "100%", placeholder = "leave blank to detect (names like control, mock, untreated, WT)"),
            textInput(ns("p_pairs"), "Treated:control pairs (advanced)", "", width = "100%", placeholder = "e.g. Heat:Control, Heat2:Control2"),
            checkboxInput(ns("p_baseline"), "The first time point is a baseline (before the treatment): shown, but never called", TRUE),
            fluidRow(
              column(6, numericInput(ns("p_log2fc"), "Minimum |log2 fold change|", 1, min = 0, step = 0.1)),
              column(6, numericInput(ns("p_q"), "Maximum q-value", 0.05, min = 0, max = 1, step = 0.01))),
            fluidRow(
              column(6, numericInput(ns("p_minfpkm"), "Minimum mean expression in some sample (CPM for counts, FPKM otherwise)", 2, min = 0, step = 0.5)),
              column(6, numericInput(ns("p_mintp"), "Main set: changed at least this many time points", 2, min = 1, step = 1))),
            textInput(ns("p_drop"), "Samples to leave out (e.g. outliers seen in the PCA)", "", width = "100%", placeholder = "sample names, comma separated"),
            pathUI(ns("trajgenes"), "Genes for the trajectory plots (optional)", type = "file",
                   help = "A text or CSV file with one gene id or name per line. Without it, the trajectories show the genes found by the entropy steps."),
            conditionalPanel("input.kind_bulk == 'fastq'", ns = ns,
              div(class = "opt-group", "Raw reads"),
              checkboxInput(ns("p_tux_fresh"), "Start the Tux step from scratch (ignore the finished parts of an earlier run)", FALSE)),
            div(class = "opt-group", "Regulatory network (GENIST)"),
            numericInput(ns("p_genist_max"), "Most genes in the network (0 = no limit)", 500, min = 0, step = 50),
            pathUI(ns("genist_tf"), "Transcription factor list (optional)", type = "file",
                   help = "A text or CSV file with one gene id per line (first column). Without it every gene may regulate every other gene, which is slow for more than about 100 genes: a list of transcription factors makes the network run many times faster.")),
          conditionalPanel("input.type == 'single_cell'", ns = ns,
            div(class = "opt-group", "Single cell / single nucleus (Seurat)"),
            fluidRow(
              column(6, numericInput(ns("p_sc_minfeat"), "Minimum genes per cell", 200, min = 0, step = 50)),
              column(6, numericInput(ns("p_sc_mincells"), "Minimum cells per gene", 3, min = 0, step = 1))),
            fluidRow(
              column(6, numericInput(ns("p_sc_mincounts"), "Minimum UMIs per cell", 500, min = 0, step = 100)),
              column(6, numericInput(ns("p_sc_maxcounts"), "Maximum UMIs per cell (0 = no limit)", 0, min = 0, step = 1000))),
            fluidRow(
              column(6, numericInput(ns("p_sc_res"), "Clustering resolution", 1, min = 0.1, step = 0.1)),
              column(6, numericInput(ns("p_sc_npcs"), "Principal components (0 = automatic)", 0, min = 0, step = 1))),
            selectInput(ns("p_sc_integ"), "Integrate the samples with", c("CCA (Seurat default)" = "cca", "RPCA (faster for big data)" = "rpca", "Harmony" = "harmony", "No integration (only merge)" = "none"), width = "100%"),
            p(class = "hint", "The cell filters depend on your protocol: look at the QC plot that StRym saves and rerun if needed.")))),
      div(class = "card", h4("Try a sample experiment"),
        actionButton(ns("demo"), "Download the sample experiment and fill everything in"),
        p(class = "hint", "Saves a small made-up heat-shock study (control vs heat, 3 time points, 2 replicates) in a ", code("sample_experiment"),
          " folder next to StRym, then fills in every field and file for you."))),
    column(6,
      div(class = "card", h4("Checks"), uiOutput(ns("validation"))),
      div(class = "card", h4("5 · Analysis plan"),
        plan_rows, uiOutput(ns("grn_plan"))),
      div(class = "card", h4("6 · Run"), job_panel_ui(ns, "Run the selected steps")),
      div(class = "card", h4("Results"), uiOutput(ns("results")),
        uiOutput(ns("browser_ui")))))
}

design_server <- function(id, roots, app_dir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns; root <- session$rootScope()
    steps <- strym_steps(); ids <- vapply(steps, `[[`, "", "id"); mods <- strym_modules()
    proj <- pathServer("dir", roots, "dir")
    data <- pathServer("data", roots, "any", c("csv", "tsv", "txt", "xlsx", "xls", "rds", "h5", "h5ad", "loom", "gz", "fastq", "fq"), multiple = TRUE)
    ss <- pathServer("ss", roots, "file")
    dif <- pathServer("diff", roots, "any", c("diff", "txt", "tsv", "csv", "zip"), multiple = TRUE)
    gtf <-pathServer("gtf", roots, "file", c("gtf", "gff", "gff3")); fasta <- pathServer("fasta", roots, "file", c("fa", "fasta", "fna"))
    observeEvent(input$use_ref, { gtf$set(file.path(app_dir, "reference", "genes.gtf")); fasta$set(file.path(app_dir, "reference", "genome.fa")) })

    chosen <- function(opts) Filter(function(i) isTRUE(input[[paste0("opt_", i)]]), opts)     # which of these options are ticked
    traj <- pathServer("trajgenes", roots, "file", c("txt", "csv", "tsv", "xlsx", "xls"))
    adapter <- pathServer("adapter", roots, "file", c("fa", "fasta", "fna", "txt"))
    gtf_list <- pathServer("genist_tf", roots, "file", c("txt", "csv", "tsv", "xlsx", "xls"))
    kind <- reactive(if (identical(input$type, "single_cell")) input$kind_sc else input$kind_bulk)
    params_in <- reactive(list(control = input$p_control, pairs = input$p_pairs, baseline_first = isTRUE(input$p_baseline),
      min_log2fc = input$p_log2fc, q_cutoff = input$p_q, min_fpkm = input$p_minfpkm, min_sig_tp = input$p_mintp, drop_reps = input$p_drop,
      genist_max_genes = input$p_genist_max, genist_tf_file = gtf_list$path(), tux_fresh = isTRUE(input$p_tux_fresh),
      trajectory_genes = traj$path(), sc_min_features = input$p_sc_minfeat, sc_min_cells = input$p_sc_mincells, sc_min_counts = input$p_sc_mincounts,
      sc_max_counts = input$p_sc_maxcounts, sc_resolution = input$p_sc_res, sc_npcs = input$p_sc_npcs, sc_integration = input$p_sc_integ))
    dsg_base <- reactive(strym_normalize_design(list(                          # the design without its sample sheet
      project_name = input$name, project_dir = proj$path(), data_type = input$type,
      data_kind = kind(),
      gtf_file = gtf$path(), fasta_file = fasta$path(), read_type = input$read_type, read_platform = if (identical(kind(), "fastq") && !is.null(input$platform)) input$platform else "short", adapter_file = adapter$path(), threads = whole_threads(input$threads),
      data_path = data$path(), fastq_files = if (identical(kind(), "fastq")) data$files() else character(),
      expr_unit = if (is.null(input$unit)) "auto" else input$unit,
      diff_path = if (!identical(kind(), "counts")) character() else if (length(dif$files())) dif$files() else dif$path(),
      sample_sheet = "", n_conditions = input$ncond, condition_names = input$cnames,
      n_timepoints = input$ntime, timepoint_names = input$tnames,
      de_methods = chosen(c("deseq2", "shannon_entropy", "modified_shannon_entropy")),
      plots = chosen(c("qc", "pca", "umap", "tsne", "trajectory")), grn_methods = chosen(c("genist", "rtpstar")),
      interactive_browser = input$browser, params = params_in())))

    # ---- the sample sheet: built from the inputs, shown for confirmation, saved in the project folder ----
    safe_name <- reactive({ n <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9._-]+", "_", trimws(input$name))); if (nzchar(n)) n else "project" })
    rv_sheet <- reactiveVal(NULL); sheet_note <- reactiveVal(""); last_key <- reactiveVal(""); pending_sheet <- reactiveVal(NULL); saved_sheet <- reactiveVal(NULL)
    own_sheet <- reactive(ss$path())                                           # a sheet the researcher brought along
    sheet_key <- reactive({ d <- dsg_base()
      as.character(jsonlite::toJSON(list(d$data_kind, d$data_path, d$fastq_files, d$read_type, d$condition_names, d$timepoint_names, d$n_conditions, d$n_timepoints), auto_unbox = TRUE)) })
    observe({                                                                  # rebuild only when something the sheet depends on really changed
      key <- sheet_key(); if (identical(key, isolate(last_key()))) return()
      last_key(key); pend <- isolate(pending_sheet())
      if (is.data.frame(pend)) { pending_sheet(NULL); rv_sheet(pend); sheet_note("Restored from the sample sheet saved with this project."); return() }
      a <- strym_auto_sample_sheet(isolate(dsg_base())); rv_sheet(a$sheet); sheet_note(a$note)
    })
    shown_sheet <- reactive({                                                  # what the page shows and the run uses
      if (nzchar(own_sheet()) && file.exists(own_sheet())) strym_read_sample_sheet(own_sheet()) else rv_sheet() })
    sheet_path <- reactive({                                                   # file the checks read (a preview copy until the run starts)
      if (nzchar(own_sheet())) return(own_sheet())
      s <- rv_sheet(); if (!is.data.frame(s)) return("")
      f <- file.path(JOB_DIR, "sheet_preview.csv"); utils::write.csv(s, f, row.names = FALSE); f })
    dsg <- reactive({ d <- dsg_base(); d$sample_sheet <- sheet_path(); d })
    plan <- reactive(strym_plan(dsg()))

    # ---- keep the step checkboxes in line with the design (user toggles survive unrelated edits) ----
    prev <- reactiveVal(NULL)
    observeEvent(plan(), {
      p <- plan(); pv <- prev()
      for (i in seq_len(nrow(p))) {
        id <- p$id[i]
        changed <- is.null(pv) || !identical(pv$default[pv$id == id], p$default[i]) || !identical(pv$applicable[pv$id == id], p$applicable[i])
        if (changed) updateCheckboxInput(session, paste0("sel_", id), value = p$default[i])
        session$sendCustomMessage("tux-toggle", list(id = ns(paste0("sel_", id)), enable = p$applicable[i]))
      }
      prev(p)
    })
    selected <- reactive({
      p <- plan(); ids[vapply(ids, function(i) isTRUE(input[[paste0("sel_", i)]]) && isTRUE(p$applicable[p$id == i]), logical(1))]
    })
    blocked <- reactive(strym_blocked(dsg(), selected()))

    for (s in steps) local({
      id <- s$id
      output[[paste0("badge_", id)]] <- renderUI({
        p <- plan(); a <- isTRUE(p$applicable[p$id == id])
        if (!a) span(class = "badge-s badge-na", "not applicable")
        else if (id %in% selected()) { b <- blocked()[[id]]; if (nzchar(b)) span(class = "badge-s badge-blocked", "blocked") else span(class = "badge-s badge-ready", "ready") }
      })
      output[[paste0("note_", id)]] <- renderUI({
        p <- plan(); note <- p$note[p$id == id]
        b <- if (id %in% selected()) blocked()[[id]] else ""
        tagList(if (nzchar(b)) div(class = "warn", style = "font-size:12px", "⚠ ", b),
                if (nzchar(note) && !nzchar(b)) div(class = "hint", style = "font-style:italic", note))
      })
    })
    mplan <- reactive(strym_module_plan(dsg()))
    for (m in mods) local({
      id <- m$id; tab <- m$tab
      observeEvent(input[[paste0("open_", id)]], updateNavbarPage(root, "module", selected = tab))
    })

    # ---- "What would you like?": every option says whether it fits this design, and why not when it is grayed out ----
    opt_map <- c(deseq2 = "deseq2", shannon_entropy = "shannon_entropy", modified_shannon_entropy = "modified_shannon_entropy",
                 qc = "qc_overview", pca = "pca", trajectory = "trajectory", umap = "seurat_clustering", tsne = "seurat_clustering",
                 genist = "genist_network", rtpstar = "rtpstar")               # option id -> step / module id
    opt_info <- function(id) {
      p <- plan(); m <- mplan(); r <- if (id %in% p$id) p[p$id == id, ] else m[m$id == id, ]
      list(app = isTRUE(r$applicable), note = r$note, fit = isTRUE(r$fit), why = r$why)
    }
    opt_note <- function(i) if (!i$app) div(class = "why-na", "Not available for this design: ", i$note)   # options that work just are not grayed out
    for (o in setdiff(names(opt_map), c("umap", "tsne"))) local({ k <- opt_map[[o]]; output[[paste0("optnote_", o)]] <- renderUI(opt_note(opt_info(k))) })
    output$optnote_seurat <- renderUI(opt_note(opt_info("seurat_clustering")))
    observe({                                                                  # gray out what cannot work, and untick it
      p <- plan(); m <- mplan(); app <- c(setNames(p$applicable, p$id), setNames(m$applicable, m$id))
      for (o in names(opt_map)) {
        a <- isTRUE(app[[opt_map[[o]]]])
        session$sendCustomMessage("tux-toggle", list(id = ns(paste0("opt_", o)), enable = a))
        if (!a && isTRUE(isolate(input[[paste0("opt_", o)]]))) updateCheckboxInput(session, paste0("opt_", o), value = FALSE)
      }
    })
    grn_names <- c(genist = "GENIST", rtpstar = "RTP-STAR")
    output$grn_plan <- renderUI({                                              # chosen network methods appear in the plan
      g <- intersect(dsg()$grn_methods, "rtpstar"); m <- mplan(); if (!length(g)) return(NULL)     # GENIST is an automatic step now (see the plan above)
      tagList(div(class = "plan-group", "Network inference (set up in its own tab)"),
        lapply(g, function(id) { r <- m[m$id == id, ]
          div(class = "plan-row",
            div(class = "plan-main", tags$b(r$label), div(class = "hint", r$desc),
                div(class = "hint", style = "font-style:italic", "Opens its own tab: it needs the gene list and expression files you choose there. ",
                    "Save its results in the ", span(class = "path", "06_regulatory_networks"), " folder inside your output folder.")),
            actionButton(ns(paste0("open_", id)), paste("Open the", grn_names[[id]], "tab"), class = "btn-default btn-sm"))
        }))
    })

    output$sheet_note <- renderUI({
      if (nzchar(own_sheet())) div(class = "hint", "Using your own sample sheet (shown below, read-only).")
      else if (!is.data.frame(rv_sheet())) div(class = "hint", sheet_note())
      else div(class = "ok", style = "font-weight:400", sheet_note())
    })
    output$sheet_tbl <- renderDT({
      s <- shown_sheet(); shiny::validate(shiny::need(is.data.frame(s) && nrow(s) > 0, "No sample sheet yet."))
      own <- nzchar(own_sheet())
      datatable(s, rownames = FALSE, selection = "none",
                editable = if (own) FALSE else list(target = "cell", disable = list(columns = c(0, which(names(s) == "file") - 1L))),
                options = list(dom = "t", paging = FALSE, scrollY = "280px", scrollCollapse = TRUE, ordering = FALSE))
    })
    observeEvent(input$sheet_tbl_cell_edit, {
      info <- input$sheet_tbl_cell_edit; s <- rv_sheet(); req(is.data.frame(s), !nzchar(own_sheet()))
      r <- info$row; ci <- info$col + 1L; v <- trimws(as.character(info$value))
      if (identical(names(s)[ci], "replicate")) {
        n <- suppressWarnings(as.integer(v)); if (is.na(n)) return(showNotification("Replicate must be a whole number.", type = "warning"))
        s[r, ci] <- n
      } else s[r, ci] <- v
      rv_sheet(s)
    })
    save_sheet_file <- function(s) {                                           # write the sheet into the project folder; returns the path
      p <- proj$path(); dir.create(p, recursive = TRUE, showWarnings = FALSE)
      f <- file.path(p, paste0("sample_sheet_", safe_name(), ".csv")); utils::write.csv(s, f, row.names = FALSE); f
    }
    observeEvent(input$save_sheet, {
      if (nzchar(own_sheet())) return(showNotification("You are using your own sample sheet, so there is nothing to save.", type = "message"))
      s <- rv_sheet(); if (!is.data.frame(s)) return(showNotification("There is no sample sheet yet: fill in your input data above first.", type = "warning"))
      if (!nzchar(proj$path())) return(showNotification("Choose the output folder first (section 1): the sample sheet is saved there.", type = "warning"))
      f <- tryCatch(save_sheet_file(s), error = function(e) NULL)
      if (is.null(f)) showNotification("Could not save into that folder (is the drive connected and writable?).", type = "error")
      else { saved_sheet(f); showNotification(paste("Saved:", to_windows_path(f)), duration = 8) }
    })
    output$sheet_status <- renderUI({ f <- saved_sheet()
      if (!is.null(f) && file.exists(f)) div(class = "hint", "Saved in the project folder: ", span(class = "path", to_windows_path(f))) })
    output$dl_sheet <- downloadHandler(function() paste0("sample_sheet_", safe_name(), ".csv"),
      function(file) { s <- shown_sheet(); utils::write.csv(if (is.data.frame(s)) s else data.frame(), file, row.names = FALSE) })

    # ---- checks ----
    check <- reactive({
      safe <- function(expr, fallback) tryCatch(expr, error = function(e) fallback(conditionMessage(e)))      # an error in a check must never end the session (grey page)
      v <- safe(strym_validate_design(dsg()), function(m) list(errors = paste("The design could not be checked:", m), warnings = character(), info = character()))
      sel <- safe(strym_check_selection(dsg(), selected()), function(m) list(errors = character(), warnings = paste("The selection could not be checked:", m)))
      errs <- c(v$errors, sel$errors); warns <- c(v$warnings, sel$warnings, strym_check_networks(dsg(), selected()))
      if (any(c("shannon_entropy", "modified_shannon_entropy", "trajectory") %in% selected())) {           # the treated-vs-control comparisons, before the run
        cc <- safe(strym_check_comparisons(shown_sheet(), dsg()$params, dsg()$timepoint_names), function(m) list(problems = paste("The comparisons could not be checked:", m), pairs = NULL))
        warns <- c(warns, cc$problems)
        if (!is.null(cc$pairs)) v$info <- c(v$info, paste0("Comparisons: ", paste(paste0(cc$pairs$treated, " vs ", cc$pairs$control), collapse = "; "), "."))
      }
      for (p in c(proj$path(), data$path(), ss$path(), gtf$path(), fasta$path())) if (nzchar(p)) { u <- unmounted_drive_problem(p); if (nzchar(u)) errs <- c(errs, u) }
      list(errors = unique(errs), warnings = unique(warns), info = v$info)
    })
    output$validation <- renderUI({
      v <- check()
      if (!length(v$errors) && !length(v$warnings))
        return(tagList(div(class = "ok", "✔ Everything looks good."), lapply(v$info, function(i) div(class = "hint", i))))
      tagList(lapply(v$errors, function(e) div(class = "bad", "✘ ", e)), lapply(v$warnings, function(w) div(class = "warn", "⚠ ", w)),
              lapply(v$info, function(i) div(class = "hint", i)))
    })
    can_run <- reactive(!length(check()$errors))
    design_ok <- reactive(!length(strym_validate_design(dsg())$errors))

    # ---- run ----
    tick <- reactiveVal(0)
    job <- job_controller("analysis", "run_analysis_cli.R", input, output, session, can_run, app_dir,
                          on_finished = function(st) tick(tick() + 1))
    start <- function(steps_to_run) {
      d <- dsg()
      if (!nzchar(own_sheet()) && is.data.frame(rv_sheet())) {                 # the sheet StRym built goes into the project folder and is used from there
        d$sample_sheet <- save_sheet_file(rv_sheet()); saved_sheet(d$sample_sheet)
      }
      job$start(list(design = d, steps = steps_to_run), d$project_dir, file.path(d$project_dir, "index.html"))
    }
    pending <- reactiveVal(character())
    ask_run <- function(steps_to_run) {                                        # nothing starts until the researcher confirms
      d <- dsg(); pending(steps_to_run)
      est <- tryCatch(strym_estimate_run(d, steps_to_run), error = function(e) list(seconds = 0, need = NA_real_, free = NA_real_, lines = character()))
      sh <- shown_sheet(); win <- function(p) if (nzchar(p)) to_windows_path(p) else "(none)"
      row <- function(k, v) tags$tr(tags$th(k, style = "width:150px"), tags$td(v))
      kind_txt <- c(counts = "a counts table", fastq = "raw reads (FASTQ)", `10x` = "a 10x folder / .h5 file", seurat = "a Seurat object")[[d$data_kind]]
      facts <- tags$table(class = "table table-condensed tbl-sm",
        row("Project", d$project_name), row("Output folder", span(class = "path", win(d$project_dir))),
        row("Input", tagList(span(class = "path", win(d$data_path)), if (length(d$fastq_files)) sprintf(" (%d chosen files)", length(d$fastq_files)), paste0(" — ", kind_txt))),
        if (length(d$diff_path)) row("Comparisons", tagList(span(class = "path", win(d$diff_path[1])), if (length(d$diff_path) > 1) sprintf(" (+%d more)", length(d$diff_path) - 1L), " (Cuffdiff results from another run)")),
        row("Data type", if (d$data_type == "bulk") "Whole organ/tissue RNA-seq" else "Single cell/single nucleus RNA-seq"),
        row("Conditions", paste0(d$n_conditions, if (length(d$condition_names)) paste0(" (", paste(d$condition_names, collapse = ", "), ")") else "")),
        row("Time points", paste0(d$n_timepoints, if (length(d$timepoint_names)) paste0(" (", paste(d$timepoint_names, collapse = ", "), ")") else "")),
        if (any(c("shannon_entropy", "modified_shannon_entropy", "trajectory") %in% steps_to_run)) {
          cc <- strym_check_comparisons(sh, d$params, d$timepoint_names)
          if (!is.null(cc$pairs)) row("Compared", paste(paste0(cc$pairs$treated, " vs ", cc$pairs$control, " (control)"), collapse = "; "))
        })
      labs <- vapply(steps, `[[`, "", "label")[match(steps_to_run, ids)]
      space <- if (is.finite(est$need)) {
        low <- is.finite(est$free) && est$free < est$need
        div(class = if (low) "bad" else "ok", style = "font-weight:400", sprintf("About %s needed on the output drive; %s free.", fmt_bytes(est$need), fmt_bytes(est$free)))
      } else div(class = "hint", "Small: a few MB. Your reads are not copied or changed.")
      tm <- if (est$seconds >= 60) tagList(
        div(HTML(sprintf("Estimated time: <b>%s</b>. Allow up to <b>%s</b> (the estimate + 10%%).", strym_fmt_duration(est$seconds), strym_fmt_duration(est$seconds * 1.1, up = TRUE)))),
        lapply(est$lines, function(l) div(class = "hint", l)),
        div(class = "hint", "Based on a previous run on this computer; it varies with the number of reads and with other programs running."))
        else div("Estimated time: under a minute.")
      sheet_ui <- if (is.data.frame(sh) && nrow(sh)) {
        tp <- if ("time_point" %in% names(sh) && any(nzchar(sh$time_point))) sh$time_point else rep("(one time point)", nrow(sh))
        tab <- table(Condition = sh$condition, `Time point` = tp)
        shown <- utils::head(sh, 40)
        tagList(div(sprintf("%d samples, %d condition%s, %d time point%s. Is this what you expect?", nrow(sh), length(unique(sh$condition)), if (length(unique(sh$condition)) == 1) "" else "s",
                            length(unique(tp)), if (length(unique(tp)) == 1) "" else "s")),
          tags$table(class = "table table-condensed tbl-sm", style = "width:auto",
            tags$thead(tags$tr(tags$th("Condition"), lapply(colnames(tab), tags$th))),
            tags$tbody(lapply(rownames(tab), function(r) tags$tr(tags$td(tags$b(r)), lapply(tab[r, ], function(n) tags$td(paste0(n, " sample", if (n == 1) "" else "s"))))))),
          tags$details(tags$summary("Show every sample"), tags$table(class = "table table-condensed tbl-sm",
            tags$thead(tags$tr(lapply(names(shown), tags$th))), tags$tbody(lapply(seq_len(nrow(shown)), function(i) tags$tr(lapply(shown[i, ], function(v) tags$td(as.character(v))))))),
            if (nrow(sh) > 40) div(class = "hint", sprintf("… and %d more.", nrow(sh) - 40))))
      } else div(class = "hint", "No sample sheet: samples are grouped as described in the checks.")
      warns <- check()$warnings
      showModal(modalDialog(title = "Ready to start? Please check", size = "l", easyClose = FALSE,
        h5("Your experiment"), facts,
        h5("What will run"), tags$ol(lapply(labs, tags$li)),
        h5("Disk space and time"), space, tm,
        h5("Sample sheet (this is how StRym understood your experiment)"), sheet_ui,
        if (length(warns)) tagList(h5("Things to be aware of"), lapply(warns, function(w) div(class = "warn", "⚠ ", w))),
        footer = tagList(modalButton("Go back and edit"), actionButton(ns("run_confirm"), "Yes, start the run", class = "btn-primary"))))
    }
    observeEvent(input$run, { req(can_run(), !job$running()); ask_run(selected()) })
    observeEvent(input$run_confirm, { removeModal(); req(!job$running()); st <- pending(); if (length(st)) { pending(character()); start(st) } })
    for (s in steps) local({
      id <- s$id
      observeEvent(input[[paste0("one_", id)]], { req(design_ok(), !job$running())
      why <- tryCatch(strym_blocked(dsg(), id)[[id]], error = function(e) "")
      if (nzchar(why)) return(showNotification(paste0(s$label, ": ", why), type = "error", duration = 12))        # not a run that ends "Finished" with the step skipped
      ask_run(id) })
    })
    observe({                                                  # per-step buttons follow the design and the job state
      p <- plan(); run <- job$running(); ok <- design_ok()
      for (i in seq_len(nrow(p))) session$sendCustomMessage("tux-toggle", list(id = ns(paste0("one_", p$id[i])), enable = ok && p$applicable[i] && !run))
    })

    # ---- results ----
    summary_df <- reactive({
      tick(); f <- file.path(proj$path(), "run_summary.csv")
      if (nzchar(proj$path()) && file.exists(f)) tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    })
    observe({ tick(); p <- sub("/+$", "", proj$path()); if (nzchar(p) && dir.exists(p)) tryCatch(addResourcePath("strymproj", p), error = function(e) NULL) })   # lets the page show plot previews (never ends the session)
    output$results <- renderUI({
      s <- summary_df()
      if (is.null(s)) return(div(class = "hint", "No results in the project folder yet. Describe your experiment above and press Run."))
      pd <- proj$path()
      rows <- lapply(seq_len(nrow(s)), function(i) tags$tr(
        tags$td(tags$b(s$label[i]), br(), span(class = "path", paste0(s$folder[i], "/"))),
        tags$td(span(class = paste0("badge-s badge-", s$status[i]), s$status[i])),
        tags$td(class = "hint", s$message[i])))
      imgs <- list.files(pd, pattern = "\\.(png|svg|jpe?g)$", recursive = TRUE, ignore.case = TRUE)
      imgs <- imgs[!grepl("^07_interactive_browser/", imgs)]
      tagList(
        div(class = "path", "Project folder: ", to_windows_path(pd)),
        tags$a(href = "strymproj/index.html", target = "_blank", "Open the summary page ↗"), br(), br(),
        tags$table(class = "table table-condensed", tags$thead(tags$tr(tags$th("Step"), tags$th("Status"), tags$th("Note"))), tags$tbody(rows)),
        if (length(imgs)) tagList(h5("Plots"), div(class = "gallery", lapply(imgs, function(f) tags$figure(
          tags$a(href = paste0("strymproj/", URLencode(f)), target = "_blank", tags$img(src = paste0("strymproj/", URLencode(f)))), tags$figcaption(f))))))
    })

    # ---- interactive browser launcher ----
    bstate_file <- file.path(JOB_DIR, "browser_state.json")
    btick <- reactiveTimer(2000); bversion <- reactiveVal(0)
    browser_alive <- reactive({ btick(); bversion()
      st <- if (file.exists(bstate_file)) tryCatch(jsonlite::fromJSON(bstate_file), error = function(e) NULL)
      !is.null(st) && pid_alive_pattern(st$pid, "runApp") })
    browser_dir <- reactive({ tick(); d <- file.path(proj$path(), "07_interactive_browser"); if (file.exists(file.path(d, "app.R"))) d })
    output$browser_ui <- renderUI({
      tick(); d <- browser_dir(); if (is.null(d)) return(NULL)
      alive <- browser_alive()
      html_f <- file.path(d, "strym_browser.html")
      tagList(hr(), h5("Interactive browser"),
        if (file.exists(html_f)) div(
          p(class = "hint", "A single .html file with your plots and tables: open it, save it or e-mail it. It works without StRym or an internet connection."),
          tags$a(href = "strymproj/07_interactive_browser/strym_browser.html", target = "_blank", class = "btn btn-primary btn-sm", "Open the browser (.html) ↗"), " ",
          tags$a(href = "strymproj/07_interactive_browser/strym_browser.html", download = "strym_browser.html", class = "btn btn-default btn-sm", "Download the .html"),
          div(class = "path", style = "margin-top:6px", to_windows_path(html_f)), br()),
        if (alive) tagList(div(class = "ok", "Running — ", tags$a(href = "http://localhost:3939", target = "_blank", "open http://localhost:3939 ↗")),
                           actionButton(ns("stop_browser"), "Stop the browser", class = "btn-default btn-sm"))
        else tagList(div(class = "hint", "Optional live version (a template that reads the project folder while StRym is running)."),
                     actionButton(ns("launch_browser"), "Launch the live browser", class = "btn-default btn-sm")))
    })
    observeEvent(input$launch_browser, {
      d <- browser_dir(); req(d)
      if (browser_alive()) return()
      code <- sprintf("shiny::runApp(%s, port = 3939, host = \"0.0.0.0\", launch.browser = FALSE)", deparse(normalizePath(d)))
      p <- processx::process$new("Rscript", c("-e", code), stdout = file.path(JOB_DIR, "browser.log"), stderr = "2>&1",
                                 cleanup = FALSE, cleanup_tree = FALSE)
      jsonlite::write_json(list(pid = p$get_pid(), dir = d), bstate_file, auto_unbox = TRUE); bversion(bversion() + 1)
    })
    observeEvent(input$stop_browser, {
      st <- if (file.exists(bstate_file)) tryCatch(jsonlite::fromJSON(bstate_file), error = function(e) NULL)
      if (!is.null(st) && pid_alive_pattern(st$pid, "runApp")) kill_tree(st$pid); bversion(bversion() + 1)
    })

    # ---- demo, and restoring the last design when the page is reopened ----
    apply_design <- function(d) {
      updateTextInput(session, "name", value = d$project_name); proj$set(d$project_dir)
      updateRadioButtons(session, "type", selected = d$data_type)
      if (identical(d$data_type, "single_cell")) updateSelectInput(session, "kind_sc", selected = d$data_kind)
      else updateSelectInput(session, "kind_bulk", selected = d$data_kind)
      data$set(d$data_path); data$set_files(d$fastq_files)
      updateSelectInput(session, "unit", selected = d$expr_unit)
      if (length(d$diff_path) > 1) dif$set_files(d$diff_path) else dif$set(if (length(d$diff_path)) d$diff_path else "")
      auto_saved <- nzchar(d$sample_sheet) && file.exists(d$sample_sheet) && grepl("^sample_sheet_.*\\.csv$", basename(d$sample_sheet)) &&
        identical(normalizePath(dirname(d$sample_sheet), mustWork = FALSE), normalizePath(d$project_dir, mustWork = FALSE))
      if (auto_saved) {                                                        # the sheet StRym built earlier: show it again (with its edits)
        pending_sheet(strym_read_sample_sheet(d$sample_sheet)); ss$set(""); saved_sheet(d$sample_sheet)
      } else { pending_sheet(NULL); ss$set(d$sample_sheet) }
      gtf$set(d$gtf_file); fasta$set(d$fasta_file)
      updateRadioButtons(session, "read_type", selected = d$read_type); updateNumericInput(session, "threads", value = d$threads)
      updateSelectInput(session, "platform", selected = d$read_platform); adapter$set(d$adapter_file)
      updateNumericInput(session, "ncond", value = d$n_conditions); updateNumericInput(session, "ntime", value = d$n_timepoints)
      updateTextInput(session, "cnames", value = paste(d$condition_names, collapse = ", "))
      updateTextInput(session, "tnames", value = paste(d$timepoint_names, collapse = ", "))
      for (o in names(opt_map)) {                                              # tick the options the design asks for
        want <- o %in% c(d$de_methods, d$plots, d$grn_methods)
        updateCheckboxInput(session, paste0("opt_", o), value = want)
      }
      updateCheckboxInput(session, "browser", value = d$interactive_browser)
      pp <- strym_normalize_params(d$params)                                  # analysis settings
      updateTextInput(session, "p_control", value = pp$control); updateTextInput(session, "p_pairs", value = pp$pairs); updateTextInput(session, "p_drop", value = pp$drop_reps)
      updateCheckboxInput(session, "p_baseline", value = pp$baseline_first)
      for (k in c("p_log2fc", "p_q", "p_minfpkm", "p_mintp", "p_sc_minfeat", "p_sc_mincells", "p_sc_mincounts", "p_sc_maxcounts", "p_sc_res", "p_sc_npcs"))
        updateNumericInput(session, k, value = pp[[c(p_log2fc = "min_log2fc", p_q = "q_cutoff", p_minfpkm = "min_fpkm", p_mintp = "min_sig_tp", p_sc_minfeat = "sc_min_features",
                                                    p_sc_mincells = "sc_min_cells", p_sc_mincounts = "sc_min_counts", p_sc_maxcounts = "sc_max_counts", p_sc_res = "sc_resolution", p_sc_npcs = "sc_npcs")[[k]]]])
      updateSelectInput(session, "p_sc_integ", selected = pp$sc_integration); traj$set(pp$trajectory_genes)
      updateNumericInput(session, "p_genist_max", value = pp$genist_max_genes); gtf_list$set(pp$genist_tf_file); updateCheckboxInput(session, "p_tux_fresh", value = isTRUE(pp$tux_fresh))
    }
    observeEvent(input$demo, {                                                 # download (create) the sample experiment, then fill in every field
      base <- file.path(app_dir, "sample_experiment")
      files <- tryCatch(strym_make_demo_data(file.path(base, "data")), error = function(e) NULL)
      if (is.null(files)) return(showNotification("Could not save the sample experiment next to StRym (is the folder writable?).", type = "error"))
      writeLines(c("StRym sample experiment (made-up data)", "",
                   "data/demo_counts.csv         300 genes x 12 samples of raw counts",
                   "data/demo_sample_sheet.csv   the matching sample sheet (StRym also builds this itself)", "",
                   "Design: control vs heat, 3 time points (0h, 1h, 4h), 2 replicates each.",
                   "results/ is where StRym saves what it produces when you run the sample experiment."), file.path(base, "README.txt"))
      apply_design(strym_normalize_design(list(
        project_name = "Sample experiment (heat shock)", project_dir = file.path(base, "results"),
        data_type = "bulk", data_kind = "counts", data_path = files$counts, sample_sheet = "",
        n_conditions = 2, condition_names = "control, heat", n_timepoints = 3, timepoint_names = "0h, 1h, 4h",
        de_methods = c("deseq2", "modified_shannon_entropy"), plots = c("pca", "trajectory"), interactive_browser = TRUE)))
      showNotification(sprintf("Sample experiment saved in %s and every field filled in. Check the sample sheet, then press Run.", to_windows_path(base)), duration = 12)
    })
    # ---- clear all fields / start a new project (nothing on disk is deleted) ----
    observeEvent(input$clear_all, {
      if (job$running()) return(showNotification("A run is in progress: wait for it to finish (or cancel it) before clearing.", type = "warning"))
      showModal(modalDialog(title = "Clear all fields?", easyClose = TRUE,
        "This empties the form so you can start a new project. Nothing on disk is deleted: your output folder, results and input files stay exactly as they are.",
        footer = tagList(modalButton("Keep my entries"), actionButton(ns("clear_confirm"), "Clear all fields", class = "btn-primary"))))
    })
    observeEvent(input$clear_confirm, {
      removeModal()
      if (job$running()) return()
      updateTextInput(session, "name", value = ""); proj$set(""); data$set(""); ss$set(""); gtf$set(""); fasta$set("")
      rv_sheet(NULL); pending_sheet(NULL); saved_sheet(NULL); last_key(""); sheet_note("")
      updateRadioButtons(session, "type", selected = "bulk"); updateSelectInput(session, "kind_bulk", selected = "counts"); updateSelectInput(session, "kind_sc", selected = "10x")
      updateRadioButtons(session, "read_type", selected = "single"); updateSelectInput(session, "platform", selected = "short")
      updateSelectInput(session, "unit", selected = "auto"); dif$set(""); adapter$set("")                                     # the table-unit choice and the Cuffdiff results box
      updateNumericInput(session, "threads", value = max(1L, min(4L, parallel::detectCores() - 1L)))
      updateNumericInput(session, "ncond", value = 2); updateNumericInput(session, "ntime", value = 1)
      updateTextInput(session, "cnames", value = ""); updateTextInput(session, "tnames", value = "")
      for (o in names(opt_map)) updateCheckboxInput(session, paste0("opt_", o), value = identical(o, "qc"))         # the Data overview stays ticked, as on a fresh page
      for (i in ids) updateCheckboxInput(session, paste0("sel_", i), value = FALSE)
      updateCheckboxInput(session, "browser", value = FALSE)
      for (k in c("p_control", "p_pairs", "p_drop")) updateTextInput(session, k, value = ""); traj$set(""); updateCheckboxInput(session, "p_baseline", value = TRUE)
      df <- strym_default_params(); updateNumericInput(session, "p_genist_max", value = df$genist_max_genes); gtf_list$set(""); updateCheckboxInput(session, "p_tux_fresh", value = FALSE)
      updateNumericInput(session, "p_log2fc", value = df$min_log2fc); updateNumericInput(session, "p_q", value = df$q_cutoff); updateNumericInput(session, "p_minfpkm", value = df$min_fpkm)
      updateNumericInput(session, "p_mintp", value = df$min_sig_tp); updateNumericInput(session, "p_sc_minfeat", value = df$sc_min_features); updateNumericInput(session, "p_sc_mincells", value = df$sc_min_cells)
      updateNumericInput(session, "p_sc_mincounts", value = df$sc_min_counts); updateNumericInput(session, "p_sc_maxcounts", value = df$sc_max_counts)
      updateNumericInput(session, "p_sc_res", value = df$sc_resolution); updateNumericInput(session, "p_sc_npcs", value = df$sc_npcs); updateSelectInput(session, "p_sc_integ", selected = df$sc_integration)
      prev(NULL)
      unlink(c(job_state_file("analysis"), file.path(JOB_DIR, "analysis_config.json")))   # so the page does not bring the old project back
      job$reset(); tick(tick() + 1)
      showNotification("All fields cleared. Ready for a new project.", duration = 5)
    })
    observeEvent(TRUE, {
      st <- job$state(); cf <- file.path(JOB_DIR, "analysis_config.json")
      cfg <- if (!is.null(st) && file.exists(cf)) tryCatch(jsonlite::fromJSON(cf), error = function(e) NULL)
      if (!is.null(cfg$design)) apply_design(strym_normalize_design(cfg$design))
    }, once = TRUE)
  })
}
