design <- function(...) { d <- strym_default_design(); a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]; d }
plan_row <- function(p, id) p[p$id == id, ]

# ---- the plan: which steps apply, which are ticked ----
test_that("steps are ticked from what the user asked for, if they apply", {
  p <- strym_plan(design(data_kind = "fastq", n_conditions = 2, n_timepoints = 3, de_methods = c("deseq2", "shannon_entropy"),
                         plots = c("pca"), interactive_browser = TRUE))
  expect_true(plan_row(p, "deseq2")$default); expect_true(plan_row(p, "shannon_entropy")$default)
  expect_false(plan_row(p, "modified_shannon_entropy")$default)         # not requested
  expect_true(plan_row(p, "pca")$default); expect_false(plan_row(p, "seurat_clustering")$default)
  expect_false(plan_row(p, "trajectory")$default); expect_true(plan_row(p, "web_browser")$default)
  p <- strym_plan(design(data_type = "single_cell", data_kind = "10x", plots = c("umap", "pca")))       # single cell: Seurat instead of PCA
  expect_true(plan_row(p, "seurat_clustering")$default); expect_false(plan_row(p, "pca")$applicable)
})

test_that("one condition and one time point leave nothing to compare", {
  p <- strym_plan(design(n_conditions = 1, n_timepoints = 1, de_methods = c("deseq2", "modified_shannon_entropy")))
  for (id in c("deseq2", "shannon_entropy", "modified_shannon_entropy")) expect_false(plan_row(p, id)$applicable)
  expect_false(plan_row(p, "deseq2")$default)                           # requested but impossible -> not ticked
  expect_match(plan_row(p, "deseq2")$note, "at least 2")
})

test_that("Seurat clustering is grayed out for whole organ/tissue data, with the reason", {
  p <- strym_plan(design(data_type = "bulk", n_timepoints = 3, plots = c("umap", "tsne")))
  expect_false(plan_row(p, "seurat_clustering")$applicable); expect_false(plan_row(p, "seurat_clustering")$default)
  expect_match(plan_row(p, "seurat_clustering")$note, "single cell/single nucleus")
  p <- strym_plan(design(data_type = "single_cell", plots = "tsne"))
  expect_true(plan_row(p, "seurat_clustering")$applicable); expect_true(plan_row(p, "seurat_clustering")$default)
  expect_true(plan_row(p, "seurat_clustering")$fit)
})

test_that("the entropy methods need Tux results and several time points; the reason is given", {
  fq <- function(...) design(data_kind = "fastq", n_conditions = 2, ...)
  p <- strym_plan(fq(n_timepoints = 1, de_methods = "modified_shannon_entropy"))
  expect_false(plan_row(p, "modified_shannon_entropy")$applicable); expect_false(plan_row(p, "modified_shannon_entropy")$default)
  expect_match(plan_row(p, "modified_shannon_entropy")$note, "at least 3 time points")
  expect_false(plan_row(p, "shannon_entropy")$applicable); expect_match(plan_row(p, "shannon_entropy")$note, "at least 2 time points")
  p <- strym_plan(fq(n_timepoints = 2))                                     # 2 time points: Shannon yes, modified (outliers) not yet
  expect_true(plan_row(p, "shannon_entropy")$applicable); expect_false(plan_row(p, "modified_shannon_entropy")$applicable)
  expect_true(plan_row(strym_plan(fq(n_timepoints = 3)), "modified_shannon_entropy")$applicable)
  p <- strym_plan(design(n_conditions = 2, n_timepoints = 4, expr_unit = "normalized"))   # normalized values (FPKM/TPM): no q-values can be calculated
  expect_false(plan_row(p, "shannon_entropy")$applicable); expect_match(plan_row(p, "shannon_entropy")$note, "start from raw reads")
  expect_false(plan_row(p, "modified_shannon_entropy")$applicable)
  expect_true(plan_row(strym_plan(design(n_conditions = 2, n_timepoints = 4)), "shannon_entropy")$applicable)   # raw counts: q-values are calculated with DESeq2
  p <- strym_plan(design(data_type = "single_cell", data_kind = "10x", n_conditions = 2, n_timepoints = 4))
  expect_match(plan_row(p, "shannon_entropy")$note, "single cell/single nucleus")
})

test_that("each option says whether it is a best fit for the design", {
  one <- strym_plan(design(data_type = "bulk", n_conditions = 2, n_timepoints = 1))
  expect_true(plan_row(one, "deseq2")$fit); expect_false(plan_row(one, "shannon_entropy")$fit)
  expect_false(plan_row(one, "modified_shannon_entropy")$fit)
  tc <- strym_plan(design(data_type = "bulk", data_kind = "fastq", n_conditions = 2, n_timepoints = 4))
  expect_true(plan_row(tc, "shannon_entropy")$fit); expect_true(plan_row(tc, "modified_shannon_entropy")$fit)
  expect_true(all(nzchar(tc$why[tc$id %in% c("deseq2", "shannon_entropy", "modified_shannon_entropy", "pca", "trajectory")])))
  expect_false(plan_row(strym_plan(design(n_conditions = 1)), "deseq2")$fit)   # not applicable -> never a best fit
})

test_that("with one condition there is no control to compare with: entropy methods and DESeq2 do not apply", {
  p <- strym_plan(design(data_kind = "fastq", n_conditions = 1, n_timepoints = 4))
  expect_false(plan_row(p, "shannon_entropy")$applicable); expect_match(plan_row(p, "shannon_entropy")$note, "at least 2 conditions")
  expect_false(plan_row(p, "deseq2")$applicable)
})

test_that("trajectories need 3 time points of whole organ/tissue data; PCA is for bulk data", {
  expect_false(plan_row(strym_plan(design(n_timepoints = 2, plots = "trajectory")), "trajectory")$applicable)
  expect_true(plan_row(strym_plan(design(n_timepoints = 3, plots = "trajectory")), "trajectory")$applicable)
  sc <- strym_plan(design(data_type = "single_cell", data_kind = "10x", n_timepoints = 3, plots = "trajectory"))
  expect_false(plan_row(sc, "trajectory")$applicable); expect_match(plan_row(sc, "trajectory")$note, "whole organ/tissue")
  expect_true(plan_row(strym_plan(design()), "pca")$applicable)
  expect_false(plan_row(sc, "pca")$applicable); expect_match(plan_row(sc, "pca")$note, "Seurat")
})

test_that("blocked steps are reported with the reason", {
  d <- design(n_timepoints = 4)
  b <- strym_blocked(d, c("trajectory"));            expect_match(b[["trajectory"]], "Needs one of")
  expect_equal(strym_blocked(d, c("trajectory", "shannon_entropy"))[["trajectory"]], "")
  expect_match(strym_blocked(d, c("trajectory", "seurat_clustering"))[["trajectory"]], "Needs one of")     # Seurat is no longer a source of genes
  d2 <- d; d2$params <- list(trajectory_genes = "/some/genes.txt")                                          # ... but a gene list needs no earlier step
  expect_equal(strym_blocked(d2, "trajectory")[["trajectory"]], "")
  expect_match(strym_blocked(d, "web_browser")[["web_browser"]], "at least one")
  expect_equal(strym_blocked(d, c("web_browser", "pca"))[["web_browser"]], "")
  expect_match(strym_blocked(design(n_timepoints = 1), c("trajectory", "pca"))[["trajectory"]], "3 time points")
})

test_that("selection check: nothing selected is an error, blocked steps are warnings", {
  expect_match(strym_check_selection(design(), character())$errors, "at least one")
  r <- strym_check_selection(design(n_timepoints = 4), c("trajectory", "pca"))
  expect_length(r$errors, 0); expect_match(r$warnings, "Trajectory plots")
  expect_match(strym_check_selection(design(n_timepoints = 4), "trajectory")$errors, "None of the selected")
})

test_that("network modules are offered according to the design", {
  m <- strym_module_plan(design(n_timepoints = 2, data_kind = "counts"))
  expect_false(m$applicable[m$id == "genist"]); expect_true(m$applicable[m$id == "rtpstar"]); expect_false(m$applicable[m$id == "tux"])
  m <- strym_module_plan(design(n_timepoints = 5, data_kind = "fastq"))
  expect_true(all(m$applicable))
})

test_that("network inference choices are kept in the design and need a gene-selection method", {
  d <- strym_normalize_design(design(grn_methods = c("genist", "rtpstar")))
  expect_equal(d$grn_methods, c("genist", "rtpstar")); expect_length(strym_normalize_design(design())$grn_methods, 0)
  expect_length(strym_check_networks(design(), "pca"), 0)                                        # nothing chosen: nothing to say
  w <- strym_check_networks(design(grn_methods = "genist"), "pca"); expect_match(w, "GENIST needs a list of genes")
  expect_match(strym_check_networks(design(grn_methods = c("genist", "rtpstar")), "pca"), "GENIST and RTP-STAR")
  expect_length(strym_check_networks(design(grn_methods = "genist"), c("pca", "shannon_entropy")), 0)
  m <- strym_module_plan(design(n_timepoints = 1)); expect_false(m$applicable[m$id == "genist"]); expect_match(m$note[m$id == "genist"], "at least 3")
  expect_true(m$fit[m$id == "rtpstar"]); expect_false(m$fit[m$id == "genist"])
})

# ---- the design checks ----
test_that("design validation explains what is missing", {
  v <- strym_validate_design(design(), check_paths = TRUE)
  expect_true(any(grepl("name", v$errors))); expect_true(any(grepl("project folder", v$errors))); expect_true(any(grepl("expression data", v$errors)))
  v <- strym_validate_design(design(project_name = "x", project_dir = tempfile(), n_conditions = 0, n_timepoints = 0), FALSE)
  expect_true(any(grepl("conditions must be", v$errors))); expect_true(any(grepl("time points must be", v$errors)))
})

test_that("design validation compares names and the sample sheet with the numbers given", {
  d <- tempfile("v"); dir.create(d); demo <- strym_make_demo_data(d)
  ok <- design(project_name = "p", project_dir = file.path(d, "proj"), data_path = demo$counts, sample_sheet = demo$sample_sheet,
               n_conditions = 2, n_timepoints = 3)
  v <- strym_validate_design(ok); expect_length(v$errors, 0); expect_length(v$warnings, 0)
  expect_true(any(grepl("12 samples, 2 conditions", v$info))); expect_true(any(grepl("table: 300 rows", v$info)))
  v <- strym_validate_design(modifyList(ok, list(n_conditions = 3, n_timepoints = 2, condition_names = "a, b")))
  expect_true(any(grepl("sample sheet has 2 conditions but you set 3", v$warnings)))
  expect_true(any(grepl("sample sheet has 3 time points but you set 2", v$warnings)))
  expect_true(any(grepl("listed 2 condition names", v$warnings)))
  v <- strym_validate_design(modifyList(ok, list(sample_sheet = "")))
  expect_true(any(grepl("No sample sheet", v$warnings)))
  v <- strym_validate_design(modifyList(ok, list(data_path = file.path(d, "missing.csv"))))
  expect_true(any(grepl("not found", v$errors)))
})

test_that("sample sheet template has one row per condition x time point x replicate", {
  t <- strym_sample_sheet_template(design(n_conditions = 2, condition_names = "control, heat", n_timepoints = 3), 2)
  expect_equal(nrow(t), 12); expect_equal(unique(t$condition), c("control", "heat")); expect_equal(unique(t$time_point), c("t1", "t2", "t3"))
  expect_named(t, c("sample", "condition", "time_point", "replicate", "file")); expect_false(anyDuplicated(t$sample) > 0)
})

# ---- running ----
demo_design <- function(root, ...) {
  demo <- strym_make_demo_data(file.path(root, "data"))
  design(project_name = "Heat shock demo", project_dir = file.path(root, "project"), data_path = demo$counts,
         sample_sheet = demo$sample_sheet, n_conditions = 2, condition_names = "control, heat", n_timepoints = 3,
         timepoint_names = "0h, 1h, 4h", ...)
}

test_that("a full run writes the organized project folder and its summary", {
  root <- tempfile("run"); dir.create(root)
  demo <- strym_make_demo_data(file.path(root, "data")); glist <- file.path(root, "genes.txt"); writeLines(c("GENE001", "GENE002", "GENE003"), glist)
  d <- demo_design(root, de_methods = "deseq2", plots = c("pca", "trajectory"), interactive_browser = TRUE, params = list(trajectory_genes = glist))
  res <- quiet(strym_run(d, app_dir = here))
  pd <- file.path(root, "project")
  expect_equal(res$summary$step, c("deseq2", "pca", "trajectory", "web_browser"))
  expect_equal(res$summary$status, c("ok", "ok", "ok", "ok"))                       # every step is real now
  for (f in c("00_design/design.json", "00_design/plan.csv", "00_design/sample_sheet.csv", "run_summary.csv", "index.html",
              "02_differential_expression/deseq2/deseq2_DEGs.csv", "02_differential_expression/deseq2/normalized_counts.csv", "04_pca/pca_plot.png", "04_pca/pca_scores.csv", "05_trajectories/trajectories.png",
              "05_trajectories/trajectories.pdf", "07_interactive_browser/app.R", "07_interactive_browser/README.txt", "07_interactive_browser/strym_browser.html"))
    expect_true(file.exists(file.path(pd, f)), info = f)
  expect_true(dir.exists(file.path(pd, "06_regulatory_networks")))
  back <- jsonlite::fromJSON(file.path(pd, "00_design", "design.json"))
  expect_equal(back$n_timepoints, 3); expect_equal(back$condition_names, c("control", "heat"))
  expect_no_error(parse(file.path(pd, "07_interactive_browser", "app.R")))        # the template app is valid R
  html <- paste(readLines(file.path(pd, "index.html")), collapse = "\n")
  expect_match(html, "Heat shock demo"); expect_match(html, "Trajectory plots"); expect_match(html, "DESeq2 differential expression")
  expect_equal(back$params$trajectory_genes, glist)                                                          # the settings are kept with the design
  expect_equal(read.csv(file.path(pd, "run_summary.csv"))$step, res$summary$step)
  expect_match(html, '<img alt="StRym"[^>]*src="data:image/png;base64,'); expect_gt(nchar(html), 50000)       # the logo is embedded in the results page
  expect_match(paste(readLines(file.path(pd, "07_interactive_browser", "strym_browser.html"), warn = FALSE), collapse = "\n"), '<header class="top"><img alt="StRym"')
})

test_that("the logo is embedded when the app folder has it and left out quietly when it does not", {
  expect_match(strym_logo_img(here), "^<img alt=\"StRym\".*src=\"data:image/png;base64,iVBOR"); expect_match(strym_logo_img(here, 30), "height:30px")
  expect_equal(strym_logo_img(NULL), ""); expect_equal(strym_logo_img(""), ""); expect_equal(strym_logo_img(tempdir()), "")
  d <- tempfile("noapp"); dir.create(d); sm <- data.frame(step = "pca", label = "PCA", status = "ok", seconds = 1, folder = "04_pca", files = "", message = "", stringsAsFactors = FALSE)
  f <- strym_write_index(d, design(project_name = "No logo"), sm)                                            # no app folder given: the page still works
  h <- paste(readLines(f, warn = FALSE), collapse = "\n"); expect_match(h, "No logo"); expect_false(grepl("data:image/png", h, fixed = TRUE))
  expect_true(all(file.exists(file.path(here, "www", c("strym-logo.png", "strym-logo-white.png", "favicon.png", "favicon-192.png")))))
})

test_that("a single step can be run on its own, leaving other folders untouched", {
  root <- tempfile("run"); dir.create(root)
  d <- demo_design(root)
  res <- quiet(strym_run(d, steps = "pca", app_dir = here))
  expect_equal(res$summary$step, "pca")
  pd <- file.path(root, "project")
  expect_true(file.exists(file.path(pd, "04_pca", "pca_plot.png")))
  expect_false(dir.exists(file.path(pd, "02_differential_expression")))
})

test_that("re-running a step replaces its earlier results, and leaves the other steps' folders alone", {
  root <- tempfile("run"); dir.create(root)
  d <- demo_design(root, de_methods = "deseq2", plots = "pca")
  quiet(strym_run(d, app_dir = here)); pd <- file.path(root, "project")
  writeLines("old", file.path(pd, "04_pca", "pca_placeholder.svg")); writeLines("keep", file.path(pd, "02_differential_expression", "deseq2", "mine.txt"))
  quiet(strym_run(d, steps = "pca", app_dir = here))
  expect_false(file.exists(file.path(pd, "04_pca", "pca_placeholder.svg"))); expect_true(file.exists(file.path(pd, "04_pca", "pca_plot.png")))
  expect_true(file.exists(file.path(pd, "02_differential_expression", "deseq2", "mine.txt")))
})

test_that("running one step after a full run updates only its row of the summary", {
  root <- tempfile("run"); dir.create(root)
  d <- demo_design(root, de_methods = "deseq2", plots = "pca")
  quiet(strym_run(d, app_dir = here))                                         # deseq2, pca
  res <- quiet(strym_run(d, steps = "pca", app_dir = here))
  expect_equal(res$summary$step, "pca")                                      # this run: one step ...
  s <- read.csv(file.path(root, "project", "run_summary.csv"))
  expect_equal(s$step, c("deseq2", "pca"))                                   # ... but the summary still lists both, in order
  expect_match(paste(readLines(file.path(root, "project", "index.html")), collapse = "\n"), "DESeq2 differential expression")
})

test_that("a step that cannot run is skipped with the reason; the rest still run", {
  root <- tempfile("run"); dir.create(root)
  d <- demo_design(root)
  res <- quiet(strym_run(d, steps = c("trajectory", "pca"), app_dir = here))     # trajectory lacks its prerequisite
  s <- setNames(res$summary$status, res$summary$step)
  expect_equal(s[["trajectory"]], "skipped"); expect_equal(s[["pca"]], "ok")
  expect_match(res$summary$message[res$summary$step == "trajectory"], "Needs one of")
  expect_false(dir.exists(file.path(root, "project", "05_trajectories")))
})

test_that("a failing step is recorded, dependants are skipped, other steps still run", {
  app <- tempfile("app"); dir.create(app); file.copy(file.path(here, "analysis"), app, recursive = TRUE)
  writeLines('strym_step_shannon_entropy <- function(ctx) stop("boom")', file.path(app, "analysis", "shannon_entropy.R"))
  root <- tempfile("run"); dir.create(root)
  d <- synth_design(synth_tux(root, orig = FALSE))
  res <- quiet(strym_run(d, steps = c("shannon_entropy", "trajectory", "pca"), app_dir = app))
  s <- setNames(res$summary$status, res$summary$step)
  expect_equal(s[["shannon_entropy"]], "failed"); expect_equal(s[["trajectory"]], "skipped"); expect_equal(s[["pca"]], "ok")
  expect_match(res$summary$message[res$summary$step == "shannon_entropy"], "boom")
})

test_that("a real script dropped into analysis/ replaces the placeholder and receives the context", {
  app <- tempfile("app"); dir.create(app); file.copy(file.path(here, "analysis"), app, recursive = TRUE)
  writeLines(c('strym_step_pca <- function(ctx) {',
               '  write.csv(data.frame(n_cond = ctx$design$n_conditions, data = basename(ctx$data_path)), file.path(ctx$out_dir, "pca_scores.csv"), row.names = FALSE)',
               '  list(status = "ok", files = "pca_scores.csv", message = "real script")', '}'), file.path(app, "analysis", "pca.R"))
  root <- tempfile("run"); dir.create(root)
  res <- quiet(strym_run(demo_design(root), steps = "pca", app_dir = app))
  expect_equal(res$summary$status, "ok")
  expect_equal(read.csv(file.path(root, "project", "04_pca", "pca_scores.csv"))$data, "demo_counts.csv")
})

test_that("run refuses an invalid design or unknown steps", {
  expect_error(strym_run(design(), steps = "pca", app_dir = here), "name")
  root <- tempfile("run"); dir.create(root)
  expect_error(strym_run(demo_design(root), steps = "nope", app_dir = here), "Unknown step")
  expect_error(strym_run(demo_design(root), steps = character(), app_dir = here), "No analysis steps|name")
})

test_that("command-line runner: DONE on success, ERROR (exit 1) when a step fails", {
  root <- tempfile("cli"); dir.create(root); d <- demo_design(root, plots = "pca")
  cfg <- file.path(root, "cfg.json")
  jsonlite::write_json(list(design = d, steps = "pca"), cfg, auto_unbox = TRUE)
  res <- suppressWarnings(system2("Rscript", c(shQuote(file.path(here, "run_analysis_cli.R")), shQuote(cfg)), stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("^DONE: .*index.html", res))); expect_null(attr(res, "status")); expect_true(any(grepl("\\[PROGRESS 100\\]", res)))
  bad <- d; bad$project_name <- ""
  jsonlite::write_json(list(design = bad, steps = "pca"), cfg, auto_unbox = TRUE)
  res <- suppressWarnings(system2("Rscript", c(shQuote(file.path(here, "run_analysis_cli.R")), shQuote(cfg)), stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("^ERROR: ", res))); expect_equal(attr(res, "status"), 1)
})

test_that("demo data has the promised shape", {
  d <- tempfile("demo"); r <- strym_make_demo_data(d)
  m <- read.csv(r$counts, check.names = FALSE); ss <- read.csv(r$sample_sheet)
  expect_equal(dim(m), c(300, 13)); expect_equal(nrow(ss), 12); expect_equal(names(m)[-1], ss$sample)
  expect_true(all(m[, -1] >= 0))
})

# ---- the exportable interactive browser (.html) ----
test_that("the interactive browser is ONE self-contained html with the plots and tables inside", {
  source(file.path(here, "analysis", "web_browser.R"), local = TRUE)
  pd <- tempfile("proj"); dir.create(file.path(pd, "04_pca"), recursive = TRUE); dir.create(file.path(pd, "02_de"))
  writeLines('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><rect width="10" height="10"/></svg>', file.path(pd, "04_pca", "pca.svg"))
  write.csv(data.frame(gene = c("A", "B</script>", 'say "hi"'), log2FC = c(1.5, -2, 0)), file.path(pd, "02_de", "genes.csv"), row.names = FALSE)
  writeLines("secret", file.path(pd, "run_summary.csv"))                                         # not part of the tables
  out <- file.path(pd, "07_interactive_browser", "strym_browser.html"); dir.create(dirname(out))
  res <- strym_browser_html(pd, design(project_name = "Demo <b>", n_conditions = 2, n_timepoints = 3), out)
  expect_equal(c(res$plots, res$tables), c(1, 1))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_false(grepl("src=\"http|href=\"http|<script src|<link ", html))                       # nothing is loaded from the internet
  expect_match(html, "Jordan Reed"); expect_match(html, "Sozzani lab")
  json <- sub('(?s).*<script id="data" type="application/json">(.*?)</script>.*', "\\1", html, perl = TRUE)
  expect_false(grepl("</script", json, fixed = TRUE))                                           # cell contents cannot close the tag
  dat <- jsonlite::fromJSON(gsub("<\\/", "</", json, fixed = TRUE), simplifyVector = FALSE)
  expect_equal(dat$title, "Demo <b>"); expect_match(dat$plots[[1]]$src, "^data:image/svg\\+xml;base64,")
  expect_equal(dat$tables[[1]]$name, "02_de/genes.csv"); expect_equal(unlist(dat$tables[[1]]$cols), c("gene", "log2FC"))
  expect_equal(vapply(dat$tables[[1]]$rows, function(r) r[[1]], ""), c("A", "B</script>", 'say "hi"'))
  expect_false(any(grepl("run_summary|00_design", c(dat$plots[[1]]$name, dat$tables[[1]]$name))))
})

test_that("the browser step reports its own result and an empty project still gives a valid page", {
  root <- tempfile("run"); dir.create(root)
  d <- demo_design(root, de_methods = "deseq2", interactive_browser = TRUE)
  res <- quiet(strym_run(d, steps = c("deseq2", "web_browser"), app_dir = here))
  r <- res$summary[res$summary$step == "web_browser", ]
  expect_equal(r$status, "ok"); expect_match(r$files, "strym_browser.html"); expect_match(r$message, "\\([1-9][0-9]* plots, [1-9][0-9]* tables\\)")   # DESeq2's volcano plot and tables are embedded
  f <- file.path(root, "project", "07_interactive_browser", "strym_browser.html"); expect_true(file.exists(f)); expect_gt(file.size(f), 3000)
})
