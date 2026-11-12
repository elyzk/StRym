#!/usr/bin/env Rscript
# Command-line runner for the automated analysis. Used by the GUI as a background job, usable directly:
#   Rscript run_analysis_cli.R config.json
# config.json: { "design": {...}, "steps": ["deseq2", "pca", ...] }   (steps optional: default = the plan's defaults)
.unesc <- function(x) gsub("~+~", " ", x, fixed = TRUE)       # Rscript writes a space inside a path or argument as "~+~" (a folder such as "StRym 3")
args <- .unesc(commandArgs(trailingOnly = TRUE))
if (length(args) != 1) stop("usage: Rscript run_analysis_cli.R config.json")
here <- dirname(normalizePath(.unesc(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
source(file.path(here, "io_utils.R")); source(file.path(here, "tuxedo.R")); source(file.path(here, "pipeline.R"))
cfg <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)
status <- tryCatch({
  res <- withCallingHandlers(                                  # R's "package X was built under R version Y" notes are not the researcher's business
    strym_run(cfg$design, cfg$steps, app_dir = here, progress = function(pct, label) emit_progress(pct, label)),
    warning = function(w) if (grepl("was built under R version", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
  failed <- res$summary$step[res$summary$status == "failed"]
  if (length(failed)) { message("ERROR: step(s) failed: ", paste(failed, collapse = ", "), " (see the summary in ", res$project_dir, ")"); 1L }
  else { message("DONE: ", res$index); 0L }
}, error = function(e) { message("ERROR: ", conditionMessage(e)); 1L })
quit(status = status)
