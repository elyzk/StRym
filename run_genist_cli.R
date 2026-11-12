#!/usr/bin/env Rscript
# Command-line runner for GENIST. Used by the GUI as a background job, and usable directly:
#   Rscript run_genist_cli.R config.json     (keys = arguments of genist_run(), see genist.R)
.unesc <- function(x) gsub("~+~", " ", x, fixed = TRUE)       # Rscript writes a space inside a path or argument as "~+~" (a folder such as "StRym 3")
args <- .unesc(commandArgs(trailingOnly = TRUE))
if (length(args) != 1) stop("usage: Rscript run_genist_cli.R config.json")
here <- dirname(normalizePath(.unesc(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
source(file.path(here, "io_utils.R")); source(file.path(here, "genist.R"))
cfg <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)
status <- tryCatch({
  res <- do.call(genist_run, cfg)
  message("DONE: ", res$output)
  0L
}, error = function(e) { message("ERROR: ", conditionMessage(e)); 1L })
quit(status = status)
