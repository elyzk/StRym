#!/usr/bin/env Rscript
# Command-line runner for the TUX pipeline. Used by the GUI (app.R) as a background
# process, and usable directly:  Rscript run_cli.R config.json
# config.json keys = arguments of run_tuxedo() (see tuxedo.R).
.unesc <- function(x) gsub("~+~", " ", x, fixed = TRUE)       # Rscript writes a space inside a path or argument as "~+~" (a folder such as "StRym 3")
args <- .unesc(commandArgs(trailingOnly = TRUE))
if (length(args) != 1) stop("usage: Rscript run_cli.R config.json")
here <- dirname(normalizePath(.unesc(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
source(file.path(here, "tuxedo.R"))
cfg <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)
if (is.null(cfg$adapter_fasta) || !nzchar(cfg$adapter_fasta))
  cfg$adapter_fasta <- file.path(here, "IlluminaAdaptorSeq.fasta")
if (!is.null(cfg$tool_dir) && !nzchar(cfg$tool_dir)) cfg$tool_dir <- NULL
res <- tryCatch({
  do.call(run_tuxedo, cfg)
  0L
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  1L
})
quit(status = res)
