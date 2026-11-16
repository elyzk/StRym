`%||%` <- function(a, b) if (is.null(a)) b else a
here <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(here, "genist.R"))) here <- normalizePath("../..")
source(file.path(here, "io_utils.R")); source(file.path(here, "genist.R")); source(file.path(here, "rtpstar.R"))
source(file.path(here, "tuxedo.R")); source(file.path(here, "pipeline.R")); strym_load_registry(here)
source(file.path(here, "drives.R"))
suppressPackageStartupMessages(library(shiny)); source(file.path(here, "filebrowser.R")); source(file.path(here, "gui_helpers.R"))   # Browse, whole_threads(), path fields
quiet <- function(expr) suppressMessages(expr)
with_platform <- function(p, code) { old <- options(strym.platform = p); on.exit(options(old)); force(code) }
