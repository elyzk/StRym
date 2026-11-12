# The program folder may have a space in its name (macOS names a second download "StRym 3"). Rscript then writes the space as "~+~" in the script's
# own path and in its arguments; the background runners must put it back, or they cannot find their own files.
test_that("Rscript really does hide a space as ~+~ (the reason the runners decode their own path)", {
  skip_if(!nzchar(Sys.which("Rscript")))
  d <- file.path(tempfile("sp"), "StRym 3"); dir.create(d, recursive = TRUE); f <- file.path(d, "show.R")
  writeLines('cat(grep("^--file=", commandArgs(FALSE), value = TRUE), commandArgs(TRUE), sep = "\\n")', f)
  out <- system2(Sys.which("Rscript"), c(shQuote(f), shQuote(file.path(d, "my config.json"))), stdout = TRUE, stderr = TRUE)
  expect_true(any(grepl("StRym~+~3", out, fixed = TRUE)) || any(grepl("StRym 3", out, fixed = TRUE)))     # either spelling: what matters is that the runners cope with both
})

test_that("the background runners start from a folder whose name has a space, and read a config file whose path has one", {
  skip_if(!nzchar(Sys.which("Rscript")))
  root <- file.path(tempfile("sp"), "StRym 3"); dir.create(file.path(root, "analysis"), recursive = TRUE)
  file.copy(list.files(here, pattern = "\\.R$", full.names = TRUE), root); file.copy(list.files(file.path(here, "analysis"), pattern = "\\.R$", full.names = TRUE), file.path(root, "analysis"))
  cfg <- file.path(root, "my config.json"); writeLines('{"genes_file":"nope.csv","time_course_file":"nope2.csv","output_folder":"out"}', cfg)
  run <- function(script) system2(Sys.which("Rscript"), c(shQuote(file.path(root, script)), shQuote(cfg)), stdout = TRUE, stderr = TRUE)
  out <- run("run_genist_cli.R")
  expect_true(any(grepl("^ERROR:.*cannot be found", out)), info = paste(out, collapse = " | "))             # it got as far as GENIST itself, which reports the missing file
  expect_false(any(grepl("~+~|cannot open (the connection|file)|No such file", out)))
  for (s in c("run_rtpstar_cli.R", "run_cli.R", "run_analysis_cli.R")) {                                       # the other runners find their own files too
    o <- run(s); expect_false(any(grepl("~+~|cannot open (the connection|file)", o)), info = paste(s, paste(o, collapse = " | ")))
    expect_true(any(grepl("^ERROR:", o)), info = paste(s, paste(o, collapse = " | ")))                           # ... and fail only because this config is not a real one
  }
})
