# Run all tests:  Rscript tests/run_tests.R
library(testthat)
args <- commandArgs(FALSE); f <- gsub("~+~", " ", sub("^--file=", "", args[grep("^--file=", args)]), fixed = TRUE)   # "~+~" = a space in the folder name
setwd(file.path(dirname(normalizePath(f)), "testthat"))
res <- test_dir(".", reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat("\nPASSED:", sum(df$passed), " FAILED:", sum(df$failed), " ERRORS:", sum(df$error), " SKIPPED:", sum(df$skipped), "\n")
quit(status = as.integer(sum(df$failed) + sum(df$error) > 0))
