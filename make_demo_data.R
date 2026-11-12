# Regenerates the demo experiment used by the "Load the demo experiment" button:
#   Rscript make_demo_data.R
# 300 genes, 2 conditions (control / heat) x 3 time points (0h, 1h, 4h) x 2 replicates = 12 samples.
here <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)))   # "~+~" = a space in the folder name
source(file.path(here, "io_utils.R")); source(file.path(here, "pipeline.R"))
print(strym_make_demo_data(file.path(here, "example_data", "strym_demo")))
