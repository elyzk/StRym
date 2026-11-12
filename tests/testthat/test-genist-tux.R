# GENIST inputs built from Tux (Cuffdiff) results: the right genes, the time points in the right order, clear problems.

# A small made-up Tux result inside a StRym project folder. Cuffdiff lists the groups in the order of the sample sheet
# (lateSE first), which is NOT the time order (earlyInit, midProto, lateSE): the builder must put them in time order.
tux_toy <- function(design = TRUE, tps = c("earlyInit", "midProto", "lateSE"), flat = FALSE) {
  root <- tempfile("toy"); dd <- if (flat) file.path(root, "diff_out") else file.path(root, "01_quantification", "tux", "diff_out")
  dir.create(dd, recursive = TRUE)
  if (design) { dir.create(file.path(root, "00_design"))
    writeLines(jsonlite::toJSON(list(condition_names = "-Fe", timepoint_names = tps, n_timepoints = length(tps)), auto_unbox = TRUE),
               file.path(root, "00_design", "design.json")) }
  #            gene         FPKM late, mid, early      q (late-mid, late-early, mid-early)          log2FC (same order)           status
  g <- list(AT1G00010 = list(c(100, 20, 5),   c(.001, .001, .001), c(-2.3, -4.3, -2.0), "OK"),     # changes a lot, significant: chosen
            AT1G00020 = list(c(10, 12, 11),   c(.9, .9, .9),       c(.2, .1, -.1),     "OK"),      # no change
            AT1G00030 = list(c(50, 5, 5),     c(.2, .2, .9),       c(-3, -3, 0),       "OK"),      # large but not significant
            AT1G00040 = list(c(10, 12, 14),   c(.01, .01, .01),    c(.3, .5, .2),      "OK"),      # significant but small
            AT1G00050 = list(c(.5, .1, .2),   c(.001, .001, .001), c(-2, -2, 1),       "OK"),      # changes but hardly expressed
            AT1G00060 = list(c(80, 8, 8),     c(.001, .001, .001), c(-3, -3, 0),       "NOTEST"),  # not tested
            AT1G00070 = list(c(0, 8, 9),      c(.001, .001, .9),   c(Inf, Inf, 0),     "OK"),      # infinite fold change: chosen
            AT1G00080 = list(c(30, 28, 3),    c(.9, .001, .001),   c(0, -3.2, -3.1),   "OK"))      # only some pairs differ: chosen
  pairs <- list(c("-Fe_lateSE", "-Fe_midProto"), c("-Fe_lateSE", "-Fe_earlyInit"), c("-Fe_midProto", "-Fe_earlyInit"))
  de <- do.call(rbind, lapply(names(g), function(n) do.call(rbind, lapply(1:3, function(i)
    data.frame(test_id = paste0("X_", n), gene_id = paste0("X_", n), gene = n, locus = "1:1-2", sample_1 = pairs[[i]][1], sample_2 = pairs[[i]][2],
               status = g[[n]][[4]], value_1 = 1, value_2 = 1, `log2(fold_change)` = g[[n]][[3]][i], test_stat = 0, p_value = 0.01,
               q_value = g[[n]][[2]][i], significant = "yes", check.names = FALSE, stringsAsFactors = FALSE)))))
  de <- rbind(de, data.frame(test_id = "X_nn", gene_id = "X_nn", gene = "-", locus = "1:1-2", sample_1 = "-Fe_lateSE", sample_2 = "-Fe_midProto", status = "OK",
                             value_1 = 1, value_2 = 1, `log2(fold_change)` = -5, test_stat = 0, p_value = 0, q_value = 0, significant = "yes", check.names = FALSE))
  de$`log2(fold_change)` <- ifelse(is.infinite(de$`log2(fold_change)`), "inf", de$`log2(fold_change)`)
  utils::write.table(de, file.path(dd, "gene_exp.diff"), sep = "\t", quote = FALSE, row.names = FALSE)
  tr <- data.frame(tracking_id = paste0("X_", names(g)), class_code = "-", nearest_ref_id = "-", gene_id = paste0("X_", names(g)), gene_short_name = names(g),
                   tss_id = "-", locus = "1:1-2", length = "-", coverage = "-", check.names = FALSE, stringsAsFactors = FALSE)
  for (grp in c("lateSE", "midProto", "earlyInit")) { j <- match(grp, c("lateSE", "midProto", "earlyInit"))
    tr[[paste0("-Fe_", grp, "_FPKM")]] <- vapply(g, function(x) x[[1]][j], 0); tr[[paste0("-Fe_", grp, "_conf_lo")]] <- 0
    tr[[paste0("-Fe_", grp, "_conf_hi")]] <- 0; tr[[paste0("-Fe_", grp, "_status")]] <- "OK" }
  extra <- tr[1, ]; extra$tracking_id <- "X_dup"; extra$gene_id <- "X_dup"; extra[, grep("_FPKM$", names(extra))] <- 1       # a second locus with the same gene name, less expressed
  nn <- tr[1, ]; nn$tracking_id <- "X_nn"; nn$gene_id <- "X_nn"; nn$gene_short_name <- "-"; nn[, grep("_FPKM$", names(nn))] <- 999     # a locus without a gene name
  utils::write.table(rbind(tr, extra, nn), file.path(dd, "genes.fpkm_tracking"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(root = root, diff = dd, out = file.path(root, "inputs"))
}

test_that("the Tux results are found from the project folder, the tux folder or diff_out; anything else gets a plain message", {
  fx <- tux_toy()
  for (p in c(fx$root, file.path(fx$root, "01_quantification"), file.path(fx$root, "01_quantification", "tux"), fx$diff)) {
    f <- tux_find_results(paste0(p, "/")); expect_equal(f$diff_dir, normalizePath(fx$diff, winslash = "/"))
  }
  expect_equal(tux_find_results(fx$root)$project_dir, normalizePath(fx$root, winslash = "/"))
  expect_equal(tux_find_results(fx$root)$design$timepoint_names, c("earlyInit", "midProto", "lateSE"))
  expect_error(tux_find_results(tempdir()), "No Tux results found")
  expect_error(tux_find_results(file.path(fx$root, "nowhere")), "Folder not found")
  expect_null(tux_find_results(tux_toy(design = FALSE)$diff)$design)
})

test_that("the time order comes from the design's time point names, not from Tux's own order", {
  fx <- tux_toy(); g <- tux_groups(tux_find_results(fx$root))
  expect_false(g$guessed); expect_equal(g$default, "-Fe")
  expect_equal(g$series[["-Fe"]], c(earlyInit = "-Fe_earlyInit", midProto = "-Fe_midProto", lateSE = "-Fe_lateSE"))
  expect_equal(g$available[1], "-Fe_lateSE")                                       # Tux's own order starts with lateSE
  g2 <- tux_groups(tux_find_results(tux_toy(design = FALSE)$diff))                  # no design record: Tux's order, flagged as a guess
  expect_true(g2$guessed); expect_equal(unname(g2$series[[1]]), c("-Fe_lateSE", "-Fe_midProto", "-Fe_earlyInit"))
})

test_that("a control condition is not chosen as the time course by default", {
  fx <- tux_toy(); f <- file.path(fx$diff, "genes.fpkm_tracking")
  tr <- read.delim(f, check.names = FALSE, quote = "", stringsAsFactors = FALSE)       # add a control condition's groups, listed FIRST
  ctl <- tr[0]; for (grp in c("lateSE", "midProto", "earlyInit")) ctl[[paste0("control_", grp, "_FPKM")]] <- rep(1, nrow(tr))
  utils::write.table(cbind(tr[1:9], ctl, tr[-(1:9)]), f, sep = "\t", quote = FALSE, row.names = FALSE)
  g <- tux_groups(tux_find_results(fx$root))
  expect_equal(names(g$series), c("control", "-Fe")); expect_equal(g$default, "-Fe")
  expect_equal(names(g$series[["control"]]), c("earlyInit", "midProto", "lateSE"))
})

test_that("the genes file and the time course file are built correctly", {
  fx <- tux_toy(); r <- genist_inputs_from_tux(fx$root, fx$out)
  tc <- read.csv(r$time_course_file, check.names = FALSE); gn <- read.csv(r$genes_file)
  expect_equal(names(tc), c("gene", "earlyInit", "midProto", "lateSE"))                # time order, whatever order Tux used
  expect_equal(tc$gene, c(paste0("AT1G000", c("10", "20", "30", "40", "50", "60", "70", "80"))))   # one row per named gene; "-" loci and the duplicate are gone
  expect_equal(unlist(tc[tc$gene == "AT1G00010", -1]), c(earlyInit = 5, midProto = 20, lateSE = 100))   # the most expressed locus of a repeated name
  expect_equal(gn$gene, c("AT1G00010", "AT1G00070", "AT1G00080"))                      # significant, >= 1 log2 change, expressed, tested, named
  expect_equal(r$time_points, c("earlyInit", "midProto", "lateSE")); expect_false(r$guessed); expect_false(r$capped); expect_equal(r$n_genes, 3)
  expect_match(r$summary, "3 genes chosen"); expect_match(r$summary, "earlyInit → midProto → lateSE")
  # looser / stricter filters
  expect_equal(read.csv(genist_inputs_from_tux(fx$root, fx$out, min_log2fc = 0.25)$genes_file)$gene, c("AT1G00010", "AT1G00040", "AT1G00070", "AT1G00080"))
  expect_equal(read.csv(genist_inputs_from_tux(fx$root, fx$out, min_fpkm = 0)$genes_file)$gene, c("AT1G00010", "AT1G00050", "AT1G00070", "AT1G00080"))
  expect_equal(read.csv(genist_inputs_from_tux(fx$root, fx$out, q_cutoff = 0.5)$genes_file)$gene, c("AT1G00010", "AT1G00030", "AT1G00070", "AT1G00080"))
  # the limit keeps the strongest changes (AT1G00010: 4.3, AT1G00070: infinite, AT1G00080: 3.2)
  r2 <- genist_inputs_from_tux(fx$root, fx$out, max_genes = 2); expect_true(r2$capped); expect_equal(r2$n_pass, 3)
  expect_equal(read.csv(r2$genes_file)$gene, c("AT1G00010", "AT1G00070")); expect_match(r2$summary, "kept the 2 strongest of 3")
  expect_false(genist_inputs_from_tux(fx$root, fx$out, max_genes = 0)$capped)
  # the files can be written to a folder that does not exist yet, and nothing in the Tux folder changed
  before <- tools::md5sum(file.path(fx$diff, c("genes.fpkm_tracking", "gene_exp.diff")))
  genist_inputs_from_tux(fx$diff, file.path(fx$root, "a", "b", "c")); expect_true(file.exists(file.path(fx$root, "a", "b", "c", "genist_genes_from_tux.csv")))
  expect_equal(tools::md5sum(file.path(fx$diff, c("genes.fpkm_tracking", "gene_exp.diff"))), before)
})

test_that("the time order can be edited, and a wrong or too short one is refused", {
  fx <- tux_toy()
  r <- genist_inputs_from_tux(fx$root, fx$out, order = c("lateSE", " midProto", "earlyInit "))
  expect_equal(names(read.csv(r$time_course_file, check.names = FALSE)), c("gene", "lateSE", "midProto", "earlyInit"))
  expect_error(genist_inputs_from_tux(fx$root, fx$out, order = c("earlyInit", "midProto", "late")), "not in the Tux results: late")
  expect_error(genist_inputs_from_tux(fx$root, fx$out, order = c("earlyInit", "midProto")), "at least 3 time points")
  expect_error(genist_inputs_from_tux(fx$root, fx$out, order = c("earlyInit", "midProto", "midProto")), "listed twice")
  expect_error(genist_inputs_from_tux(fx$root, fx$out, q_cutoff = 0), "between 0 and 1")
  expect_error(genist_inputs_from_tux(fx$root, fx$out, q_cutoff = 1e-9), "Only 0 gene")
  expect_error(genist_inputs_from_tux(tempdir(), fx$out), "No Tux results found")
  rg <- genist_inputs_from_tux(tux_toy(design = FALSE)$diff, fx$out); expect_true(rg$guessed)                   # without a design record: Tux's order, flagged
  expect_equal(rg$time_points, c("-Fe_lateSE", "-Fe_midProto", "-Fe_earlyInit"))
  two <- tux_toy(tps = c("earlyInit", "midProto")); expect_error(genist_inputs_from_tux(two$root, two$out), "at least 3 time points")
})

test_that("GENIST runs on the files built from Tux results", {
  fx <- tux_toy(); r <- genist_inputs_from_tux(fx$root, fx$out, min_log2fc = 0.25, min_fpkm = 0)
  out <- file.path(fx$root, "net")
  quiet(genist_run(r$genes_file, r$time_course_file, output_folder = out, fc_thresh = 1.1, time_percent = 0.1, n_levels = 3, low_conn = 0.05))
  expect_true(file.exists(file.path(out, "cytoscape_table.txt"))); expect_true(file.exists(file.path(out, "clusters.xlsx")))
  expect_equal(names(read.delim(file.path(out, "cytoscape_table.txt"))), c("source", "target", "width", "sign"))
})
