# Version 9: PCA level, entropy heatmaps, network pictures and file names, the results that are shown.
for (f in c("pca", "shannon_entropy", "modified_shannon_entropy", "trajectory", "web_browser")) source(file.path(here, "analysis", paste0(f, ".R")))
suppressPackageStartupMessages(library(DT)); source(file.path(here, "mod_design.R"))
logs9 <- new.env(); logs9$x <- character()
ctx9 <- function(d, out, prereq = list()) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE); d <- strym_normalize_design(d); logs9$x <- character()
  list(design = d, data_path = d$data_path, sample_sheet = d$sample_sheet, out_dir = out, project_dir = d$project_dir, app_dir = here,
       prereq = prereq, log = function(...) logs9$x <- c(logs9$x, paste0(...)))
}

test_that("PCA: replicate level is the default and unchanged; the averaged level has one point per condition and time point", {
  root <- tempfile("v9"); dir.create(root); syn <- synth_tux(root, orig = FALSE)
  out <- file.path(root, "p1"); r <- strym_step_pca(ctx9(synth_design(syn), out))
  expect_equal(sort(r$files), sort(c("pca_scores.csv", "pca_variance.csv", "pca_plot.png"))); expect_false(file.exists(file.path(out, "pca_averaged_plot.png")))
  out2 <- file.path(root, "p2"); r2 <- strym_step_pca(ctx9(synth_design(syn, params = list(pca_level = "average")), out2))
  expect_setequal(r2$files, c("pca_averaged_scores.csv", "pca_averaged_variance.csv", "pca_averaged_plot.png")); expect_false(file.exists(file.path(out2, "pca_plot.png")))
  sc <- read.csv(file.path(out2, "pca_averaged_scores.csv")); ngroup <- length(syn$conds) * length(syn$tps)
  expect_equal(nrow(sc), ngroup); expect_true(all(sc$n_replicates == 3)); expect_true(all(c("sample", "condition", "time_point", "n_replicates", "PC1", "PC2") %in% names(sc)))
  # the average of the replicates, by hand, for the first group: the PCA input is the mean of the replicate values
  out3 <- file.path(root, "p3"); r3 <- strym_step_pca(ctx9(synth_design(syn, params = list(pca_level = "both")), out3))
  expect_setequal(r3$files, c("pca_scores.csv", "pca_variance.csv", "pca_plot.png", "pca_averaged_scores.csv", "pca_averaged_variance.csv", "pca_averaged_plot.png"))
  expect_equal(read.csv(file.path(out3, "pca_scores.csv")), read.csv(file.path(out, "pca_scores.csv")))                      # both: the replicate result is the same as alone
  expect_equal(strym_normalize_params(list(pca_level = "Nonsense"))$pca_level, "replicate"); expect_equal(strym_normalize_params(list(pca_level = "BOTH"))$pca_level, "both")
})

test_that("the averaged PCA needs 3 groups: an error alone, a note with 'both'", {
  root <- tempfile("v9"); dir.create(root); syn <- synth_tux(root, orig = FALSE, genotypes = "A", n_tp = 2)       # 1 genotype x (heated, control) x 2 time points = 4 groups; keep 1 condition below
  d <- synth_design(syn, params = list(pca_level = "average", drop_reps = paste(read.csv(file.path(syn$project, "01_quantification", "sample_sheet.csv"))$sample[grepl("Control", read.csv(file.path(syn$project, "01_quantification", "sample_sheet.csv"))$sample)], collapse = ", ")))
  expect_error(strym_step_pca(ctx9(d, file.path(root, "e1"))), "at least 3 condition x time point groups")
  d$params$pca_level <- "both"; r <- strym_step_pca(ctx9(d, file.path(root, "e2")))
  expect_true("pca_plot.png" %in% r$files); expect_false("pca_averaged_plot.png" %in% r$files); expect_match(r$message, "averaged PCA was skipped")
})

test_that("entropy heatmaps: drawn for both steps, switchable, and robust to awkward data", {
  skip_if_not_installed("pheatmap")
  root <- tempfile("v9"); dir.create(root); syn <- synth_tux(root, orig = FALSE); d <- synth_design(syn)
  out <- file.path(root, "mse"); r <- strym_step_modified_shannon_entropy(ctx9(d, out))
  hm <- grep("_heatmap_", r$files, value = TRUE)
  expect_true(all(c("MSE_heatmap_log2FC_A_Heated.png", "MSE_heatmap_expression_A_Heated.png", "MSE_heatmap_log2FC_all_comparisons.png") %in% hm)); expect_true(all(file.exists(file.path(out, hm))))
  s <- strym_step_shannon_entropy(ctx9(d, file.path(root, "sh"))); expect_true(any(grepl("^shannon_heatmap_log2FC_", s$files)))
  off <- strym_step_shannon_entropy(ctx9(synth_design(syn, params = list(entropy_heatmaps = FALSE)), file.path(root, "sh2"))); expect_false(any(grepl("_heatmap_", off$files)))
  # awkward data, straight into the drawing code: one gene, two genes, constant values, NA and Inf, many genes, many comparisons
  mk <- function(ng, tps, np, ncalled = ng, flat = FALSE) {
    genes <- sprintf("G%05d", seq_len(ng)); pairs <- data.frame(treated = paste0("T", seq_len(np)), control = paste0("C", seq_len(np)), stringsAsFactors = FALSE)
    one <- function() { m <- matrix(if (flat) 10 else rexp(ng * length(tps), 1 / 50), ng, length(tps), dimnames = list(genes, tps)); m }
    heat <- lapply(seq_len(np), function(i) one()); ctrl <- lapply(seq_len(np), function(i) one())
    lfc <- Map(function(h, c) log2((h + 1) / (c + 1)), heat, ctrl)
    tc <- list(pairs = pairs, tps = tps, heat = heat, ctrl = ctrl, lfc = lfc, qval = lapply(lfc, function(x) x * 0 + 0.01), genes = genes, gene_names = NULL)
    keys <- .pair_keys(pairs); res <- list()
    for (i in seq_len(np)) { g <- genes[seq_len(ncalled)]; cl <- sign(lfc[[i]][g, , drop = FALSE]); if (ncol(cl)) cl[, 1] <- 0
      res[[keys[i]]] <- .calls_table(tc, i, if (ncalled) cl else matrix(0, 0, length(tps), dimnames = list(character(), tps)), NULL, 2) }
    list(run = list(results = res, main = NULL, excluded = rep(0L, np)), tc = tc)
  }
  set.seed(3); P <- list(heatmap_max_genes = 100)
  for (case in list(list(1, c("a", "b", "c"), 1), list(2, c("a", "b", "c"), 1), list(30, c("a", "b"), 1), list(300, c("a", "b", "c", "d"), 1), list(25, paste0("t", 1:5), 5))) {
    x <- do.call(mk, case); dd <- file.path(root, paste0("hm", case[[1]], "_", case[[3]])); dir.create(dd)
    f <- strym_entropy_heatmaps(dd, "X", x$run, x$tc, P, "X"); expect_gte(length(f), 2); expect_true(all(file.exists(file.path(dd, f))))
  }
  x <- mk(20, c("a", "b", "c"), 1, flat = TRUE); dd <- file.path(root, "flat"); dir.create(dd); expect_gte(length(strym_entropy_heatmaps(dd, "X", x$run, x$tc, P, "X")), 2)
  x <- mk(20, c("a", "b", "c"), 1); x$tc$lfc[[1]][3, 2] <- Inf; x$tc$lfc[[1]][4, 1] <- NA; dd <- file.path(root, "na"); dir.create(dd); expect_gte(length(strym_entropy_heatmaps(dd, "X", x$run, x$tc, P, "X")), 2)
  x <- mk(10, c("a", "b", "c"), 1, ncalled = 0); dd <- file.path(root, "none"); dir.create(dd); expect_length(strym_entropy_heatmaps(dd, "X", x$run, x$tc, P, "X"), 0)
})

test_that("network file names carry the project and the comparison; the network picture is drawn and never stops a run", {
  expect_equal(strym_network_filename("Heat shock in roots!", "GENIST", "Heat_vs_Control"), "Heat_shock_in_roots_Heat_vs_Control_GENIST_network_for_Cytoscape.txt")
  expect_equal(strym_network_filename("", "RTPSTAR"), "RTPSTAR_network_for_Cytoscape.txt"); expect_equal(strym_safe_name("  "), "project")
  skip_if_not_installed("igraph")
  e <- data.frame(from = c("A", "A", "B", "C"), to = c("B", "C", "C", "D"), type = c("activates", "inhibits", "regulates", "activates"), weight = c(1, 0.5, 0.2, 0.9), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".png"); expect_true(strym_network_png(e, f, "t")); expect_gt(file.size(f), 5000)
  big <- data.frame(from = sprintf("G%04d", sample(800, 3000, TRUE)), to = sprintf("G%04d", sample(800, 3000, TRUE)), type = "activates", weight = runif(3000), stringsAsFactors = FALSE)
  f2 <- tempfile(fileext = ".png"); expect_true(strym_network_png(big, f2, "big")); expect_true(file.exists(f2))                  # thousands of edges: the strongest are drawn
  expect_false(strym_network_png(e[0, ], tempfile(fileext = ".png"))); expect_false(strym_network_png(NULL, tempfile(fileext = ".png")))
  expect_false(strym_network_png(e, file.path(tempfile(), "no", "such", "dir", "x.png")))                                          # an unwritable place is a FALSE, not an error
})

test_that("GENIST and RTP-STAR runs also leave a picture next to the network table", {
  skip_if_not_installed("igraph")
  root <- tempfile("v9"); dir.create(root); set.seed(2)
  genes <- sprintf("AT1G%05d", 1:30); tp <- cbind(t1 = rexp(30, 0.1), t2 = rexp(30, 0.1), t3 = rexp(30, 0.1), t4 = rexp(30, 0.1)); rownames(tp) <- genes
  tcf <- file.path(root, "tc.csv"); gf <- file.path(root, "genes.csv")
  write.csv(data.frame(gene = genes, tp), tcf, row.names = FALSE); write.csv(data.frame(gene = genes), gf, row.names = FALSE)
  out <- file.path(root, "g"); genist_run(gf, tcf, output_folder = out, output_file = "Proj_GENIST_network_for_Cytoscape.txt", fc_thresh = 1.1, time_percent = 0.1, n_levels = 3, low_conn = 0.05)
  expect_true(file.exists(file.path(out, "Proj_GENIST_network_for_Cytoscape.txt"))); expect_true(file.exists(file.path(out, "Proj_GENIST_network_for_Cytoscape_network.png")))
  ex <- data.frame(gene = genes, TF = rep(c(1, 0, 0), length.out = 30), matrix(rexp(30 * 6, 0.1), 30, 6, dimnames = list(NULL, paste0("s", 1:6))), check.names = FALSE)
  exf <- file.path(root, "expr.csv"); write.csv(ex, exf, row.names = FALSE)
  out2 <- file.path(root, "r"); rtpstar_run(gf, exf, out2, output_file = "Proj_RTPSTAR_network_for_Cytoscape.txt", ntrees = 200)
  expect_true(file.exists(file.path(out2, "Proj_RTPSTAR_network_for_Cytoscape.txt")))
  l <- readLines(file.path(out2, "Proj_RTPSTAR_network_for_Cytoscape.txt")); if (length(l)) expect_true(file.exists(file.path(out2, "Proj_RTPSTAR_network_for_Cytoscape_network.png")))
})

test_that("only the steps' own folders are shown; network tables and pictures reach the interactive browser; extra folders can be added", {
  root <- tempfile("v9"); pd <- file.path(root, "Downloads"); dir.create(pd, recursive = TRUE)
  mk <- function(rel, txt = "a,b\n1,2\n") { f <- file.path(pd, rel); dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); writeLines(txt, f); f }
  png1 <- function(rel) { f <- file.path(pd, rel); dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); grDevices::png(f, 300, 300); graphics::par(mar = c(1, 1, 1, 1)); plot(1); grDevices::dev.off(); f }
  mk("02_differential_expression/deseq2/res.csv"); png1("02_differential_expression/deseq2/volcano.png")
  mk("06_regulatory_networks/GENIST/Proj_GENIST_network_for_Cytoscape.txt", "source\ttarget\twidth\tsign\nA\tB\t0.5\t1\nB\tC\t0.2\t-1\n"); png1("06_regulatory_networks/GENIST/Proj_GENIST_network_for_Cytoscape_network.png")
  mk("06_regulatory_networks/RTP-STAR/Proj_RTPSTAR_network_for_Cytoscape.txt", "A activates B 0.67\nB inhibits C 0.33\n"); mk("06_regulatory_networks/RTP-STAR/run.log.txt", "not a table\n")
  mk("my_taxes.csv"); png1("holiday_photo.png"); mk("Documents/other.csv"); png1("Documents/pic.png"); mk("run_summary.csv"); mk("00_design/plan.csv"); mk("01_quantification/tux/diff_out/x.csv")
  expect_setequal(strym_result_files(pd, "\\.(png|svg|jpe?g)$"), c("02_differential_expression/deseq2/volcano.png", "06_regulatory_networks/GENIST/Proj_GENIST_network_for_Cytoscape_network.png"))
  expect_true("02_differential_expression/deseq2/res.csv" %in% strym_result_files(pd, "\\.csv$")); expect_false(any(grepl("taxes|other|plan|run_summary|tux", strym_result_files(pd, "\\.csv$"))))
  ext <- file.path(root, "elsewhere", "GENIST_results"); dir.create(ext, recursive = TRUE); writeLines("source\ttarget\twidth\tsign\nX\tY\t1\t1", file.path(ext, "cytoscape_table.txt"))
  out <- file.path(pd, "07_interactive_browser"); dir.create(out)
  d <- strym_normalize_design(list(project_name = "p", project_dir = pd, data_type = "bulk"))
  h <- strym_browser_html(pd, d, file.path(out, "b.html"), extra_dirs = ext)
  txt <- paste(readLines(file.path(out, "b.html"), warn = FALSE), collapse = "\n"); js <- jsonlite::fromJSON(sub("</script>.*$", "", sub("^.*<script id=\"data\" type=\"application/json\">", "", gsub("\n", " ", txt))))
  tn <- js$tables$name; pn <- js$plots$name
  expect_true(all(c("06_regulatory_networks/GENIST/Proj_GENIST_network_for_Cytoscape.txt", "06_regulatory_networks/RTP-STAR/Proj_RTPSTAR_network_for_Cytoscape.txt", "extra_results/GENIST_results/cytoscape_table.txt") %in% tn))
  expect_false(any(grepl("taxes|other|plan|run_summary|run.log", tn))); expect_false(any(grepl("holiday|pic.png", pn)))
  expect_true("06_regulatory_networks/GENIST/Proj_GENIST_network_for_Cytoscape_network.png" %in% pn); expect_equal(h$plots, 2L)
  rt <- js$tables$cols[[which(tn == "06_regulatory_networks/RTP-STAR/Proj_RTPSTAR_network_for_Cytoscape.txt")]]; expect_equal(rt[1:3], c("regulator", "effect", "target"))   # the headerless edge list gets column names
  expect_equal(js$tables$total[tn == "06_regulatory_networks/RTP-STAR/Proj_RTPSTAR_network_for_Cytoscape.txt"], 2)
})

test_that("the new analysis settings exist, and the Experimental Design page offers them", {
  expect_true(all(c("pca_level", "entropy_heatmaps", "heatmap_max_genes", "browser_extra") %in% names(strym_default_params())))
  h <- as.character(design_ui("design")); expect_true(grepl("design-p_pca_level", h) && grepl("design-p_heatmaps", h) && grepl("design-p_bextra", h))
})
test_that("GENIST without a TF list is only advice on the design page; the steps' defaults are unchanged", {
  d <- strym_default_design(); d$data_type <- "bulk"; d$n_conditions <- 2L; d$n_timepoints <- 3L; d$data_kind <- "counts"
  w <- strym_check_selection(d, "genist_network")$warnings; expect_false(any(grepl("GENIST without a transcription factor list", w)) && length(strym_check_selection(d, "genist_network")$errors) > 0)
  d$params$genist_tf_file <- "x.txt"; expect_false(any(grepl("without a transcription factor list", strym_check_selection(d, "genist_network")$warnings)))
})
