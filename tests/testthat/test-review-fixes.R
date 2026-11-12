# Eleven problems found by an outside review of StRym, each reproduced first and fixed. One test per problem.
source(file.path(here, "analysis", "seurat_clustering.R"))
sc_ctx <- function(path, out, kind = "10x", ..., sample_sheet = "") {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  d <- strym_default_design(); d$data_type <- "single_cell"; d$data_kind <- kind; d$data_path <- path; d$project_dir <- dirname(out)
  a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]
  logs <- new.env(); logs$x <- character()
  list(design = strym_normalize_design(d), data_path = path, sample_sheet = sample_sheet, out_dir = out, project_dir = dirname(out), app_dir = here, prereq = list(),
       log = function(...) logs$x <- c(logs$x, paste0(...)), logs = logs)
}

test_that("1. a Seurat object without a sample column (a merged object) is analyzed, and every sample gets ITS condition", {
  skip_if_not_installed("Seurat")
  root <- tempfile("r1"); dir.create(root); syn <- synth_10x(file.path(root, "s"), cells = 40)
  objs <- lapply(seq_along(syn$dirs), function(i) Seurat::CreateSeuratObject(counts = .read_10x(syn$dirs[i]), project = syn$samples[i]))
  mg <- merge(objs[[1]], objs[-1], add.cell.ids = syn$samples); expect_false("sample" %in% colnames(mg@meta.data))
  rds <- file.path(root, "merged.rds"); saveRDS(mg, rds)
  ctx <- sc_ctx(rds, file.path(root, "o"), kind = "seurat", plots = "umap", n_conditions = 2L, params = list(sc_integration = "none"))
  r <- suppressWarnings(strym_step_seurat_clustering(ctx)); expect_equal(r$status, "ok")
  obj <- readRDS(file.path(root, "o", "seurat_object.rds")); tab <- table(obj$sample, obj$condition)
  expect_equal(rownames(tab), syn$samples); expect_equal(colnames(tab), c("Control", "Treated"))
  expect_equal(unname(tab["Control-1", ]) > 0, c(TRUE, FALSE)); expect_equal(unname(tab["Treated-2", ]) > 0, c(FALSE, TRUE))     # no cell has another sample's condition
})

test_that("2. Cell Ranger folders (Sample/outs/filtered_feature_bc_matrix) are named after the sample, and the filtered matrix wins over the raw one", {
  root <- tempfile("r2"); dir.create(root); syn <- synth_10x(file.path(root, "s"), samples = c("Control-1", "Treated-1"), cells = 20, n_genes = 100)
  mk <- function(sample, kinds) for (k in kinds) { d <- file.path(root, "cr", sample, "outs", k); dir.create(d, recursive = TRUE); file.copy(list.files(syn$dirs[1], full.names = TRUE), d) }
  mk("Control-1", c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")); mk("Treated-1", "raw_feature_bc_matrix")
  f <- .find_10x(file.path(root, "cr")); f <- f[order(f$sample), ]
  expect_equal(f$sample, c("Control-1", "Treated-1")); expect_equal(basename(f$dir), c("filtered_feature_bc_matrix", "raw_feature_bc_matrix"))   # the raw one only where it is all there is
  expect_equal(.find_10x(file.path(root, "cr", "Control-1", "outs", "filtered_feature_bc_matrix"))$sample, "Control-1")
  expect_equal(.find_10x(file.path(root, "cr", "Treated-1", "outs"))$sample, "Treated-1")
})

test_that("3. a data folder whose path contains ( ) + [ ] is searched normally", {
  root <- file.path(tempfile("r3"), "Data (copy)+1 [old]"); dir.create(root, recursive = TRUE); syn <- synth_10x(root, cells = 20, n_genes = 100)
  expect_equal(sort(.find_10x(root)$sample), sort(syn$samples))
  none <- file.path(root, "no (matrices)"); dir.create(none); expect_equal(nrow(.find_10x(none)), 0); expect_equal(names(.find_10x(none)), c("sample", "dir"))   # nothing found: an empty answer, not an error
})

test_that("4. a tab / comma file written by write.table (header one field shorter) keeps its gene ids", {
  df <- data.frame(t1 = c(1, 2, 3), t2 = c(2, 3, 4), t3 = c(5, 5, 6), row.names = c("AT1G01", "AT1G02", "AT1G03")); d <- tempfile("r4"); dir.create(d)
  for (sep in c("\t", ",")) {
    f <- file.path(d, paste0("t", if (sep == "\t") ".txt" else ".csv")); write.table(df, f, sep = sep, quote = FALSE)
    x <- read_table_any(f); expect_equal(x[[1]], rownames(df)); expect_equal(ncol(x), 4); expect_equal(x$t3, df$t3); expect_equal(read_id_column(f), rownames(df))
    m <- table_to_matrix(x); expect_equal(rownames(m), rownames(df)); expect_equal(unname(m[, "t2"]), df$t2)
  }
  g <- file.path(d, "plain.csv"); write.csv(data.frame(gene = c("A", "B"), t1 = 1:2, t2 = 3:4, t3 = 5:6), g, row.names = FALSE)                  # an ordinary table is read as before
  expect_equal(names(read_table_any(g)), c("gene", "t1", "t2", "t3")); w <- file.path(d, "w.csv"); write.csv(data.frame(t1 = 1:2, t2 = 3:4), w)   # write.csv: header has an empty first name
  expect_equal(ncol(read_table_any(w)), 3)
})

# a small treated-vs-control time course for the entropy steps
tc_toy <- function(pairs, heat_var = list(1:2), lfc_rows = list(1:2), tps = c("t1", "t2", "t3")) {
  genes <- paste0("G", 1:6); mk <- function(v = 10) matrix(v, 6, 3, dimnames = list(genes, tps)); n <- nrow(pairs)
  heat <- lapply(seq_len(n), function(i) { m <- mk(); m[heat_var[[i]], ] <- rep(c(1, 20, 20), each = length(heat_var[[i]])); m })
  lfc <- lapply(seq_len(n), function(i) { m <- mk(0); m[lfc_rows[[i]], 2:3] <- 3; m }); q <- lapply(seq_len(n), function(i) { m <- mk(1); m[lfc_rows[[i]], 2:3] <- 0.001; m })
  list(pairs = pairs, tps = tps, heat = heat, ctrl = rep(list(mk()), n), lfc = lfc, qval = q, universe = genes, genes = genes, gene_names = NULL, unit = "FPKM", samples = NULL)
}

test_that("5. Shannon entropy and MSE write empty tables, not an error, when no gene is called", {
  tc <- tc_toy(data.frame(treated = "Heat", control = "Control", stringsAsFactors = FALSE)); tc$lfc[[1]][] <- 0; tc$qval[[1]][] <- 1; p <- strym_normalize_params(list())
  for (run in list(strym_run_shannon, strym_run_mse)) {
    r <- run(tc, p); expect_equal(names(r$results), "Heat"); t <- r$results$Heat
    expect_equal(nrow(t), 0); expect_true(all(c("gene", "gene_name", "comparison", "n_sig_timepoints", "main_set", "direction", "call_t1", "log2FC_t3", "q_t2") %in% names(t)))
    out <- tempfile("r5"); dir.create(out); f <- strym_write_calls(out, "x", r, tc, "test"); expect_true(all(file.exists(file.path(out, f))))
    expect_equal(read.csv(file.path(out, "x_comparisons_used.csv"))$genes_called, 0)
  }
})

test_that("6. two comparisons with the same treated condition keep separate results and files", {
  pairs <- data.frame(treated = c("Heat", "Heat"), control = c("Control", "Mock"), stringsAsFactors = FALSE)
  expect_equal(.pair_keys(pairs), c("Heat_vs_Control", "Heat_vs_Mock")); expect_equal(.pair_keys(data.frame(treated = c("Heat", "Cold"), control = "Control")), c("Heat", "Cold"))   # unique treated: names as before
  tc <- tc_toy(pairs, heat_var = list(1:2, 4:5), lfc_rows = list(1:2, 4:5)); p <- strym_normalize_params(list()); r <- strym_run_shannon(tc, p)
  expect_equal(names(r$results), c("Heat_vs_Control", "Heat_vs_Mock")); expect_equal(sort(r$results$Heat_vs_Control$gene), c("G1", "G2")); expect_equal(sort(r$results$Heat_vs_Mock$gene), c("G4", "G5"))
  expect_equal(unique(r$results$Heat_vs_Mock$comparison), "Heat_vs_Mock"); expect_equal(sort(r$main$gene), c("G1", "G2", "G4", "G5"))
  out <- tempfile("r6"); dir.create(out); f <- strym_write_calls(out, "shannon", r, tc, "Shannon entropy")
  expect_true(all(c("shannon_Heat_vs_Control.csv", "shannon_Heat_vs_Mock.csv") %in% f)); expect_equal(read.csv(file.path(out, "shannon_comparisons_used.csv"))$genes_called, c(2, 2))
  one <- tc_toy(pairs[1, ]); expect_equal(names(strym_run_shannon(one, p)$results), "Heat")                                                      # a single pair: the old file name
  m <- strym_run_mse(tc_toy(pairs, heat_var = list(1:2, 4:5), lfc_rows = list(1:2, 4:5)), p); expect_equal(names(m$results), c("Heat_vs_Control", "Heat_vs_Mock"))
})

test_that("7. the same FASTQ name in different condition folders gives different samples, files and conditions", {
  fq <- tempfile("r7"); for (g in c("ctrl", "heat")) { dir.create(file.path(fq, g), recursive = TRUE); for (r in 1:2) writeLines(c("@r", "ACGT", "+", "IIII"), file.path(fq, g, sprintf("rep%d.fastq", r))) }
  s <- infer_samples(fq); expect_equal(s$sample, c("ctrl_rep1", "ctrl_rep2", "heat_rep1", "heat_rep2")); expect_false(anyDuplicated(s$sample) > 0)
  flat <- tempfile("r7f"); dir.create(flat); for (n in c("A1", "A2", "B1")) writeLines(c("@r", "ACGT", "+", "IIII"), file.path(flat, paste0(n, ".fastq")))
  expect_equal(infer_samples(flat)$sample, c("A1", "A2", "B1"))                                                                                    # no clash: names unchanged
  expect_equal(.unique_sample_names(c("a", "a", "b"), c("g1", "g2", "g1")), c("g1_a", "g2_a", "b")); expect_equal(.unique_sample_names(c("a", "a"), c("g", "g")), c("g_a", "g_a_1"))
  d <- strym_default_design(); d$data_type <- "bulk"; d$data_kind <- "fastq"; d$data_path <- fq; d$read_type <- "single"; d$n_conditions <- 2L
  a <- strym_auto_sample_sheet(d)$sheet; expect_equal(a$sample, s$sample); expect_equal(a$file, c("ctrl/rep1.fastq", "ctrl/rep2.fastq", "heat/rep1.fastq", "heat/rep2.fastq"))
  ss <- file.path(fq, "sheet.csv"); write.csv(a, ss, row.names = FALSE); d$sample_sheet <- ss
  q <- strym_quantify_samples(d); expect_equal(sub(paste0(fq, "/"), "", q$file, fixed = TRUE), a$file); expect_equal(q$group, c("ctrl", "ctrl", "heat", "heat"))   # each sample finds ITS file
})

test_that("8. text and annotation columns of a counts table (featureCounts: Chr, Start, End, Strand, Length; gene_name) are not samples", {
  set.seed(1); cnt <- data.frame(Geneid = paste0("G", 1:6), Chr = "1", Start = 100 * (1:6), End = 100 * (1:6) + 50, Strand = "+", Length = 50 + 10 * (1:6), gene_name = paste0("n", 1:6),
                                 c1 = rpois(6, 50), c2 = rpois(6, 50), h1 = rpois(6, 80), h2 = rpois(6, 80))
  f <- tempfile(fileext = ".csv"); write.csv(cnt, f, row.names = FALSE)
  msgs <- testthat::capture_messages(r <- .counts_raw(f)); expect_equal(colnames(r$m), c("c1", "c2", "h1", "h2")); expect_match(msgs, "left out 6 column\\(s\\).*Chr, Start, End, Strand, Length, gene_name")
  expect_equal(r$samples$sample, c("c1", "c2", "h1", "h2")); expect_equal(unname(r$m["G3", ]), unname(unlist(cnt[3, c("c1", "c2", "h1", "h2")])))
  expect_equal(strym_counts_samples(f), c("c1", "c2", "h1", "h2"))                                                                                 # the draft sample sheet lists only real samples
  cnt$weird_notes <- "see paper"; g <- tempfile(fileext = ".csv"); write.csv(cnt, g, row.names = FALSE); expect_equal(colnames(suppressMessages(.counts_raw(g))$m), c("c1", "c2", "h1", "h2"))   # any text column
  ss <- tempfile(fileext = ".csv"); write.csv(data.frame(sample = c("c1", "c2", "Length"), condition = "x"), ss, row.names = FALSE)                  # a sheet that names it keeps it
  expect_true("Length" %in% colnames(suppressMessages(.counts_raw(f, ss))$m))
  e <- tempfile(fileext = ".csv"); write.csv(data.frame(Geneid = "G1", Chr = "1", gene_name = "x"), e, row.names = FALSE); expect_error(.counts_raw(e), "no numeric sample columns")
})

test_that("9. integration on small samples falls back instead of crashing (CCA / RPCA need at least 30 cells per sample unless limited)", {
  skip_if_not_installed("Seurat"); skip_if_not_installed("harmony")
  root <- tempfile("r9"); dir.create(root); syn <- synth_10x(root, cells = 22)
  for (meth in c("cca", "rpca")) {
    ctx <- sc_ctx(syn$root, file.path(root, paste0("o", meth)), plots = "umap", n_conditions = 2L, params = list(sc_integration = meth, sc_min_counts = 50, sc_min_features = 5))
    r <- suppressWarnings(strym_step_seurat_clustering(ctx)); expect_equal(r$status, "ok", info = meth)
    expect_true(any(grepl("Integrated 4 samples with (CCA|RPCA|HARMONY)", ctx$logs$x)) || any(grepl("only merged", ctx$logs$x)), info = meth)
  }
})

test_that("10. a symbol / description column in a time-course file is not a time point", {
  set.seed(2); G <- paste0("AT1G", sprintf("%05d", 1:30)); root <- tempfile("r10"); dir.create(root)
  tcm <- data.frame(gene = G, symbol = paste0("SYM", 1:30), t1 = rpois(30, 10), t2 = rpois(30, 12), t3 = rpois(30, 20))
  tcf <- file.path(root, "tc.csv"); write.csv(tcm, tcf, row.names = FALSE); gf <- file.path(root, "genes.csv"); write.csv(data.frame(gene = G), gf, row.names = FALSE)
  m <- table_to_matrix(read_table_any(tcf)); expect_equal(colnames(m), c("t1", "t2", "t3")); expect_equal(attr(m, "dropped"), "symbol"); expect_false(anyNA(m))
  msgs <- character(); withCallingHandlers(res <- genist_run(gf, tcf, output_folder = file.path(root, "out"), cores = 1), message = function(x) { msgs <<- c(msgs, conditionMessage(x)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("time course file has non-numeric column\\(s\\) \\(symbol\\).*ignored", msgs))); expect_true(file.exists(file.path(root, "out", "cytoscape_table.txt")))
  two <- file.path(root, "two.csv"); write.csv(tcm[, c("gene", "symbol", "t1", "t2")], two, row.names = FALSE)                                   # symbol + 2 numbers = only 2 time points
  expect_error(suppressMessages(genist_run(gf, two, output_folder = file.path(root, "out2"), cores = 1)), "at least 3 time points")
})

test_that("11. results can never overwrite an input file (GENIST and RTP-STAR)", {
  root <- tempfile("r11"); dir.create(root); gf <- file.path(root, "genes.csv"); tcf <- file.path(root, "tc.csv"); cl <- file.path(root, "clusters.xlsx")
  for (f in c(gf, tcf)) writeLines("x", f); file.create(cl)
  p <- genist_check_inputs(gf, tcf, "", "", "", root, "genes.csv"); expect_length(p, 1); expect_match(p, "would overwrite your genes file \\(genes\\.csv\\).*never changes your input files")
  expect_match(genist_check_inputs(gf, tcf, "", "", "", root, "tc.csv"), "time course file \\(tc\\.csv\\)")
  expect_length(genist_check_inputs(gf, tcf, "", "", "", file.path(root, "results"), "genes.csv"), 0)                                             # another folder: fine
  expect_length(genist_check_inputs(gf, tcf, "", "", "", root, "cytoscape_table.txt"), 0)
  expect_match(genist_check_inputs(gf, tcf, cl, "", "", root, "net.txt"), "clustering data file \\(clusters\\.xlsx\\)")                           # the fixed name clusters.xlsx counts too
  expect_error(genist_run(gf, tcf, output_folder = root, output_file = "genes.csv"), "would overwrite your genes file"); expect_equal(readLines(gf), "x")
  expect_true(same_file("/mnt/d/Data/Genes.csv", "/mnt/d/data/genes.csv")); expect_false(same_file("/tmp/A.csv", "/tmp/a.csv")); expect_false(same_file("", "/tmp/a"))   # Windows drives ignore case
  expect_match(rtpstar_check_inputs(gf, tcf, output_folder = root, output_file = "genes.csv"), "genes file \\(genes\\.csv\\)")
  expect_length(rtpstar_check_inputs(gf, tcf, output_folder = file.path(root, "r"), output_file = "genes.csv"), 0)
  file.create(file.path(root, "clusters.csv")); expect_match(rtpstar_check_inputs(gf, tcf, clustering_file = file.path(root, "clusters.csv"), output_folder = root, output_file = "b.txt"), "clustering data file \\(clusters\\.csv\\)")
})
