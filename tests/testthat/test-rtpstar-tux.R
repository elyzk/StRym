# RTP-STAR inputs (a genes file and an expression file) made from Tux (Cuffdiff) results.
suppressPackageStartupMessages(library(DT)); source(file.path(here, "mod_genist.R")); source(file.path(here, "mod_rtpstar.R"))

# A small made-up Tux run in StRym's layout: two groups (ctrl, treat) x 3 replicates, 30 genes.
#   GENE01-05 up in treat and significant  -> chosen          GENE06 significant but hardly expressed -> not chosen (FPKM filter)
#   GENE07 large change, q = 0.2           -> not chosen      GENE08 significant, |log2FC| = 0.3      -> not chosen (small change)
#   GENE09 status NOTEST                   -> not chosen      locus without a name ("-")              -> cannot be named, left out
#   a second locus named GENE01, less expressed               -> only the most expressed one is kept
tux_rtp_fixture <- function(named = TRUE, reps = 3, groups = c("ctrl", "treat"), write_rg = TRUE, de_genes = 1:5) {
  root <- tempfile("rtpt"); td <- file.path(root, "01_quantification", "tux"); dd <- file.path(td, "diff_out"); dir.create(dd, recursive = TRUE)
  set.seed(4)
  ng <- 30; ids <- sprintf("XLOC_%03d", seq_len(ng)); nm <- if (named) sprintf("GENE%02d", seq_len(ng)) else rep("-", ng)
  if (named) nm[20] <- "-"                                                # the unnamed locus
  mu <- exp(rnorm(ng, log(20), 0.5)); mu[6] <- 0.2
  rows <- list()
  for (g in groups) for (r in seq_len(reps)) {
    eff <- rep(1, ng); if (g == groups[2]) eff[de_genes] <- 6
    f <- mu * eff * exp(rnorm(ng, 0, 0.1)); f[6] <- 0.3
    rows[[length(rows) + 1]] <- data.frame(tracking_id = ids, condition = g, replicate = r - 1L, raw_frags = round(f * 10), internal_scaled_frags = round(f * 10),
                                           external_scaled_frags = round(f * 10), FPKM = f, effective_length = 1000, status = "OK", stringsAsFactors = FALSE)
  }
  RG <- do.call(rbind, rows)
  if (named) { dup <- RG[RG$tracking_id == ids[1], ]; dup$tracking_id <- "XLOC_DUP"; dup$FPKM <- 1; RG <- rbind(RG, dup) }
  wr <- function(df, f) utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  if (write_rg) wr(RG, file.path(dd, "genes.read_group_tracking"))
  tr <- data.frame(tracking_id = ids, class_code = "-", nearest_ref_id = "-", gene_id = ids, gene_short_name = nm, tss_id = "-", locus = "1:1-2", length = "-", coverage = "-", stringsAsFactors = FALSE)
  for (g in groups) tr[[paste0(g, "_FPKM")]] <- tapply(RG$FPKM[RG$condition == g & RG$tracking_id %in% ids], RG$tracking_id[RG$condition == g & RG$tracking_id %in% ids], mean)[ids]
  if (named) tr <- rbind(tr, transform(tr[1, ], tracking_id = "XLOC_DUP", gene_id = "XLOC_DUP"))
  wr(tr, file.path(dd, "genes.fpkm_tracking"))
  q <- rep(0.9, ng); q[c(de_genes, 6, 8, 9, 20)] <- 0.001; q[7] <- 0.2
  lfc <- rep(0.05, ng); lfc[de_genes] <- 2.5; lfc[6] <- 3; lfc[7] <- 3; lfc[8] <- 0.3; lfc[9] <- 3; lfc[20] <- 3
  st <- rep("OK", ng); st[9] <- "NOTEST"
  de <- data.frame(test_id = ids, gene_id = ids, gene = nm, locus = "1:1-2", sample_1 = groups[1], sample_2 = groups[2], status = st, value_1 = 1, value_2 = 1,
                   `log2(fold_change)` = lfc, test_stat = 0, p_value = 0.01, q_value = q, significant = "yes", check.names = FALSE, stringsAsFactors = FALSE)
  wr(de, file.path(dd, "gene_exp.diff"))
  used <- do.call(rbind, lapply(groups, function(g) data.frame(group = g, k = seq_len(reps), file = sprintf("/data/%s_%d.fastq.gz", g, seq_len(reps)), mate = NA, sample = sprintf("%s_s%d", g, seq_len(reps)))))
  utils::write.csv(used, file.path(td, "samples_used.csv"), row.names = FALSE)
  list(root = root, tux = td, diff = dd, out = file.path(root, "rtp_in"), genes = nm, used = used)
}

test_that("the genes file and the expression file are made from the Tux results, ready for RTP-STAR", {
  fx <- tux_rtp_fixture(); r <- tux_rtpstar_inputs(fx$diff, fx$out)
  g <- read.csv(r$genes_file); expect_equal(names(g), "gene"); expect_equal(g$gene, sprintf("GENE%02d", 1:5))             # the planted genes only
  ex <- read.csv(r$expression_file, check.names = FALSE)
  expect_equal(names(ex)[1], "gene"); expect_equal(names(ex)[-1], c("ctrl_s1", "ctrl_s2", "ctrl_s3", "treat_s1", "treat_s2", "treat_s3"))   # one column per sample, the names of samples_used.csv
  expect_true(all(g$gene %in% ex$gene)); expect_false(any(duplicated(ex$gene))); expect_false("-" %in% ex$gene)
  expect_false("GENE06" %in% ex$gene)                                                                                       # hardly expressed anywhere: not in the file
  # the values are the replicate FPKMs, and the duplicated name keeps its most expressed locus
  rg <- read.delim(file.path(fx$diff, "genes.read_group_tracking")); one <- rg[rg$tracking_id == "XLOC_001" & rg$condition == "treat" & rg$replicate == 1, "FPKM"]
  expect_equal(ex[ex$gene == "GENE01", "treat_s2"], round(one, 4)); expect_true(all(ex[ex$gene == "GENE01", -1] > 5))
  expect_true(all(ex[ex$gene %in% g$gene, "treat_s1"] > 3 * ex[ex$gene %in% g$gene, "ctrl_s1"] * 0.8))
  expect_true(file.exists(r$readme)); expect_match(r$summary, "5 differentially expressed genes"); expect_equal(r$n_genes, 5); expect_equal(r$samples, names(ex)[-1])
  expect_length(r$warnings, 0)
})

test_that("the cutoffs, the cap and the gene names behave", {
  fx <- tux_rtp_fixture()
  expect_equal(read.csv(tux_rtpstar_inputs(fx$diff, fx$out, min_log2fc = 0.2)$genes_file)$gene, sprintf("GENE%02d", c(1:5, 8)))            # a smaller change lets GENE08 in
  expect_equal(read.csv(tux_rtpstar_inputs(fx$diff, fx$out, q_cutoff = 0.5)$genes_file)$gene, sprintf("GENE%02d", c(1:5, 7)))              # a looser q lets GENE07 in
  expect_true("GENE06" %in% read.csv(tux_rtpstar_inputs(fx$diff, fx$out, min_fpkm = 0)$genes_file)$gene)
  r2 <- tux_rtpstar_inputs(fx$diff, fx$out, max_genes = 3); expect_true(r2$capped); expect_equal(r2$n_pass, 5); expect_length(read.csv(r2$genes_file)$gene, 3)
  expect_false(tux_rtpstar_inputs(fx$diff, fx$out, max_genes = 0)$capped)
  nn <- tux_rtp_fixture(named = FALSE); rn <- tux_rtpstar_inputs(nn$diff, nn$out)                                                         # an annotation without gene names: the Tux ids
  expect_equal(read.csv(rn$genes_file)$gene, sprintf("XLOC_%03d", c(1:5, 20))); expect_match(readLines(rn$readme)[grepl("named", readLines(rn$readme))][1], "Tux ids")
})

test_that("results that cannot give useful files say why", {
  fx <- tux_rtp_fixture()
  expect_error(tux_rtpstar_inputs(fx$diff, fx$out, q_cutoff = 0), "between 0 and 1"); expect_error(tux_rtpstar_inputs(fx$diff, fx$out, min_log2fc = -1), "0 or more")
  expect_error(tux_rtpstar_inputs(fx$diff, fx$out, q_cutoff = 1e-9), "Only 0 gene"); expect_error(tux_rtpstar_inputs(fx$diff, fx$out, min_log2fc = 10), "Only 0 gene")
  few <- tux_rtp_fixture(reps = 1); expect_error(tux_rtpstar_inputs(few$diff, few$out), "needs several samples")                            # 2 samples
  one <- tux_rtp_fixture(de_genes = 1); expect_error(tux_rtpstar_inputs(one$diff, one$out), "Only 1 gene")
  nr <- tux_rtp_fixture(write_rg = FALSE); expect_error(tux_rtpstar_inputs(nr$diff, nr$out), "genes.read_group_tracking")
  expect_error(tux_rtpstar_inputs(tempdir(), fx$out), "gene_exp.diff")
  single <- tux_rtp_fixture(reps = 1, groups = c("a", "b", "c")); w <- tux_rtpstar_inputs(single$diff, single$out)$warnings                    # 3 samples, one replicate each: allowed, with a note
  expect_true(any(grepl("single replicate", w))); expect_true(any(grepl("Only 3 samples", w)))
})

test_that("a Tux run leaves the files next to its results, or a note; a rerun never leaves stale files", {
  fx <- tux_rtp_fixture(); out <- fx$tux                                                       # the Tux output folder holds diff_out/
  r <- tux_write_rtpstar_inputs(out); d <- file.path(out, "rtpstar_inputs")
  expect_true(all(file.exists(file.path(d, c("rtpstar_genes_from_tux.csv", "rtpstar_expression_from_tux.csv", "README_rtpstar_inputs.txt"))))); expect_false(is.null(r))
  expect_false(file.exists(file.path(d, "rtpstar_inputs_NOT_made.txt")))
  # the same folder, but the results now give too few genes: the old files must go, a note must say why, and nothing may raise an error
  expect_error(tux_write_rtpstar_inputs(out, q_cutoff = 1e-9), NA)
  expect_false(file.exists(file.path(d, "rtpstar_genes_from_tux.csv"))); expect_true(file.exists(file.path(d, "rtpstar_inputs_NOT_made.txt")))
  expect_match(paste(readLines(file.path(d, "rtpstar_inputs_NOT_made.txt")), collapse = " "), "Only 0 gene")
  expect_error(tux_write_rtpstar_inputs(tempfile("nothing")), NA)                              # no results at all: still no error
  expect_true("rtpstar_inputs" %in% .OUTPUT_DIRS)                                              # cleared with the other outputs
})

test_that("the files go straight into RTP-STAR and give a network (project folder, tux folder or diff_out as the starting point)", {
  fx <- tux_rtp_fixture()
  for (p in c(fx$root, fx$tux, fx$diff)) { r <- rtpstar_inputs_from_tux(p, file.path(fx$out, basename(p))); expect_equal(r$n_genes, 5); expect_true(file.exists(r$expression_file)) }
  expect_equal(rtpstar_inputs_from_tux(fx$root, fx$out)$project_dir, NULL)                      # (no 00_design here: not a StRym project record)
  r <- rtpstar_inputs_from_tux(fx$root, fx$out)
  outn <- file.path(fx$root, "net"); rtpstar_run(r$genes_file, r$expression_file, outn, output_file = "net.txt", ntrees = 100)
  expect_true(file.exists(file.path(outn, "net.txt"))); l <- readLines(file.path(outn, "net.txt")); expect_true(all(grepl("^GENE[0-9]+ (activates|inhibits|regulates) GENE[0-9]+", l)))
  expect_length(rtpstar_check_inputs(r$genes_file, r$expression_file, output_folder = outn, output_file = "net.txt"), 0)              # the checks of the RTP-STAR tab pass on them
})

test_that("the RTP-STAR tab offers to start from Tux results", {
  h <- as.character(rtpstar_ui("rtp")); expect_true(grepl("Start from your Tux results", h)); expect_true(grepl("rtp-tux-txt", h)); expect_true(grepl("tux_opts", h))
  expect_true("from_tux" %in% names(formals(rtpstar_server)))
})

test_that("the RTP-STAR results of Tux-made inputs go to the project's network folder, or beside the inputs, never inside them", {
  fx <- tux_rtp_fixture(); r <- tux_rtpstar_inputs(fx$diff, file.path(fx$tux, "rtpstar_inputs"))
  expect_equal(basename(rtp_default_out(r$genes_file)), "RTP-STAR_results"); expect_equal(dirname(rtp_default_out(r$genes_file)), normalizePath(fx$tux, winslash = "/"))
  dir.create(file.path(fx$root, "00_design")); writeLines("{}", file.path(fx$root, "00_design", "design.json"))               # now it is a StRym project
  expect_equal(rtp_default_out(r$genes_file), file.path(normalizePath(fx$root, winslash = "/"), "06_regulatory_networks", "RTP-STAR"))
})