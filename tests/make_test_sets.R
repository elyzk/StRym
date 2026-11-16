#!/usr/bin/env Rscript
# make_test_sets.R -- builds the version 9/10 test-case folders (made-up data, deterministic for a given seed).
#
#   Rscript tests/make_test_sets.R <output folder> [<StRym folder>]
#
# Creates five folders, each with a README.txt (what to do, what a correct run shows) and the input files:
#   V10_1_drive_reads_tux_to_rtpstar     FASTQ + GFF3 annotation on a mounted drive: Tux, then RTP-STAR from the Tux results
#   V10_2_drive_timecourse_clutter       time-course counts on a mounted drive; output folder that already holds other files; sample-sheet
#                                        changes before a run; PCA levels; entropy heatmaps; GENIST with and without a TF list
#   V10_3_heatmap_pca_edge_cases         three tables (tiny, wide, many conditions) for the heatmaps and the PCA levels
#   V10_4_cuffdiff_single_file_modes     one gene_exp.diff file as input; GENIST and RTP-STAR tabs started from Tux results; preset clusters
#   V10_5_network_tabs_browser           GENIST and RTP-STAR tabs on their own files: file names, pictures, optional TF list, extra result
#                                        folders in the interactive browser
# Base R only (the reads come from make_synthetic.R). Folders 1 and 2 are meant to be copied to a USB stick or external drive.

argv <- commandArgs(TRUE)
if (!length(argv)) stop("usage: Rscript tests/make_test_sets.R <output folder> [<StRym folder>]")
out <- normalizePath(argv[1], mustWork = FALSE); dir.create(out, recursive = TRUE, showWarnings = FALSE)
app <- if (length(argv) >= 2) argv[2] else getwd()
set.seed(20261110)

# ---- helpers --------------------------------------------------------------------------------------------------------------------------------
wr <- function(lines, f) writeLines(lines, f, useBytes = TRUE)
wcsv <- function(df, f, ...) utils::write.csv(df, f, row.names = FALSE, ...)

# Counts for a time-course experiment. Column names: <condition>_<time point>_rep<k>. The first condition is the control. In every other condition a
# share of the genes responds after the first time point (up, down, or a transient peak); everything else is flat. Returns list(counts, sheet, truth).
sim_counts <- function(n_genes, conds, tps, reps, responders = 0.1, depth = 80, disp = 0.08, id_prefix = "GENE") {
  ids <- sprintf("%s%05d", id_prefix, seq_len(n_genes)); base <- exp(rnorm(n_genes, log(depth), 1.1))
  kind <- rep("flat", n_genes); kind[sample(n_genes, round(n_genes * responders))] <- sample(c("up", "down", "peak"), round(n_genes * responders), TRUE)
  cols <- list(); sheet <- list()
  for (cn in conds) for (ti in seq_along(tps)) for (r in seq_len(reps[[cn]][ti])) {
    eff <- rep(0, n_genes)
    if (cn != conds[1] && ti > 1) { eff[kind == "up"] <- 2.2; eff[kind == "down"] <- -2.2; if (ti == 2) eff[kind == "peak"] <- 2.5 }
    mu <- base * 2^eff * 1.0
    nm <- sprintf("%s_%s_rep%d", cn, tps[ti], r)
    cols[[nm]] <- rnbinom(n_genes, mu = mu, size = 1 / disp)
    sheet[[nm]] <- data.frame(sample = nm, condition = cn, time_point = tps[ti], replicate = r, stringsAsFactors = FALSE)
  }
  counts <- data.frame(gene = ids, as.data.frame(cols, check.names = FALSE), check.names = FALSE)
  list(counts = counts, sheet = do.call(rbind, sheet), truth = data.frame(gene = ids, kind = kind))
}
reps_all <- function(conds, tps, n) setNames(lapply(conds, function(c) rep(n, length(tps))), conds)

# A synthetic Cuffdiff result (what a Tux run leaves in diff_out/) for a time course: groups <condition>_<time point>, FPKM per replicate, q-values
# from a t-test on log2(FPKM + 1) (Benjamini-Hochberg over genes), Cuffdiff's own file layout and column names.
sim_tux <- function(dir, n_genes, conds, tps, reps = 3, responders = 0.15, named = TRUE) {
  s <- sim_counts(n_genes, conds, tps, reps_all(conds, tps, reps), responders)
  M <- as.matrix(s$counts[-1]); fp <- sweep(M, 2, colSums(M) / 1e6, "/") / 2                # FPKM-like
  ids <- sprintf("XLOC_%06d", seq_len(n_genes)); gname <- if (named) sprintf("AT1G%05d", seq_len(n_genes)) else rep("-", n_genes)
  dd <- file.path(dir, "01_quantification", "tux", "diff_out"); dir.create(dd, recursive = TRUE)
  grp <- paste(s$sheet$condition, s$sheet$time_point, sep = "_"); sheet <- s$sheet; sheet$group <- grp
  used <- data.frame(group = grp, k = sheet$replicate, file = paste0(sheet$sample, ".fastq.gz"), mate = NA, sample = sheet$sample)
  wcsv(used, file.path(dir, "01_quantification", "tux", "samples_used.csv"))
  wt <- function(df, f) utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  rg <- do.call(rbind, lapply(seq_len(ncol(fp)), function(j) data.frame(tracking_id = ids, condition = grp[j], replicate = sheet$replicate[j] - 1L, raw_frags = round(M[, j]),
         internal_scaled_frags = round(M[, j]), external_scaled_frags = round(M[, j]), FPKM = round(fp[, j], 4), effective_length = 1000, status = "OK")))
  wt(rg, file.path(dd, "genes.read_group_tracking"))
  g_names <- unique(grp); tr <- data.frame(tracking_id = ids, class_code = "-", nearest_ref_id = "-", gene_id = ids, gene_short_name = gname, tss_id = "-", locus = "chr1:1-2",
                                          length = "-", coverage = "-", stringsAsFactors = FALSE)
  for (g in g_names) { tr[[paste0(g, "_FPKM")]] <- round(rowMeans(fp[, grp == g, drop = FALSE]), 4); tr[[paste0(g, "_conf_lo")]] <- 0; tr[[paste0(g, "_conf_hi")]] <- 0; tr[[paste0(g, "_status")]] <- "OK" }
  wt(tr, file.path(dd, "genes.fpkm_tracking"))
  pairs <- list(); for (ti in seq_along(tps)) pairs[[length(pairs) + 1]] <- c(paste(conds[1], tps[ti], sep = "_"), paste(conds[2], tps[ti], sep = "_"))   # control vs treated at every time point
  for (cn in conds[-1]) { a <- paste(cn, tps, sep = "_"); for (i in 1:(length(a) - 1)) for (j in (i + 1):length(a)) pairs[[length(pairs) + 1]] <- c(a[i], a[j]) }  # treated: every pair of time points
  de <- do.call(rbind, lapply(pairs, function(p) {
    x <- log2(fp[, grp == p[1], drop = FALSE] + 1); y <- log2(fp[, grp == p[2], drop = FALSE] + 1)
    pv <- vapply(seq_len(n_genes), function(i) tryCatch(stats::t.test(x[i, ], y[i, ])$p.value, error = function(e) 1), 1); pv[is.na(pv)] <- 1
    m1 <- rowMeans(fp[, grp == p[1], drop = FALSE]); m2 <- rowMeans(fp[, grp == p[2], drop = FALSE]); lf <- log2((m2 + 1e-3) / (m1 + 1e-3))
    data.frame(test_id = ids, gene_id = ids, gene = gname, locus = "chr1:1-2", sample_1 = p[1], sample_2 = p[2], status = "OK", value_1 = round(m1, 4), value_2 = round(m2, 4),
               `log2(fold_change)` = round(lf, 4), test_stat = 0, p_value = pv, q_value = stats::p.adjust(pv, "BH"), significant = "no", check.names = FALSE, stringsAsFactors = FALSE)
  }))
  de$significant <- ifelse(de$q_value < 0.05, "yes", "no")
  wt(de, file.path(dd, "gene_exp.diff"))
  list(counts = s$counts, sheet = s$sheet, fpkm = data.frame(gene = gname, as.data.frame(round(fp, 4), check.names = FALSE), check.names = FALSE), diff = de, truth = s$truth, genes = gname)
}

# GFF3 version of a GTF (genes with their exons), to test GFF3 annotation input.
gtf_to_gff3 <- function(gtf, gff3) {
  l <- readLines(gtf); l <- l[nzchar(l)]; f <- strsplit(l, "\t", fixed = TRUE)
  gid <- vapply(f, function(x) sub('.*gene_id "([^"]+)".*', "\\1", x[9]), ""); out <- "##gff-version 3"
  for (g in unique(gid)) { r <- f[gid == g]; s <- as.integer(vapply(r, `[`, "", 4)); e <- as.integer(vapply(r, `[`, "", 5))
    out <- c(out, sprintf("%s\tsynthetic\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s", r[[1]][1], min(s), max(e), r[[1]][7], g, g),
             sprintf("%s\tsynthetic\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.1;Parent=%s", r[[1]][1], min(s), max(e), r[[1]][7], g, g),
             vapply(seq_along(r), function(k) sprintf("%s\tsynthetic\texon\t%d\t%d\t.\t%s\t.\tID=%s.1.exon%d;Parent=%s.1", r[[k]][1], s[k], e[k], r[[k]][7], g, k, g), "")) }
  writeLines(out, gff3)
}

readme <- function(dir, title, body) wr(c(title, strrep("=", nchar(title)), "", body), file.path(dir, "README.txt"))

# ======================================================================================================================================================
# 1. FASTQ on a drive -> Tux -> RTP-STAR from the Tux results
# ======================================================================================================================================================
d1 <- file.path(out, "V10_1_drive_reads_tux_to_rtpstar"); dir.create(d1, recursive = TRUE, showWarnings = FALSE)
ok <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(app, "make_synthetic.R")), "--out", shQuote(file.path(d1, "_tmp")), "--seed", "11", "--layout", "flat", "--gz", "--reps", "4",
                           "--genes-per-chr", "14", "--chr", "2", "--chr-len", "70000", "--read-len", "75", "--depth", "1"), stdout = FALSE, stderr = FALSE)
if (!identical(ok, 0L) && !file.exists(file.path(d1, "_tmp", "genome.fa"))) stop("make_synthetic.R failed (give the StRym folder as the second argument)")
dir.create(file.path(d1, "reads")); file.copy(list.files(file.path(d1, "_tmp", "inputs"), full.names = TRUE), file.path(d1, "reads"))
file.copy(file.path(d1, "_tmp", c("genome.fa", "genes.gtf", "truth.tsv")), d1)
gtf_to_gff3(file.path(d1, "genes.gtf"), file.path(d1, "genes.gff3")); unlink(file.path(d1, "_tmp"), recursive = TRUE)
readme(d1, "V10_1: reads on a drive, Tux, then RTP-STAR from the Tux results (uses a mounted drive)", c(
  "WHAT THIS TESTS",
  "  - a mounted drive: the whole project (reads, annotation, results) lives on a USB stick / external drive",
  "  - GFF3 annotation instead of GTF (genes.gff3 is the same annotation as genes.gtf)",
  "  - the new step at the end of a Tux run: StRym writes a genes file and a gene expression file for RTP-STAR (tux/rtpstar_inputs/)",
  "  - the TUX results page: download buttons and 'Use them in the RTP-STAR tab'; the RTP-STAR tab: 'Start from your Tux results'", "",
  "SET UP",
  "  1. Copy this whole folder to the root of a USB stick or external drive. Plug it in.",
  "  2. Open StRym > Drives. Windows: press Refresh; if the drive says 'not mounted', select it and press Mount drive. Mac / Linux: press Refresh.",
  "     The drive must show as 'ready'. Write down the name or letter you see.", "",
  "RUN 1: Experimental Design page, raw reads",
  "  - Project name: V10_1 drive reads. Output folder: a new folder ON THE DRIVE, e.g. <drive>/V10_1_results.",
  "  - Type of data: whole organ/tissue. What do you have: raw reads (FASTQ). Input data: the 'reads' folder on the drive (8 files, 4 per condition).",
  "  - Annotation: genes.gff3 (NOT the .gtf). Genome: genome.fa. Single-end, short reads. Conditions: 2 (CTL, TRT), time points: 1.",
  "  - Tick DESeq2, Data overview, PCA, Interactive web browser. Run.",
  "  EXPECTED: the pre-run window lists 8 samples in 2 groups, with a time estimate (about 3 minutes); every step ends 'ok';",
  "  01_quantification/tux/rtpstar_inputs/ holds rtpstar_genes_from_tux.csv, rtpstar_expression_from_tux.csv and README_rtpstar_inputs.txt;",
  "  the Quantify line in the results table says 'A genes file and an expression file for RTP-STAR are in tux/rtpstar_inputs/'.",
  "  The expression file has 8 sample columns; the genes file lists the differentially expressed genes (see truth.tsv: genes marked up or down).", "",
  "RUN 2: the TUX tab on the same reads (optional)",
  "  - Input: the reads folder, output: <drive>/V10_1_tux_tab. After the run the Results tab shows the card 'Inputs for RTP-STAR' with the two files.",
  "  - Press 'Use them in the RTP-STAR tab': the page jumps to RTP-STAR with the genes file and the expression file filled in. Leave the TF file empty.",
  "    EXPECTED: a note that no TF file was given (the run takes longer; it still runs). Run it: a network table named after the genes file",
  "    (rtpstar_genes_from_tux_RTPSTAR_network_for_Cytoscape.txt) and a *_network.png picture appear in the output folder.",
  "  - In the RTP-STAR tab, under 'Start from your Tux results', choose <drive>/V10_1_tux_tab, set 'Smallest change' to 2 and press Build: fewer genes.", "",
  "UNPLUG TEST (needs the drive)",
  "  - With no run going, press Release drive (Windows) and eject the drive. Plug it in again, press Refresh (and Mount drive on Windows).",
  "  - EXPECTED: the drive is listed again as ready, results are still there. While it is unplugged, typing a path on it in any box gives a plain message.", "",
  "DEBUGGING HINTS",
  "  - If Quantify fails: read the check list shown before the run, then the run log at the bottom of the page.",
  "  - Is the GFF3 refused? Try genes.gtf; if that works and the GFF3 does not, report it (the file is valid GFF3 with gene / mRNA / exon lines)."))

# ======================================================================================================================================================
# 2. Time-course counts on a drive; output folder with other files; sample-sheet changes; PCA levels; heatmaps; GENIST with / without TF list
# ======================================================================================================================================================
d2 <- file.path(out, "V10_2_drive_timecourse_clutter"); dir.create(d2, recursive = TRUE, showWarnings = FALSE)
conds <- c("control", "heat"); tps <- c("0h", "1h", "3h", "6h", "12h")
s2 <- sim_counts(4000, conds, tps, reps_all(conds, tps, 3), responders = 0.12, id_prefix = "AT1G")
wcsv(s2$counts, file.path(d2, "counts.csv")); wcsv(s2$sheet, file.path(d2, "sample_sheet.csv")); wcsv(s2$truth, file.path(d2, "truth_responders.csv"))
tf <- data.frame(gene = s2$counts$gene[sample(4000, 300)]); wcsv(tf, file.path(d2, "TF_list.csv"))
cl <- file.path(d2, "output_folder_with_other_files"); dir.create(file.path(cl, "old photos"), recursive = TRUE)
wcsv(data.frame(item = c("rent", "food"), cost = c(900, 300)), file.path(cl, "budget.csv")); wcsv(data.frame(a = 1:3, b = 4:6), file.path(cl, "results_of_something_else.csv"))
grDevices::png(file.path(cl, "holiday_photo.png"), 400, 300); graphics::plot(1:10); grDevices::dev.off()
grDevices::png(file.path(cl, "old photos", "cat.png"), 400, 300); graphics::plot(10:1); grDevices::dev.off()
wr("Not a StRym file.", file.path(cl, "notes.txt")); wr(c("source\ttarget\twidth\tsign", "X\tY\t1\t1"), file.path(cl, "some_network.txt"))
readme(d2, "V10_2: time-course counts on a drive, an output folder with other files in it (uses a mounted drive)", c(
  "WHAT THIS TESTS",
  "  - a mounted drive for input AND output; an output folder that already holds the researcher's own files (budget.csv, photos, ...): StRym must",
  "    show only its own results in the Results list and in the interactive browser, and must never touch the other files",
  "  - the sample sheet shown before the run always matches the entries on the page (names, numbers of conditions and time points)",
  "  - PCA level (replicate / averaged / both); heatmaps of the entropy results; GENIST with and without a TF list; network file names and picture",
  "  - 4000 genes, 2 conditions (control, heat) x 5 time points (0h 1h 3h 6h 12h) x 3 replicates; 12% of the genes respond to heat (truth_responders.csv)", "",
  "SET UP",
  "  1. Copy this folder to a USB stick / external drive and plug it in. Drives tab: Refresh (Windows: Mount drive if it says not mounted).", "",
  "RUN 1",
  "  - Experimental Design: project 'V10_2 drive', input counts.csv on the drive, sample sheet blank (StRym builds it).",
  "  - Conditions 2, time points 5. Condition names: control, heat. Time-point names: 0h, 1h, 3h, 6h, 12h.",
  "  - Output folder: <drive>/V10_2_drive_timecourse_clutter/output_folder_with_other_files  (the folder that already holds budget.csv and the photos).",
  "  - Tick DESeq2, Shannon entropy, Modified Shannon entropy, Data overview, PCA, Trajectories, GENIST, the interactive browser.",
  "  - Analysis settings: 'PCA: what is one point?' = Both; leave the heatmap box ticked. Leave the TF list empty and set 'Most genes in the network' to 120",
  "    (GENIST without a TF list on all ~300 genes the entropy steps find would take 15 minutes or more; with 120 it takes about a minute).",
  "  EXPECTED before running: a warning that GENIST without a TF list takes longer (not an error); the sample sheet table shows 30 samples;",
  "  the confirmation window shows 2 conditions x 5 time points with 3 samples each.",
  "  EXPECTED after: every step ok. Results list: only plots from 00_data_overview, 02_..., 04_pca, 05_trajectories, 06_regulatory_networks;",
  "  NO holiday_photo.png, no old photos/cat.png. 04_pca has pca_plot.png (30 points) and pca_averaged_plot.png (10 points: 2 conditions x 5 time points).",
  "  02_differential_expression/shannon_entropy and modified_shannon_entropy each hold *_heatmap_log2FC_heat.png and *_heatmap_expression_heat.png.",
  "  06_regulatory_networks/GENIST_auto/heat/results/ holds V10_2_drive_heat_GENIST_network_for_Cytoscape.txt (the project name is in the file name)",
  "  and a *_network.png picture. Open 07_interactive_browser/strym_browser.html: Tables must list the GENIST network table, Plots the network picture,",
  "  and NOTHING from budget.csv / some_network.txt / the photos. The other files in the output folder are unchanged.", "",
  "RUN 2: the sample sheet always matches the page",
  "  - Press 'Clear all fields', then enter the same experiment again but type Condition names 'ctrl, hs' (names that do not appear in the column names).",
  "  - EXPECTED: the sample sheet and the confirmation window describe what StRym could match; they never show the previous run's 'control / heat' labels.",
  "  - Then change the number of time points to 4 without touching the names and open the confirmation window again: it shows 4 time points.",
  "  - Edit one cell of the sample sheet by hand, then change the data file: a warning appears if the edited sheet no longer matches the data.", "",
  "RUN 3: GENIST with the TF list",
  "  - Add TF_list.csv under Analysis settings > Transcription factor list. EXPECTED: no 'takes longer' warning; the same steps run, faster GENIST.", "",
  "UNPLUG TEST",
  "  - Close the browser tab during a run (the run continues), reopen it, and see the progress. After it finishes: Release drive, eject, plug in, Refresh.",
  "  - EXPECTED: results are still on the drive; choosing the same output folder shows the earlier results in the Results list.", "",
  "DEBUGGING HINTS",
  "  - If the Results list shows files that StRym did not write, the filter strym_result_files() (pipeline.R) is the place to look.",
  "  - If the sample sheet shows old labels, look at refresh_sheet() and the sheet_key in mod_design.R."))
cat("1 and 2 written\n")

# ======================================================================================================================================================
# 3. Heatmap and PCA edge cases (three tables)
# ======================================================================================================================================================
d3 <- file.path(out, "V10_3_heatmap_pca_edge_cases"); dir.create(d3, recursive = TRUE, showWarnings = FALSE)
a <- file.path(d3, "A_tiny_few_genes"); dir.create(a)
sa <- sim_counts(60, c("control", "treated"), c("0h", "2h", "6h"), reps_all(c("control", "treated"), c("0h", "2h", "6h"), 2), responders = 0.1, depth = 150)
wcsv(sa$counts, file.path(a, "counts.csv")); wcsv(sa$sheet, file.path(a, "sample_sheet.csv"))
b <- file.path(d3, "B_wide_20000_genes"); dir.create(b)
sb <- sim_counts(20000, c("control", "heat"), c("0h", "1h", "3h", "6h"), reps_all(c("control", "heat"), c("0h", "1h", "3h", "6h"), 3), responders = 0.2)
wcsv(sb$counts, file.path(b, "counts.csv")); wcsv(sb$sheet, file.path(b, "sample_sheet.csv"))
cc <- file.path(d3, "C_many_conditions"); dir.create(cc)
cn <- c("control", paste0("treatment", 1:5)); tp <- c("0h", "1h", "4h", "8h")
sc <- sim_counts(3000, cn, tp, reps_all(cn, tp, 2), responders = 0.15)
wcsv(sc$counts, file.path(cc, "counts.csv")); wcsv(sc$sheet, file.path(cc, "sample_sheet.csv"))
dd3 <- file.path(d3, "D_only_one_gene_responds"); dir.create(dd3)
sd3 <- sim_counts(500, c("control", "treated"), c("0h", "1h", "4h"), reps_all(c("control", "treated"), c("0h", "1h", "4h"), 3), responders = 0.0)
sd3$counts[1, grep("^treated_(1h|4h)", names(sd3$counts))] <- sd3$counts[1, grep("^treated_(1h|4h)", names(sd3$counts))] * 8L + 400L; sd3$counts[1, 1] <- "ONE_RESPONDER"
wcsv(sd3$counts, file.path(dd3, "counts.csv")); wcsv(sd3$sheet, file.path(dd3, "sample_sheet.csv"))
readme(d3, "V10_3: entropy heatmaps and PCA levels on awkward tables", c(
  "WHAT THIS TESTS",
  "  - the heatmaps (pheatmap) of the genes the entropy steps call, from a few genes to many genes and conditions; the PCA level setting",
  "  - four tables, each with its own sample_sheet.csv (use it, or leave the sheet blank: StRym reads the column names)", "",
  "FOR EVERY TABLE",
  "  - Experimental Design, counts table, tick DESeq2 + Shannon + Modified Shannon + PCA + Data overview. Analysis settings: PCA = Both. Heatmap box ticked.",
  "  - Look at 04_pca and at the heatmap pictures in 02_differential_expression/shannon_entropy and .../modified_shannon_entropy.", "",
  "A_tiny_few_genes (60 genes, 2 conditions x 3 time points x 2 replicates)",
  "  EXPECTED: runs; maybe only a handful of genes called. If at least 2 genes are called a heatmap is drawn per comparison; with 0-1 called genes the",
  "  step still finishes 'ok' (no heatmap or a one-row heatmap, never an error). PCA: 12 points at replicate level, 6 averaged points.", "",
  "B_wide_20000_genes (20 000 genes, 2 conditions x 4 time points x 3 replicates)",
  "  EXPECTED: thousands of genes are called but each heatmap shows at most 100 (the setting 'Heatmap: most genes per comparison'); the figure has no gene names",
  "  (names are printed only up to 80 genes) and is tall but readable. Time: a few minutes. Try 'Heatmap: most genes' = 30: names appear, the figure is short.", "",
  "C_many_conditions (3000 genes, control + 5 treatments x 4 time points x 2 replicates = 48 samples)",
  "  EXPECTED: five treated-vs-control comparisons; each gets its own two heatmaps, plus one '..._heatmap_log2FC_all_comparisons.png' with all comparisons side by side",
  "  (20 columns, separated by gaps, a bar on top naming the comparison). PCA averaged: 24 points (6 conditions x 4 time points); replicate level: 48.", "",
  "D_only_one_gene_responds (500 genes, one gene, ONE_RESPONDER, strongly up in treated after 0h)",
  "  EXPECTED: no crash. With the default cutoffs the main set holds 1 gene (or 0 if the cutoffs miss it): the heatmaps show one row or are not made; GENIST is skipped",
  "  with a plain message (not enough genes for a network).", "",
  "SWITCHES TO TRY ON ANY TABLE",
  "  - untick 'Draw a heatmap of the genes the entropy steps call': no *_heatmap_* files.",
  "  - PCA = Replicate only: files are the same names as in earlier versions (pca_plot.png, pca_scores.csv, pca_variance.csv).",
  "  - PCA = Averaged with a table that has only 2 groups (take any table and set 1 time point): a clear message that 3 groups are needed.", "",
  "DEBUGGING HINTS",
  "  - Heatmap code: strym_entropy_heatmaps() and .hm_draw() in analysis/common.R. PCA: analysis/pca.R."))
cat("3 written\n")

# ======================================================================================================================================================
# 4. One gene_exp.diff as input; GENIST and RTP-STAR tabs from Tux results; preset clusters
# ======================================================================================================================================================
d4 <- file.path(out, "V10_4_cuffdiff_single_file_modes"); dir.create(d4, recursive = TRUE, showWarnings = FALSE)
tx <- sim_tux(file.path(d4, "tux_results_project"), 1500, c("control", "heat"), c("0h", "1h", "3h", "6h"), reps = 3, responders = 0.2)
dir.create(file.path(d4, "tux_results_project", "00_design")); jsonlite::write_json(list(project_name = "tux toy", condition_names = c("control", "heat"),
  timepoint_names = c("0h", "1h", "3h", "6h"), n_conditions = 2, n_timepoints = 4), file.path(d4, "tux_results_project", "00_design", "design.json"), auto_unbox = TRUE)
dir.create(file.path(d4, "table_with_single_diff_file")); sd <- file.path(d4, "table_with_single_diff_file")
wcsv(tx$fpkm, file.path(sd, "fpkm_table.csv"))
sh <- tx$sheet; wcsv(sh, file.path(sd, "sample_sheet.csv"))
single <- tx$diff[tx$diff$sample_1 %in% paste0("control_", c("0h", "1h", "3h", "6h")) & sub("^heat_", "control_", tx$diff$sample_2) == tx$diff$sample_1, ]
utils::write.table(single, file.path(sd, "gene_exp.diff"), sep = "\t", quote = FALSE, row.names = FALSE)
set.seed(5); sig <- tx$diff[tx$diff$q_value < 0.05, ]; genes_de <- unique(sig$gene[order(sig$q_value)]); genes_de <- genes_de[seq_len(min(150, length(genes_de)))]      # the 150 most significant genes
tp <- file.path(d4, "rtpstar_preset_clusters"); dir.create(tp)
wcsv(data.frame(gene = genes_de), file.path(tp, "genes_file.csv"))
ex <- tx$fpkm[tx$fpkm$gene %in% genes_de, ]; wcsv(ex, file.path(tp, "expression_file.csv"))
wcsv(data.frame(gene = genes_de, name = genes_de, cluster = rep(1:3, length.out = length(genes_de))), file.path(tp, "my_clusters.csv"))
readme(d4, "V10_4: one Cuffdiff file as input, and the GENIST / RTP-STAR tabs started from Tux results", c(
  "WHAT THIS TESTS (all of these were gaps in the earlier test sets)",
  "  - a SINGLE gene_exp.diff file given under 'Differential-expression results' (not a zip, not a folder)",
  "  - the GENIST tab, 'Start from your Tux results'; the RTP-STAR tab, 'Start from your Tux results'; RTP-STAR with your own clusters file", "",
  "PART A: table_with_single_diff_file",
  "  - Experimental Design. Input: fpkm_table.csv (FPKM values), unit: 'Already normalized'. Differential-expression results: gene_exp.diff (the one file).",
  "  - Conditions 2 (control, heat), time points 4 (0h, 1h, 3h, 6h); sample sheet: sample_sheet.csv or blank.",
  "  - Tick Shannon, Modified Shannon, PCA, Trajectories, GENIST, interactive browser.",
  "  EXPECTED: the comparisons are matched to the conditions and time points by name; DESeq2 is greyed out (normalized values); entropy steps run on the file's q-values;",
  "  GENIST runs (no TF list: a note says it takes longer). The 'Compared' line of the confirmation window reads heat vs control.", "",
  "PART B: the GENIST tab from Tux results",
  "  - GENIST tab > Start from your Tux results > folder: tux_results_project. EXPECTED: 'Found the Tux results in .../diff_out'; the time points are listed in the",
  "    order 0h, 1h, 3h, 6h (read from 00_design/design.json); Build the genes file and the time course file: genes and time-course boxes fill in, up to 500 genes.",
  "  - Run GENIST (TF file empty is allowed: warning only). EXPECTED: results tab shows the network; output file name starts with the genes file name",
  "    (genist_genes_from_tux_GENIST_network_for_Cytoscape.txt) and a *_network.png sits beside it.",
  "  - Repeat with the folder tux_results_project/01_quantification/tux/diff_out: StRym still finds the design record in the parent folders, so the order is the same.",
  "    To see the other case copy diff_out ALONE to a new folder (outside tux_results_project) and choose that: a note says the time order is Tux's own and must be checked.", "",
  "PART C: the RTP-STAR tab from Tux results",
  "  - RTP-STAR tab > Start from your Tux results > tux_results_project. EXPECTED: 'The Tux run already made these files' is NOT shown (this fixture has no",
  "    rtpstar_inputs folder); press Build: genes file (differentially expressed genes) and expression file (24 sample columns) fill in.",
  "  - Run it (no TF file: warning only). EXPECTED: network table + picture, results tab shows the network.", "",
  "PART D: RTP-STAR with your own clusters",
  "  - RTP-STAR tab: Genes file: rtpstar_preset_clusters/genes_file.csv, Expression: expression_file.csv, Advanced > 'Use my own clusters file': my_clusters.csv",
  "    (gene, name, cluster number in column 3). EXPECTED: runs; clusters.csv in the output folder lists the clusters 1-3; hubs connect the clusters.", "",
  "DEBUGGING HINTS",
  "  - Cuffdiff file matching: strym_load_diff() and strym_resolve_diff_groups() in analysis/common.R. Tux-result readers: tux_find_results() in genist.R,",
  "    tux_rtpstar_inputs() in tuxedo.R."))
cat("4 written\n")

# ======================================================================================================================================================
# 5. GENIST / RTP-STAR tabs on their own files; names, pictures, optional TF list; extra result folders in the browser
# ======================================================================================================================================================
d5 <- file.path(out, "V10_5_network_tabs_browser"); dir.create(d5, recursive = TRUE, showWarnings = FALSE)
s5 <- sim_counts(1200, c("control", "stress"), c("t1", "t2", "t3", "t4", "t5"), reps_all(c("control", "stress"), c("t1", "t2", "t3", "t4", "t5"), 3), responders = 0.2, id_prefix = "AT3G")
wcsv(s5$counts, file.path(d5, "counts.csv")); wcsv(s5$sheet, file.path(d5, "sample_sheet.csv"))
M <- as.matrix(s5$counts[-1]); cpm <- sweep(M, 2, colSums(M) / 1e6, "/")
gen <- s5$counts$gene[s5$truth$kind != "flat"][1:70]
tcm <- sapply(c("t1", "t2", "t3", "t4", "t5"), function(t) rowMeans(cpm[, grepl(paste0("^stress_", t, "_"), colnames(cpm)), drop = FALSE]))
own <- file.path(d5, "own_files"); dir.create(own)
wcsv(data.frame(gene = gen), file.path(own, "genes_file.csv")); wcsv(data.frame(gene = s5$counts$gene, round(tcm, 3)), file.path(own, "time_course_file.csv"))
wcsv(data.frame(gene = s5$counts$gene, round(cpm[, grepl("^stress_", colnames(cpm))], 3), check.names = FALSE), file.path(own, "expression_file_for_rtpstar.csv"))
wcsv(data.frame(gene = gen, TF = rep(c(1, 0, 0, 0), length.out = length(gen)), round(cpm[match(gen, s5$counts$gene), grepl("^stress_", colnames(cpm))], 3), check.names = FALSE), file.path(own, "expression_file_with_TF_column.csv"))
wcsv(data.frame(gene = s5$counts$gene[sample(1200, 120)]), file.path(own, "TF_list.csv"))
utils::write.table(data.frame(gene = s5$counts$gene, symbol = sprintf("SYM%04d", seq_len(1200))), file.path(own, "gene_names.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)   # no header, like the TuxNet gene-name files
readme(d5, "V10_5: GENIST and RTP-STAR tabs on their own files; file names, pictures, optional TF list; extra result folders in the browser", c(
  "WHAT THIS TESTS",
  "  - network table names (not cytoscape_table.txt / biograph.txt any more), the picture next to each table, and the TF list being optional (warning only)",
  "  - results saved OUTSIDE a project reaching the project's interactive browser through 'More result folders for the interactive browser'", "",
  "PART A: GENIST tab on its own files (folder own_files)",
  "  - Genes file: genes_file.csv. Time course file: time_course_file.csv. Leave the TF box EMPTY (use the link 'No TF list' if the Arabidopsis list is filled in).",
  "    Output folder: a new folder outside any StRym project, e.g. <somewhere>/network_results.",
  "  - EXPECTED: when you pick the genes file the output file name becomes genes_file_GENIST_network_for_Cytoscape.txt; a warning says that without a TF list",
  "    the network takes longer; it does NOT stop you from pressing Run. After the run: Results tab shows the network, the folder holds the table, a *_network.png",
  "    and clusters.xlsx. Typing your own name in 'Output file name' before picking the genes file keeps your name.",
  "  - Run again with TF_list.csv in the TF box: the warning disappears.", "",
  "PART B: RTP-STAR tab on its own files",
  "  - Genes file: genes_file.csv. Expression: expression_file_with_TF_column.csv (it has a TF column: no TF warning). Output: a second new folder.",
  "  - EXPECTED: output name genes_file_RTPSTAR_network_for_Cytoscape.txt; after the run a network table, a picture, clusters.csv (if clustering used).",
  "  - Run once more with expression_file_for_rtpstar.csv (no TF column) and no TF file: advisory warning only; it runs.", "",
  "PART C: the interactive browser with extra folders",
  "  - Experimental Design: counts.csv, sample_sheet.csv (conditions control, stress; time points t1..t5), output: a third new folder (the project).",
  "    Tick DESeq2, Shannon, Data overview, PCA, interactive browser. In Analysis settings > 'More result folders for the interactive browser' enter the two output",
  "    folders of parts A and B (separate with a comma). Run.",
  "  - EXPECTED: strym_browser.html > Tables and Plots list the GENIST and RTP-STAR tables and pictures under 'extra_results/<folder name>/...' next to the project's own.",
  "    Without the setting they are not in the browser (they are outside the project).",
  "  - Put the output folders of A and B INSIDE the project's 06_regulatory_networks folder instead and rebuild the browser (Run only this): they appear without the setting.", "",
  "DEBUGGING HINTS",
  "  - Names: strym_network_filename() in io_utils.R; pictures: strym_network_png() in io_utils.R; browser: analysis/web_browser.R."))
cat("5 written\n")

# ---- a short index ------------------------------------------------------------------------------------------------------------------------------
wr(c("StRym version 9/10 test cases", "", "Folders 1 and 2 are meant to be copied to a USB stick or external drive (they test mounted drives).",
     "Each folder has a README.txt: what to do, what a correct run shows, and where to look when something is wrong.", "",
     "V10_1_drive_reads_tux_to_rtpstar   reads + GFF3 on a drive; Tux; RTP-STAR from the Tux results; unplug test",
     "V10_2_drive_timecourse_clutter     counts on a drive; output folder with other files; sample sheet follows the page; PCA levels; heatmaps; GENIST",
     "V10_3_heatmap_pca_edge_cases       tiny / 20 000-gene / many-condition / one-responder tables",
     "V10_4_cuffdiff_single_file_modes   single gene_exp.diff; GENIST and RTP-STAR tabs from Tux results; own clusters",
     "V10_5_network_tabs_browser         network file names and pictures; optional TF list; extra result folders in the browser", "",
     "Made by tests/make_test_sets.R (made-up data; truth files say which genes respond)."), file.path(out, "INDEX.txt"))
cat("done:", out, "\n")
