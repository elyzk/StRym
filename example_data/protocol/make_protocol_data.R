# Regenerates the protocol example experiments (small, synthetic, Arabidopsis gene ids, known truth).
#   Rscript make_protocol_data.R [output folder]            (default: this folder)
#
#   exp1_heatshock_counts/      raw counts, control vs heat, 4 time points x 3 replicates (24 samples), 1,500 genes (120 TFs); 150 genes respond to heat.
#                               Use for: the whole Experimental Design page, DESeq2, entropy steps (q-values calculated by StRym), PCA, trajectories, GENIST.
#   exp2_fpkm_plus_cuffdiff/    the same biology as a merged FPKM table (one column per condition x time point, no replicates) plus Cuffdiff comparison
#                               files, one per time point, in folders inside one .zip, with differently written sample labels.
#                               Use for: "FPKM table + Cuffdiff results from another run".
#   exp5_messy_data/            exp1 with the problems real data has: a missing time point, empty cells, "-" and text in the table, a duplicated gene,
#                               a nearly empty library, an all-zero sample, a listed sample that is not in the table. Use for: the warnings and checks.
here <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)))
args <- commandArgs(TRUE); out <- if (length(args)) args[1] else here
app <- normalizePath(file.path(here, "..", ".."))
suppressPackageStartupMessages({ library(readxl); library(DESeq2) })
set.seed(2026)

gn <- as.data.frame(read_excel(file.path(app, "reference", "GeneNames.xlsx"), col_names = FALSE, .name_repair = "minimal"))
ids <- as.character(gn[[1]]); sym <- as.character(gn[[2]]); names(sym) <- ids
tfs <- intersect(as.character(read_excel(file.path(app, "reference", "Arabidopsis_TFs.xlsx"), col_names = FALSE, .name_repair = "minimal")[[1]]), ids)
genes <- c(sample(tfs, 120), sample(setdiff(ids, tfs), 1380)); genes <- sample(genes)
ng <- length(genes); tps <- c("0h", "1h", "3h", "6h"); reps <- 1:3

# truth: 150 heat-responsive genes (40 of them TFs), four response shapes (log2 fold change at each time point)
resp <- c(sample(intersect(genes, tfs), 40), sample(setdiff(genes, tfs), 110))
shape <- list(early = c(0, 2.0, 1.4, 0.6), late = c(0, 0.3, 1.3, 2.2), down = c(0, -1.1, -1.7, -1.7), transient = c(0, 2.4, 0.8, 0))
lfc <- matrix(0, ng, length(tps), dimnames = list(genes, tps))
kind <- sample(names(shape), length(resp), TRUE); for (i in seq_along(resp)) lfc[resp[i], ] <- shape[[kind[i]]] * runif(1, 0.8, 1.2)
base <- pmin(pmax(exp(rnorm(ng, log(150), 1.1)), 5), 20000); names(base) <- genes

grid <- expand.grid(rep = reps, tp = tps, cond = c("control", "heat"), stringsAsFactors = FALSE)
grid$sample <- sprintf("%s_%s_rep%d", grid$cond, grid$tp, grid$rep)
depth <- runif(nrow(grid), 0.75, 1.3)
mu <- sapply(seq_len(nrow(grid)), function(j) base * depth[j] * 2^(if (grid$cond[j] == "heat") lfc[, grid$tp[j]] else 0))
counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 12), ng, nrow(grid), dimnames = list(genes, grid$sample))
sheet <- data.frame(sample = grid$sample, condition = grid$cond, time_point = grid$tp, replicate = grid$rep)
w <- function(x, f) { dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); utils::write.csv(x, f, row.names = FALSE, quote = TRUE) }

# ---- exp1: counts ----
d1 <- file.path(out, "exp1_heatshock_counts")
w(data.frame(gene_id = genes, gene_name = unname(ifelse(is.na(sym[genes]), "", sym[genes])), counts, check.names = FALSE), file.path(d1, "heatshock_counts.csv"))
w(sheet, file.path(d1, "heatshock_sample_sheet.csv"))
w(data.frame(gene_id = resp, kind = kind, is_TF = resp %in% tfs), file.path(d1, "TRUTH_heat_responsive_genes.csv"))
writeLines(c("Synthetic heat-shock time course (Arabidopsis gene ids).",
             "heatshock_counts.csv       raw counts, 1,500 genes x 24 samples (gene_name is an annotation column, not a sample)",
             "heatshock_sample_sheet.csv sample, condition (control / heat), time_point (0h 1h 3h 6h), replicate",
             "TRUTH_heat_responsive_genes.csv  the 150 genes that really respond (for checking the results)"), file.path(d1, "README.txt"))

# ---- exp2: FPKM table + Cuffdiff files ----
d2 <- file.path(out, "exp2_fpkm_plus_cuffdiff")
lib <- colSums(counts); cpm <- sweep(counts, 2, lib, "/") * 1e6
groups <- paste(sheet$condition, sheet$time_point, sep = "_"); gl <- unique(groups)
gm <- sapply(gl, function(g) rowMeans(cpm[, groups == g, drop = FALSE])) / 10        # "FPKM-like" scale
xloc <- sprintf("XLOC_%06d", seq_len(ng))
w(data.frame(tracking_id = xloc, gene_short_name = genes, round(gm, 5), check.names = FALSE), file.path(d2, "merged_FPKM.csv"))
w(data.frame(sample = gl, condition = sub("_.*", "", gl), time_point = sub(".*_", "", gl), replicate = 1), file.path(d2, "fpkm_sample_sheet.csv"))
tmp <- file.path(d2, "cuffdiff_results"); dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
for (t in tps) {
  s <- sheet[sheet$time_point == t, ]; cnt <- counts[, s$sample]
  dds <- DESeqDataSetFromMatrix(cnt, data.frame(cond = factor(s$condition, levels = c("control", "heat")), row.names = s$sample), ~ cond)
  r <- results(DESeq(dds, quiet = TRUE), contrast = c("cond", "heat", "control"))
  q <- r$padj; q[is.na(q)] <- 1
  de <- data.frame(test_id = xloc, gene_id = xloc, gene = genes, locus = "chr:1-100", sample_1 = paste0(t, "_Control"), sample_2 = paste0(t, "_Heat"),
                   status = "OK", value_1 = round(gm[, paste0("control_", t)], 5), value_2 = round(gm[, paste0("heat_", t)], 5),
                   `log2(fold_change)` = r$log2FoldChange, test_stat = r$stat, p_value = r$pvalue, q_value = q, significant = ifelse(q < 0.05, "yes", "no"), check.names = FALSE)
  f <- file.path(tmp, paste0(t, "_diff_out"), "gene_exp.diff"); dir.create(dirname(f), showWarnings = FALSE); utils::write.table(de, f, sep = "\t", quote = FALSE, row.names = FALSE)
}
# a small zip writer in base R (files stored without compression), so no zip program is needed
write_zip <- function(zipfile, root) {
  files <- list.files(root, recursive = TRUE, all.files = TRUE); tab <- integer(256)
  for (n in 0:255) { c <- n; for (k in 1:8) c <- if (bitwAnd(c, 1L) == 1L) bitwXor(-306674912L, bitwShiftR(c, 1L)) else bitwShiftR(c, 1L); tab[n + 1] <- c }
  crc32 <- function(b) { c <- -1L; for (x in as.integer(b)) c <- bitwXor(tab[bitwAnd(bitwXor(c, x), 255L) + 1L], bitwShiftR(c, 8L)); bitwNot(c) }
  con <- file(zipfile, "wb"); on.exit(close(con)); cd <- list(); pos <- 0
  u2 <- function(x) writeBin(as.integer(x), con, size = 2, endian = "little"); u4 <- function(x) writeBin(as.integer(x), con, size = 4, endian = "little")
  for (f in files) {
    b <- readBin(file.path(root, f), "raw", file.size(file.path(root, f))); nm <- charToRaw(f); crc <- crc32(b)
    cd[[length(cd) + 1]] <- list(nm = nm, crc = crc, n = length(b), pos = pos)
    u4(0x04034b50); u2(20); u2(0); u2(0); u2(0); u2(33); u4(crc); u4(length(b)); u4(length(b)); u2(length(nm)); u2(0)   # stored; date 1980-01-01
    writeBin(nm, con); writeBin(b, con); pos <- pos + 30 + length(nm) + length(b)
  }
  start <- pos
  for (e in cd) { u4(0x02014b50); u2(20); u2(20); u2(0); u2(0); u2(0); u2(33); u4(e$crc); u4(e$n); u4(e$n); u2(length(e$nm)); u2(0); u2(0); u2(0); u2(0); u4(0); u4(e$pos)
                  writeBin(e$nm, con); pos <- pos + 46 + length(e$nm) }
  u4(0x06054b50); u2(0); u2(0); u2(length(cd)); u2(length(cd)); u4(pos - start); u4(start); u2(0)
}
write_zip(file.path(d2, "cuffdiff_results.zip"), tmp)
if (file.exists(file.path(d2, "cuffdiff_results.zip"))) unlink(tmp, recursive = TRUE)
writeLines(c("The same experiment as a table of FPKM-like values (one column per condition x time point) plus Cuffdiff-style comparison files.",
             "merged_FPKM.csv           tracking_id, gene_short_name, then one column per group (control_0h ... heat_6h)",
             "cuffdiff_results.zip      one gene_exp.diff per time point, each in its own folder inside the zip (labels written like 0h_Control / 0h_Heat)",
             "fpkm_sample_sheet.csv     sample = the table's column names"), file.path(d2, "README.txt"))

# ---- exp5: messy data ----
d5 <- file.path(out, "exp5_messy_data")
keep <- !(sheet$condition == "heat" & sheet$time_point == "3h")                       # a whole time point missing for one condition
cm <- counts[, sheet$sample[keep]]; sh <- sheet[keep, ]
cm[, "control_6h_rep3"] <- rbinom(ng, cm[, "control_6h_rep3"], 0.01)                  # a nearly empty library
cm[, "heat_1h_rep3"] <- 0                                                             # an all-zero sample
tab <- data.frame(gene_id = genes, notes = "", cm, check.names = FALSE, stringsAsFactors = FALSE)
tab[sample(ng, 25), sample(3:ncol(tab), 1)] <- NA
for (k in 1:12) tab[sample(ng, 1), sample(3:ncol(tab), 1)] <- "-"
tab <- rbind(tab, tab[5, ]); tab$gene_id[10] <- ""                                    # a duplicated gene and a gene without an id
w(tab, file.path(d5, "messy_counts.csv"))
w(rbind(sh, data.frame(sample = "heat_6h_rep4", condition = "heat", time_point = "6h", replicate = 4)), file.path(d5, "messy_sample_sheet.csv"))
writeLines(c("exp1 with the problems of real data (use to see how StRym warns and copes):",
             " * heat has no 3h samples at all   * empty cells and '-' in the counts   * a text column ('notes')",
             " * a gene listed twice and a row without an id   * control_6h_rep3 has ~1% of the reads   * heat_1h_rep3 is all zeros",
             " * heat_6h_rep4 is in the sample sheet but not in the table"), file.path(d5, "README.txt"))
cat("Written to", out, "\n"); print(list.files(out, recursive = TRUE))
