# A synthetic Tux experiment, written in two layouts:
#   * StRym's:  <root>/project/01_quantification/{sample_sheet.csv, counts.csv, tux/{samples_used.csv, diff_out/{genes.read_group_tracking, gene_exp.diff}}}
#               one Cuffdiff run over all groups (group = <condition>_<time point>)
#   * the original baseline scripts':  <root>/orig/genes.read_group_tracking_<genotype>_T<n>, gene_exp.diff_<genotype>_T<n>
#               one pair of files per genotype x time point (condition names contain "Heated" / "Control")
# Both hold the SAME numbers, so StRym's generalized steps can be compared with the original scripts.
# Genes: flat, sustained up / down after the baseline, a transient peak (genotype A only), control drift (varies without heat)
# and noisy genes.
synth_tux <- function(root, genotypes = c("A", "B"), n_tp = 5, n_rep = 3, n_genes = 400, seed = 7, orig = TRUE) {
  set.seed(seed)
  genes <- sprintf("XLOC_%06d", seq_len(n_genes)); gname <- sprintf("GENE%04d", seq_len(n_genes))
  base <- exp(rnorm(n_genes, log(30), 1.2))
  tps <- paste0("T", seq_len(n_tp) - 1)
  kind <- sample(c("flat", "up", "down", "transient", "drift", "noise"), n_genes, TRUE, prob = c(.3, .17, .15, .1, .13, .15))
  drift_prof <- t(sapply(seq_len(n_genes), function(i) { v <- rep(0, n_tp); v[sample(2:n_tp, 1)] <- sample(c(-1.6, 1.6), 1); v }))
  prof <- function(i, g, heated) {                                             # log2 effect per time point
    v <- rep(0, n_tp)
    if (kind[i] == "drift") v <- drift_prof[i, ]
    if (heated) {
      if (kind[i] == "up") v[3:n_tp] <- v[3:n_tp] + 1.9
      if (kind[i] == "down") v[3:n_tp] <- v[3:n_tp] - 1.9
      if (kind[i] == "transient" && g == genotypes[1]) v[3] <- v[3] + 2.2
    }
    v
  }
  conds <- c(rbind(paste0(genotypes, "_Control"), paste0(genotypes, "_Heated")))
  rows <- list(); fp <- list()
  for (g in genotypes) for (heated in c(FALSE, TRUE)) for (t in seq_len(n_tp)) for (r in seq_len(n_rep)) {
    cond <- paste0(g, if (heated) "_Heated" else "_Control")
    sd_ <- ifelse(kind == "noise", 0.7, 0.12)
    f <- base * 2^sapply(seq_len(n_genes), function(i) prof(i, g, heated)[t]) * exp(rnorm(n_genes, 0, sd_))
    rows[[length(rows) + 1]] <- data.frame(genotype = g, heated = heated, condition = cond, t = t, tp = tps[t], replicate = r - 1L, gene = genes, FPKM = f, stringsAsFactors = FALSE)
  }
  R <- do.call(rbind, rows)
  # q-values: t-test on log2(FPKM + 1) heated vs control per genotype and time point, Benjamini-Hochberg over genes
  Q <- list()
  for (g in genotypes) for (t in seq_len(n_tp)) {
    h <- R[R$genotype == g & R$t == t & R$heated, ]; c0 <- R[R$genotype == g & R$t == t & !R$heated, ]
    mh <- matrix(log2(h$FPKM + 1), nrow = n_genes); mc <- matrix(log2(c0$FPKM + 1), nrow = n_genes)
    p <- vapply(seq_len(n_genes), function(i) tryCatch(stats::t.test(mh[i, ], mc[i, ])$p.value, error = function(e) 1), 1)
    p[is.na(p)] <- 1
    Q[[paste(g, t)]] <- stats::p.adjust(p, "BH")
  }
  diffrow <- function(g, t, s1, s2) data.frame(test_id = genes, gene_id = genes, gene = gname, locus = "chr1:1-2", sample_1 = s1, sample_2 = s2, status = "OK",
                                               value_1 = 1, value_2 = 1, log2.fold_change. = 0, test_stat = 0, p_value = 0.01, q_value = Q[[paste(g, t)]], significant = "no")
  tracking <- function(x, cond_name) data.frame(tracking_id = x$gene, condition = cond_name, replicate = x$replicate, raw_frags = round(x$FPKM * 10),
                                                internal_scaled_frags = round(x$FPKM * 10), external_scaled_frags = round(x$FPKM * 10), FPKM = x$FPKM,
                                                effective_length = 1000, status = "OK")
  wr <- function(df, f) utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- StRym layout ----
  pd <- file.path(root, "project"); q1 <- file.path(pd, "01_quantification"); td <- file.path(q1, "tux", "diff_out"); dir.create(td, recursive = TRUE)
  R$group <- paste(R$condition, R$tp, sep = "_"); R$sample <- sprintf("%s_r%d", R$group, R$replicate + 1L)
  used <- unique(R[, c("group", "replicate", "sample")]); used <- data.frame(group = used$group, k = used$replicate + 1L, file = paste0(used$sample, ".fastq.gz"), mate = NA, sample = used$sample)
  utils::write.csv(used, file.path(q1, "tux", "samples_used.csv"), row.names = FALSE)
  sheet <- unique(R[, c("sample", "condition", "tp", "replicate")]); names(sheet) <- c("sample", "condition", "time_point", "replicate"); sheet$replicate <- sheet$replicate + 1L
  sheet$file <- paste0(sheet$sample, ".fastq.gz"); utils::write.csv(sheet, file.path(q1, "sample_sheet.csv"), row.names = FALSE)
  utils::write.csv(data.frame(gene = gname, x = 1), file.path(q1, "counts.csv"), row.names = FALSE)         # only has to exist
  wr(do.call(rbind, lapply(split(R, R$group), function(x) tracking(x, x$group[1]))), file.path(td, "genes.read_group_tracking"))
  dl <- list()
  for (g in genotypes) for (t in seq_len(n_tp)) dl[[length(dl) + 1]] <- diffrow(g, t, paste0(g, "_Control_", tps[t]), paste0(g, "_Heated_", tps[t]))
  extra <- diffrow(genotypes[1], 1, paste0(genotypes[1], "_Control_", tps[1]), paste0(genotypes[1], "_Control_", tps[2])); extra$q_value <- runif(n_genes)   # an unrelated comparison
  wr(do.call(rbind, c(dl, list(extra))), file.path(td, "gene_exp.diff"))

  # ---- original layout ----
  od <- NULL
  if (orig) {
    od <- file.path(root, "orig"); dir.create(od)
    for (g in genotypes) for (t in seq_len(n_tp)) {
      x <- R[R$genotype == g & R$t == t, ]
      wr(rbind(tracking(x[x$heated, ], paste0(g, "_Heated")), tracking(x[!x$heated, ], paste0(g, "_Control"))), file.path(od, sprintf("genes.read_group_tracking_%s_T%d", g, t)))
      wr(diffrow(g, t, paste0(g, "_Control"), paste0(g, "_Heated")), file.path(od, sprintf("gene_exp.diff_%s_T%d", g, t)))
    }
  }
  list(project = pd, orig = od, genes = genes, gene_names = gname, kind = setNames(kind, genes), conds = conds, tps = tps, genotypes = genotypes, data = R)
}

# A design for the synthetic project (raw reads, so the steps read the Tux results).
synth_design <- function(syn, ...) {
  d <- strym_default_design()
  d$project_name <- "synth"; d$project_dir <- syn$project; d$data_type <- "bulk"; d$data_kind <- "fastq"; d$data_path <- file.path(syn$project, "reads")
  d$n_conditions <- length(syn$conds); d$condition_names <- syn$conds; d$n_timepoints <- length(syn$tps); d$timepoint_names <- syn$tps
  a <- list(...); for (k in names(a)) d[[k]] <- a[[k]]
  d
}
