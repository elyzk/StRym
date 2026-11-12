# a small synthetic counts experiment (control vs heat, replicates, optional missing group), shared by several test files
mk_counts <- function(tps = c("0h", "1h", "3h"), reps = 3, ng = 120, drop = NULL, seed = 5) {
  set.seed(seed); dir <- tempfile("derived_"); dir.create(dir)
  g <- expand.grid(rep = seq_len(reps), tp = tps, cond = c("control", "heat"), stringsAsFactors = FALSE)
  g$sample <- sprintf("%s_%s_r%d", g$cond, g$tp, g$rep)
  de <- seq_len(20); base <- rep(100, ng)
  mu <- sapply(seq_len(nrow(g)), function(j) base * (if (g$cond[j] == "heat" && g$tp[j] != tps[1]) c(rep(6, 10), rep(1 / 6, 10), rep(1, ng - 20)) else 1))
  m <- matrix(stats::rnbinom(length(mu), mu = mu, size = 20), ng, nrow(g), dimnames = list(sprintf("G%03d", seq_len(ng)), g$sample))
  if (length(drop)) { keep <- !(g$cond == drop[1] & g$tp == drop[2]); m <- m[, keep]; g <- g[keep, ] }
  write.csv(data.frame(gene = rownames(m), m, check.names = FALSE), file.path(dir, "counts.csv"), row.names = FALSE)
  sh <- data.frame(sample = g$sample, condition = g$cond, time_point = g$tp, replicate = g$rep); write.csv(sh, file.path(dir, "sheet.csv"), row.names = FALSE)
  list(dir = dir, table = file.path(dir, "counts.csv"), sheet = file.path(dir, "sheet.csv"), m = m, sh = sh)
}

