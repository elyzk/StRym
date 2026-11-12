# ---- BDe score: compare with a literal transcription of the MATLAB loops (uses gamma(), not lgamma()) ----
matlab_bde <- function(E, t, L, ESS = 1e-15) {      # E: time x m (m = 1..3), values 0..L-1
  m <- ncol(E); q <- L^m; bde <- 0
  cfgs <- as.matrix(expand.grid(rep(list(0:(L - 1)), m)))
  for (r in seq_len(nrow(cfgs))) {
    Nj <- apply(E, 1, function(row) all(row == cfgs[r, ]))
    bde <- bde + log(gamma(ESS / q) / gamma(sum(Nj) + ESS / q))
    for (d in 0:(L - 1)) {
      Njk <- Nj & t == d
      bde <- bde + log(gamma(sum(Njk) + ESS / (L * q)) / gamma(ESS / (L * q)))
    }
  }
  bde
}

test_that("bde_score equals the MATLAB formula for 1, 2 and 3 regulators and 2/3 levels", {
  set.seed(1)
  for (L in 2:3) for (m in 1:3) {
    Tn <- 9
    E <- matrix(sample(0:(L - 1), Tn * m, TRUE), Tn, m); t <- sample(0:(L - 1), Tn, TRUE)
    expect_equal(bde_score(t(E), t, L), matlab_bde(E, t, L), tolerance = 1e-8)
  }
})

test_that("a perfectly predictive regulator scores higher than an uninformative one", {
  t <- c(0, 0, 1, 1, 0, 1, 1, 0)
  expect_gt(bde_score(matrix(t, 1), t, 2), bde_score(matrix(c(0, 1, 0, 1, 0, 1, 0, 1), 1), t, 2))
})

# ---- preprocessing on a hand-worked example ----
# A changes at t=2,3 ; B changes at t=3,4 (one step after A) ; C is flat
test_that("candidate regulators follow the MATLAB definition (lag-1 and same-time)", {
  tc <- rbind(A = c(1, 3, 9, 9), B = c(1, 1, 3, 9), C = c(5, 5, 5, 5))
  pre <- genist_preprocess(tc, fc = 1.3, pct = 0.3, n_levels = 2)
  expect_equal(unname(which(pre$reg1, arr.ind = TRUE)), matrix(c(1, 2), 1))   # only A -> B (regulator row 1, target col 2)
  expect_true(pre$reg2[1, 2]); expect_true(pre$reg2[2, 1])                  # same-time: both directions
  expect_false(any(pre$reg1[3, ])); expect_false(any(pre$reg2[3, ]))         # flat gene regulates nothing
  expect_equal(pre$set_max_reg, 1)
  expect_equal(unname(pre$disc["A", ]), c(0, 0, 1, 1))                       # >= row mean
})

test_that("TF restriction removes non-TF regulators", {
  tc <- rbind(A = c(1, 3, 9, 9), B = c(1, 1, 3, 9), C = c(5, 5, 5, 5))
  pre <- genist_preprocess(tc, 1.3, 0.3, 2, tf_pos = c(FALSE, TRUE, TRUE))
  expect_false(any(pre$reg1[1, ]))                  # A is not a TF
})

test_that("3-level discretization uses MATLAB-style quantiles", {
  tc <- rbind(G = c(1, 2, 3, 4, 5, 6, 7, 8, 9))
  pre <- genist_preprocess(rbind(tc, tc + 1), 1.3, 0.3, 3)
  expect_equal(unname(pre$disc[1, ]), c(0, 0, 0, 1, 1, 1, 2, 2, 2))
})

# ---- inference: known activation and repression with a one-step delay ----
test_that("activator and repressor are recovered with the right sign (time lapse 1)", {
  A <- c(1, 1, 4, 4, 16, 16, 64, 64)                 # up at t=3,5,7
  B <- c(1, 1, 1, 4, 4, 16, 16, 64)                  # follows A one step later
  D <- c(256, 256, 256, 64, 64, 16, 16, 4)           # falls one step after A rises
  tc <- rbind(A = A, B = B, D = D)
  pre <- genist_preprocess(tc, 1.3, 0.3, 2)
  res <- genist_infer(pre, time_lapse = 1, low_conn = 0)               # keep every scored edge
  expect_gt(res$network[1, 2], 0); expect_equal(res$sign[1, 2], 1)     # A -> B activates
  expect_gt(res$network[1, 3], 0); expect_equal(res$sign[1, 3], -1)    # A -| D represses
  expect_equal(res$network[2, 1], 1)                                   # B -> A is perfectly predictive: top weight (1)
  # the default low-connection cut-off then drops the much weaker A -> B / A -> D edges
  res2 <- genist_infer(pre, time_lapse = 1, low_conn = 0.2)
  expect_equal(res2$network[1, 2], 0); expect_equal(res2$network[2, 1], 1)
})

test_that("network weights are relative (max 1)", {
  set.seed(3)
  tc <- matrix(exp(rnorm(6 * 8)), 6, 8, dimnames = list(letters[1:6], NULL))
  pre <- genist_preprocess(tc, 1.3, 0.3, 2)
  res <- genist_infer(pre, 1, 0.2)
  expect_lte(max(abs(res$network)), 1 + 1e-12)
})

test_that("edge table lists isolated nodes with width 0 and sign 1", {
  net <- matrix(0, 3, 3); net[1, 2] <- 0.5; sg <- matrix(0, 3, 3); sg[1, 2] <- -1
  tab <- genist_edge_table(c("x", "y", "z"), net, sg)
  expect_equal(tab$source, c("x", "z")); expect_equal(tab$target, c("y", "z"))
  expect_equal(tab$width, c(0.5, 0)); expect_equal(tab$sign, c(-1, 1))
})

# ---- end to end on files ----
test_that("genist_run writes the Cytoscape table and cluster lists", {
  skip_if_not_installed("writexl")
  set.seed(5); d <- tempfile("g"); dir.create(d)
  genes <- sprintf("G%02d", 1:24)
  tc <- matrix(exp(rnorm(24 * 7, 3, .6)), 24, 7); tc[2, ] <- tc[1, ] * c(1, 1, 3, 3, 9, 9, 27)
  tc[3, ] <- c(tc[1, 1], tc[2, 1:6])
  write.csv(data.frame(gene = genes, tc), file.path(d, "tc.csv"), row.names = FALSE)
  write.csv(data.frame(gene = genes), file.path(d, "genes.csv"), row.names = FALSE)
  cl <- matrix(exp(rnorm(24 * 4)), 24, 4); cl[1:12, 1:2] <- cl[1:12, 1:2] + 20
  write.csv(data.frame(gene = genes, cl), file.path(d, "clu.csv"), row.names = FALSE)
  quiet(genist_run(file.path(d, "genes.csv"), file.path(d, "tc.csv"), file.path(d, "clu.csv"),
                   output_folder = d, output_file = "net.txt", time_lapse = "both"))
  out <- read.delim(file.path(d, "net.txt"))
  expect_named(out, c("source", "target", "width", "sign"))
  expect_true(all(out$sign %in% c(-1, 0, 1)))
  expect_true(file.exists(file.path(d, "clusters.xlsx")))
})

test_that("genist_run reports problems clearly", {
  expect_error(genist_run("nope.csv", "nope2.csv", output_folder = tempdir()), "cannot be found")
  expect_error(genist_run("a", "b", output_folder = tempdir(), time_lapse = 2), "time_lapse")
  expect_error(genist_run("a", "b", output_folder = tempdir(), fc_thresh = 0.5), "fc_thresh")
  expect_error(genist_run("a", "b", output_folder = tempdir(), n_levels = 1), "n_levels")
})
