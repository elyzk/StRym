test_that("GENIE3 weights rank the true regulator first", {
  set.seed(2); n <- 30
  x1 <- rnorm(n); x2 <- rnorm(n); y <- 2 * x1 + rnorm(n, sd = .3)
  expr <- rbind(x1 = x1, x2 = x2, y = y)
  W <- genie3_weights(expr, tf_idx = 1:2, ntrees = 500, cores = 1, seed = 1)
  expect_gt(W["x1", "y"], W["x2", "y"]); expect_equal(W["y", "x1"], 0)       # only TFs regulate
  expect_equal(unname(diag(W)), c(0, 0, 0))
})

test_that("our GENIE3 agrees with the published GENIE3 package (rank correlation)", {
  skip_if_not_installed("GENIE3")
  set.seed(4); n <- 40; g <- 12
  expr <- matrix(rnorm(g * n), g, n, dimnames = list(paste0("g", 1:g), NULL))
  expr[8, ] <- expr[1, ] + expr[2, ] * .5 + rnorm(n, sd = .2); expr[9, ] <- expr[3, ] - expr[1, ] + rnorm(n, sd = .2)
  mine <- genie3_weights(expr, 1:6, ntrees = 2000, cores = 1, seed = 1)
  ref <- suppressWarnings(GENIE3::GENIE3(expr, regulators = paste0("g", 1:6), nTrees = 2000, nCores = 1, verbose = FALSE))
  ref_full <- matrix(0, g, g, dimnames = dimnames(mine)); ref_full[rownames(ref), colnames(ref)] <- ref
  keep <- mine != 0 | ref_full != 0
  expect_gt(cor(mine[keep], ref_full[keep], method = "spearman"), 0.8)
})

test_that("edge direction follows the time course", {
  set.seed(6); n <- 20; a <- rnorm(n)
  expr <- rbind(A = a, B = 2 * a + rnorm(n, sd = .1), C = -a + rnorm(n, sd = .1))
  tc <- rbind(A = 2^(0:4), B = 2^(0:4), C = 2^(4:0))
  r <- rtp_infer(expr, c(TRUE, FALSE, FALSE), tc, c("A", "B", "C"), ntrees = 300, seed = 1)
  e <- setNames(r$edges$dir, paste(r$edges$from, r$edges$to))
  expect_equal(unname(e["A B"]), "activates"); expect_equal(unname(e["A C"]), "inhibits")
  expect_equal(r$hub, "A")
  r2 <- rtp_infer(expr, c(TRUE, FALSE, FALSE), NULL, c("A", "B", "C"), ntrees = 300, seed = 1)
  expect_true(all(r2$edges$dir == "regulates"))
})

test_that("clusters are found for both clustering types", {
  # 60 genes: the number of clusters tried is floor(n/10-5)..floor(n/10+5) = 1..11 (the original's rule)
  set.seed(8)
  centers <- rbind(c(10, 0, 0, 0), c(0, 10, 0, 0), c(0, 0, 10, 10))
  d <- do.call(rbind, lapply(1:3, function(g) matrix(rep(centers[g, ], each = 20), 20, 4) + matrix(rnorm(80, sd = .3), 20, 4)))
  rownames(d) <- sprintf("g%02d", 1:60)
  for (ty in c("S", "T")) {
    cl <- rtp_cluster(d, ty, seed = 1)
    expect_length(unique(cl), 3)
    for (g in 1:3) expect_equal(length(unique(cl[(g - 1) * 20 + 1:20])), 1)
  }
})

test_that("gene labels use symbols but stay unique", {
  sym <- data.frame(id = c("a", "b", "c"), s = c("X", "X", "Z"))
  expect_equal(rtp_labels(c("a", "b", "c", "d"), sym), c("X", "b", "Z", "d"))   # second X falls back to its id
})

test_that("rtpstar_run: one iteration, then several iterations with an edge-proportion cut-off", {
  set.seed(9); d <- tempfile("r"); dir.create(d)
  genes <- sprintf("G%02d", 1:20); tfs <- genes[1:8]
  expr <- matrix(rnorm(20 * 9), 20, 9); expr[10, ] <- expr[1, ] * 2 + rnorm(9, sd = .1)
  write.csv(data.frame(gene = genes, expr), file.path(d, "expr.csv"), row.names = FALSE)
  write.csv(data.frame(gene = genes), file.path(d, "genes.csv"), row.names = FALSE)
  write.csv(data.frame(TF = tfs), file.path(d, "tfs.csv"), row.names = FALSE)
  quiet(rtpstar_run(file.path(d, "genes.csv"), file.path(d, "expr.csv"), d, "one.txt", tf_file = file.path(d, "tfs.csv"),
                    ntrees = 100, rf_seed = 1))
  one <- readLines(file.path(d, "one.txt")); expect_gt(length(one), 0)
  expect_true(all(grepl("^G[0-9]+ (activates|inhibits|regulates) G[0-9]+$", one)))
  expect_true(all(sub(" .*", "", one) %in% tfs))                              # only TFs regulate
  quiet(rtpstar_run(file.path(d, "genes.csv"), file.path(d, "expr.csv"), d, "multi.txt", tf_file = file.path(d, "tfs.csv"),
                    numiters = 3, maxprop = 0.5, ntrees = 100))
  multi <- readLines(file.path(d, "multi.txt"))
  expect_true(all(grepl(" [0-9.]+$", multi)))
  expect_true(file.exists(file.path(d, "biograph3.txt")))
})

test_that("rtpstar_run reports problems clearly", {
  expect_error(rtpstar_run("x.csv", "y.csv", tempdir()), "cannot be found")
})
