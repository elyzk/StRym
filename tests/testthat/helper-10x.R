# Synthetic single-cell data in the 10x layout (<root>/<sample>/raw_matrix/{matrix.mtx, features.tsv, barcodes.tsv}), with known truth:
#   3 cell types, each with 30 marker genes (Gene1-30, Gene31-60, Gene61-90); the barcode starts with the true type (ct1-, ct2-, ct3-);
#   in "Treated" samples genes Gene700-Gene720 are 5x higher in cell type 2 only (a condition effect inside one cluster).
synth_10x <- function(root, samples = c("Control-1", "Control-2", "Treated-1", "Treated-2"), cells = 200, n_genes = 800, seed = 11, gz = c(TRUE, FALSE)) {
  set.seed(seed); gz <- rep(gz, length.out = length(samples)); genes <- paste0("Gene", seq_len(n_genes)); dirs <- character()
  for (si in seq_along(samples)) {
    s <- samples[si]; type <- sample(1:3, cells, TRUE, prob = c(.4, .35, .25)); treated <- grepl("^Treated", s)
    lam <- matrix(0.5, n_genes, cells)
    for (k in 1:3) lam[((k - 1) * 30 + 1):(k * 30), type == k] <- 7
    eg <- if (n_genes >= 720) 700:720 else (n_genes - 20):n_genes                      # the condition-effect genes
    if (treated) lam[eg, type == 2] <- lam[eg, type == 2] * 5
    lam <- sweep(lam, 2, runif(cells, 0.7, 1.3), "*")                                  # depth differs between cells
    m <- Matrix::Matrix(matrix(stats::rpois(length(lam), lam), n_genes, cells), sparse = TRUE)
    d <- file.path(root, s, "raw_matrix"); dir.create(d, recursive = TRUE); dirs <- c(dirs, d)
    barcodes <- sprintf("ct%d-%s%04d", type, "AAAC", seq_len(cells))                   # same barcode suffixes in every sample, as in real runs
    ext <- if (gz[si]) ".gz" else ""
    w <- function(lines, f) { con <- if (gz[si]) gzfile(file.path(d, paste0(f, ext)), "wt") else file(file.path(d, f), "wt"); writeLines(lines, con); close(con) }
    tmp <- tempfile(fileext = ".mtx"); Matrix::writeMM(m, tmp); w(readLines(tmp), "matrix.mtx")
    w(sprintf("%s\t%s\tGene Expression", genes, genes), "features.tsv"); w(barcodes, "barcodes.tsv")
  }
  list(root = root, dirs = dirs, samples = samples)
}
