#!/usr/bin/env Rscript
# Generate a tiny synthetic LONG-READ (Oxford Nanopore-like) RNA-seq dataset to test StRym's long-read path (base R only).
#
# Uses the genome.fa / genes.gtf / truth.tsv written by make_synthetic.R (so the genes and the truth are the same), and writes
#   <out>/inputs/{CTL,TRT}{1..reps}_ont.fastq.gz   full-length-ish cDNA reads (0.4-2 kb), ~7% errors (substitutions and indels), Phred ~10-15,
#                                                  either strand (unstranded cDNA, as in Plasmidsaurus-type libraries)
# Genes labelled DE in truth.tsv are 6-8x higher/lower in TRT. Deterministic given --seed.
#
# usage: Rscript make_longread_data.R --ref DIR --out DIR [--reps 3] [--reads 120] [--seed 3]

BASES <- c("A", "C", "G", "T")
revcomp <- function(s) paste(rev(strsplit(chartr("ACGTN", "TGCAN", s), "", fixed = TRUE)[[1]]), collapse = "")

read_fasta <- function(path) {
  x <- trimws(readLines(path)); x <- x[nzchar(x)]; hdr <- startsWith(x, ">")
  grp <- cumsum(hdr)
  setNames(toupper(vapply(split(x[!hdr], grp[!hdr]), paste, "", collapse = "")), sub("\\s.*$", "", substring(x[hdr], 2)))
}

read_gtf <- function(path) {
  x <- readLines(path); x <- x[nzchar(x) & !startsWith(x, "#")]
  f <- strsplit(x, "\t", fixed = TRUE); f <- f[vapply(f, function(v) v[3] == "exon", NA)]
  gid <- vapply(f, function(v) sub('.*gene_id "([^"]*)".*', "\\1", v[9]), "")
  out <- list()
  for (i in seq_along(f)) {
    v <- f[[i]]; g <- out[[gid[i]]]
    if (is.null(g)) g <- list(chrom = v[1], strand = v[7], s = integer(), e = integer())
    g$s <- c(g$s, as.integer(v[4])); g$e <- c(g$e, as.integer(v[5])); out[[gid[i]]] <- g
  }
  lapply(out, function(g) { o <- order(g$s); g$s <- g$s[o]; g$e <- g$e[o]; g })
}

read_truth <- function(path) {
  if (!file.exists(path)) return(list())
  t <- utils::read.delim(path, stringsAsFactors = FALSE)             # gene_id, chrom, de, fold_change, base_reads
  setNames(lapply(seq_len(nrow(t)), function(i) list(de = t$de[i], fc = t$fold_change[i], base = t$base_reads[i])), t$gene_id)
}

mutate <- function(s, rate) {                                       # substitutions 45%, deletions 25%, insertions 30% of the errors
  ch <- strsplit(s, "", fixed = TRUE)[[1]]; r <- runif(length(ch)); rnd <- sample(BASES, length(ch), replace = TRUE)
  sub <- r < rate * 0.45; del <- r >= rate * 0.45 & r < rate * 0.70; ins <- r >= rate * 0.70 & r < rate
  ch[sub] <- rnd[sub]; ch[ins] <- paste0(ch[ins], sample(BASES, sum(ins), replace = TRUE))
  paste(ch[!del], collapse = "")
}

main <- function(argv = commandArgs(TRUE)) {
  o <- list(ref = NA_character_, out = NA_character_, reps = 3, reads = 120, seed = 3)
  for (i in seq(1, length(argv), by = 2)) {
    k <- sub("^--", "", argv[i]); if (!k %in% names(o)) stop("unknown option: ", argv[i])
    o[[k]] <- if (is.numeric(o[[k]])) as.numeric(argv[i + 1]) else argv[i + 1]
  }
  if (is.na(o$ref) || is.na(o$out)) stop("--ref and --out are required")
  set.seed(o$seed)
  genome <- read_fasta(file.path(o$ref, "genome.fa")); genes <- read_gtf(file.path(o$ref, "genes.gtf")); truth <- read_truth(file.path(o$ref, "truth.tsv"))
  tx <- lapply(genes, function(g) { s <- paste(substring(genome[[g$chrom]], g$s, g$e), collapse = ""); if (g$strand == "-") revcomp(s) else s })
  dir.create(file.path(o$out, "inputs"), showWarnings = FALSE, recursive = TRUE)
  qchars <- strsplit("+,-./0", "")[[1]]
  for (cond in c("CTL", "TRT")) for (r in seq_len(o$reps)) {
    path <- file.path(o$out, "inputs", sprintf("%s%d_ont.fastq.gz", cond, r)); con <- gzfile(path, "w"); n <- 0
    for (gid in names(tx)) {
      seq <- tx[[gid]]; tr <- truth[[gid]]; if (is.null(tr)) tr <- list(de = "none", fc = 1, base = 120)
      mult <- tr$base / 120
      if (cond == "TRT" && tr$de != "none") mult <- mult * (if (tr$de == "up") tr$fc else 1 / tr$fc)
      for (j in seq_len(max(1, trunc(rnorm(1, o$reads * mult, o$reads * mult * 0.12))))) {
        keep <- runif(1, 0.75, 1.0); start <- if (runif(1) < 0.5) trunc(nchar(seq) * (1 - keep)) else 0   # 5' truncation, mostly full length
        frag <- substr(seq, start + 1, start + max(300, trunc(nchar(seq) * keep)))
        if (runif(1) < 0.5) frag <- revcomp(frag)                   # unstranded cDNA
        frag <- mutate(frag, 0.07)
        qual <- paste(sample(qchars, nchar(frag), replace = TRUE), collapse = "")   # Phred ~10-15
        n <- n + 1
        writeLines(sprintf("@%s%d_ont_%d runid=synthetic ch=%d\n%s\n+\n%s", cond, r, n, sample.int(512, 1), frag, qual), con)
      }
    }
    close(con); cat(path, n, "reads\n")
  }
}

if (sys.nframe() == 0L) main()
