#!/usr/bin/env Rscript
# Generate a small synthetic RNA-seq dataset for testing the TUX pipeline (base R only; deterministic for a given --seed).
#
# Creates: genome.fa, genes.gtf, truth.tsv and
#   <out>/inputs/{control,treated}/rep{1..3}[_1/_2].fastq   (single- or paired-end)
# Genes labelled DE in truth.tsv are up/down-regulated 6-8x in 'treated'.
#
# usage: Rscript make_synthetic.R --out DIR [--seed 1] [--paired] [--chr 2] [--chr-len 60000] [--genes-per-chr 6]
#                                 [--reps 3] [--read-len 75] [--layout folders|flat] [--gz] [--depth 1]

ADAPTER <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
BASES <- c("A", "C", "G", "T")

revcomp <- function(s) paste(rev(strsplit(chartr("ACGTN", "TGCAN", s), "", fixed = TRUE)[[1]]), collapse = "")
rand_seq <- function(n) paste(sample(BASES, n, replace = TRUE), collapse = "")
pick <- function(x) x[sample.int(length(x), 1)]                     # sample() on a length-1 numeric would mean 1:x

make_genes <- function(genome, per_chr, de_fraction) {
  genes <- list(); gid <- 0
  for (chrom in names(genome)) {
    slot <- nchar(genome[[chrom]]) %/% per_chr
    for (i in seq_len(per_chr)) {
      gid <- gid + 1; pos <- (i - 1) * slot + 500
      n_ex <- pick(c(2, 3, 3, 4)); s <- integer(); e <- integer()
      for (k in seq_len(n_ex)) {
        ln <- sample(250:500, 1); s <- c(s, pos); e <- c(e, pos + ln - 1)   # 0-based inclusive
        pos <- pos + ln + sample(150:400, 1)                                # intron
      }
      genes[[gid]] <- list(id = sprintf("G%03d", gid), chrom = chrom, strand = pick(c("+", "-")),
                           s = s, e = e, base = pick(c(30, 60, 120)), de = "none", fc = 1)
    }
  }
  for (i in sample.int(length(genes), max(2, floor(length(genes) * de_fraction)))) {
    genes[[i]]$de <- pick(c("up", "down")); genes[[i]]$fc <- pick(c(6, 8))
  }
  genes
}

transcript_seq <- function(genome, g) {
  s <- paste(substring(genome[[g$chrom]], g$s + 1, g$e + 1), collapse = "")
  if (g$strand == "+") s else revcomp(s)
}

qual <- function(n) {                                               # Phred ~36 falling along the read
  q <- pmax(2, pmin(40, trunc(rnorm(n, 36, 4)) - ((seq_len(n) - 1) * 12) %/% max(n, 1)))
  intToUtf8(33 + q)
}

mutate <- function(s, rate = 0.005) {
  ch <- strsplit(s, "", fixed = TRUE)[[1]]; hit <- runif(length(ch)) < rate
  ch[hit] <- sample(BASES, sum(hit), replace = TRUE); paste(ch, collapse = "")
}

sample_reads <- function(tseq, n, rl, frag_mean = 220) {
  r1 <- character(n); r2 <- character(n)
  for (i in seq_len(n)) {
    flen <- min(max(rl + 10, trunc(rnorm(1, frag_mean, 25))), nchar(tseq))
    st <- sample.int(nchar(tseq) - flen + 1, 1) - 1
    frag <- substr(tseq, st + 1, st + flen)
    if (runif(1) < 0.5) frag <- revcomp(frag)
    a <- substr(frag, 1, rl); b <- substr(revcomp(frag), 1, rl)
    if (flen < rl) a <- substr(paste0(frag, ADAPTER), 1, rl)        # adapter read-through
    r1[i] <- mutate(a); r2[i] <- mutate(b)
  }
  list(r1, r2)
}

write_fastq <- function(path, recs, prefix) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w"); on.exit(close(con))
  q <- vapply(nchar(recs), qual, "")
  writeLines(paste0("@", prefix, ".", seq_along(recs) - 1, "\n", recs, "\n+\n", q), con)
}

main <- function(argv = commandArgs(TRUE)) {
  o <- list(out = "synthetic", seed = 1, paired = FALSE, chr = 2, `chr-len` = 60000, `genes-per-chr` = 6, reps = 3,
            `read-len` = 75, layout = "folders", gz = FALSE, depth = 1)
  i <- 1
  while (i <= length(argv)) {
    k <- sub("^--", "", argv[i])
    if (!k %in% names(o)) stop("unknown option: ", argv[i])
    if (is.logical(o[[k]])) { o[[k]] <- TRUE; i <- i + 1 }
    else { o[[k]] <- if (is.numeric(o[[k]])) as.numeric(argv[i + 1]) else argv[i + 1]; i <- i + 2 }
  }
  stopifnot(o$layout %in% c("folders", "flat"))
  set.seed(o$seed)
  genome <- setNames(lapply(seq_len(o$chr), function(i) rand_seq(o$`chr-len`)), paste0("chr", seq_len(o$chr)))
  genes <- make_genes(genome, o$`genes-per-chr`, 0.34)
  dir.create(o$out, showWarnings = FALSE, recursive = TRUE)

  fa <- character()
  for (cn in names(genome)) { s <- genome[[cn]]; st <- seq(1, nchar(s), by = 60); fa <- c(fa, paste0(">", cn), substring(s, st, st + 59)) }
  writeLines(fa, file.path(o$out, "genome.fa"))

  gtf <- unlist(lapply(genes, function(g) sprintf('%s\tsynthetic\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s.1"; exon_number "%d";',
                                                   g$chrom, as.integer(g$s + 1), as.integer(g$e + 1), g$strand, g$id, g$id, seq_along(g$s))))
  writeLines(gtf, file.path(o$out, "genes.gtf"))
  writeLines(c("gene_id\tchrom\tde\tfold_change\tbase_reads",
               vapply(genes, function(g) sprintf("%s\t%s\t%s\t%s\t%s", g$id, g$chrom, g$de, format(g$fc, nsmall = 1), format(g$base, nsmall = 1)), "")),
             file.path(o$out, "truth.tsv"))

  seqs <- lapply(genes, function(g) transcript_seq(genome, g))
  flat <- o$layout == "flat"; gz <- if (o$gz) ".gz" else ""
  for (cond in c("control", "treated")) {
    d <- if (flat) file.path(o$out, "inputs") else file.path(o$out, "inputs", cond)
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    tag <- c(control = "CTL", treated = "TRT")[[cond]]
    for (r in seq_len(o$reps)) {
      r1 <- character(); r2 <- character()
      for (k in seq_along(genes)) {
        g <- genes[[k]]; mult <- 1
        if (cond == "treated") mult <- if (g$de == "up") g$fc else if (g$de == "down") 1 / g$fc else 1
        mean <- g$base * mult * o$depth * runif(1, 0.85, 1.15)
        n <- max(1, trunc(rnorm(1, mean, sqrt(mean))))
        x <- sample_reads(seqs[[k]], n, o$`read-len`); r1 <- c(r1, x[[1]]); r2 <- c(r2, x[[2]])
      }
      ord <- sample.int(length(r1)); r1 <- r1[ord]; r2 <- r2[ord]
      if (flat) { n1 <- if (o$paired) sprintf("%s%d_Fe_R1.fastq%s", tag, r, gz) else sprintf("%s%d_Fe.fastq%s", tag, r, gz); n2 <- sprintf("%s%d_Fe_R2.fastq%s", tag, r, gz) }
      else      { n1 <- if (o$paired) sprintf("A_rep%d.fastq%s", r, gz) else sprintf("rep%d.fastq%s", r, gz); n2 <- sprintf("B_rep%d.fastq%s", r, gz) }
      write_fastq(file.path(d, n1), r1, paste0(cond, r))
      if (o$paired) write_fastq(file.path(d, n2), r2, paste0(cond, r))
    }
  }
  cat(sprintf("Wrote synthetic data to %s (%d genes, %d DE)\n", o$out, length(genes), sum(vapply(genes, function(g) g$de != "none", NA))))
}

if (sys.nframe() == 0L) main()
