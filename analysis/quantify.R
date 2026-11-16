# analysis/quantify.R -- step 1 when you start from raw reads: FASTQ -> counts table.
# Runs the same Tuxedo pipeline as the TUX tab (fastq-mcf or fastp -> HISAT2 | samtools -> Cufflinks -> Cuffmerge -> Cuffdiff, see
# tuxedo.R) inside <project>/01_quantification/tux/, then writes
#   counts.csv        genes x samples raw fragment counts (what DESeq2 and the other steps read)
#   sample_sheet.csv  sample, condition, time_point, replicate, file (derived from your sheet or the file names)
# and hands both on to the steps that follow. Cuffdiff's own comparison of the groups stays in tux/diff_out/.
strym_step_quantify <- function(ctx) {
  d <- ctx$design
  if (!exists("run_tuxedo", mode = "function")) source(file.path(ctx$app_dir, "tuxedo.R"))
  samples <- strym_quantify_samples(d)
  if (identical(d$read_platform, "long")) {                        # Oxford Nanopore / PacBio: minimap2 + featureCounts instead of Tuxedo
    if (!exists("strym_quantify_long", mode = "function")) source(file.path(ctx$app_dir, "longreads.R"))
    if (!exists("run_tuxedo", mode = "function")) source(file.path(ctx$app_dir, "tuxedo.R"))
    return(strym_quantify_long(ctx, samples))
  }
  tux_out <- file.path(ctx$out_dir, "tux")
  # An earlier run that was cut short (a reboot, a crash) is picked up where it stopped: finished, intact parts are reused (see
  # tux_resume_status). "Start the Tux step from scratch" in Analysis settings turns this off and replaces everything.
  resume <- !isTRUE(strym_normalize_params(d$params)$tux_fresh)
  if (dir.exists(tux_out) && !resume) clear_outputs(tux_out)
  ctx$log(sprintf("Quantifying %d %s samples in %d groups: %s", nrow(samples), d$read_type, length(unique(samples$group)),
                  paste(unique(samples$group), collapse = ", ")))

  labels <- c("cleaning and aligning", "assembling transcripts", "merging assemblies", "differential expression")
  st <- new.env(); st$step <- 1L
  follow <- function(m) {                                          # turn the pipeline's own markers into progress
    txt <- conditionMessage(m)
    if (grepl("^\\[STEP [0-9]/4\\]", txt)) {
      st$step <- as.integer(sub("^\\[STEP ([0-9])/4\\].*", "\\1", txt))
      ctx$sub_progress((st$step - 1) / 4, labels[st$step])
    } else if (grepl("^\\[SAMPLE [0-9]+/[0-9]+\\]", txt) && st$step %in% 1:2) {
      ij <- as.integer(strsplit(sub("^\\[SAMPLE ([0-9]+/[0-9]+)\\].*", "\\1", txt), "/")[[1]])
      ctx$sub_progress(((st$step - 1) + (ij[1] - 1) / ij[2]) / 4, sprintf("%s, sample %d of %d", labels[st$step], ij[1], ij[2]))
    }
  }
  withCallingHandlers(
    run_tuxedo(inputs_folder = d$data_path, gtf_file = d$gtf_file, fasta_file = d$fasta_file, output_folder = tux_out,
               read_type = d$read_type, samples = samples[, c("group", "file", "mate")], threads = d$threads,
               adapter_fasta = if (nzchar(d$adapter_file) && file.exists(d$adapter_file)) d$adapter_file else file.path(ctx$app_dir, "IlluminaAdaptorSeq.fasta"), resume = resume),
    message = follow)

  counts <- strym_counts_from_cuffdiff(tux_out)
  # A sample whose reads were all removed by cleaning, or none of which aligned, gives an all-zero column and a run that looks fine: say so.
  tot <- colSums(counts[-1]); med <- stats::median(tot)
  if (any(tot == 0)) stop(sprintf("No reads were counted for: %s. All of their reads were removed by the cleaning step or did not align to the genome (wrong species or genome file? reads shorter than the minimum length? long-read data?). Check those files and run again.", paste(names(tot)[tot == 0], collapse = ", ")), call. = FALSE)
  if (any(tot < 0.05 * med)) ctx$log(sprintf("WARNING: very few reads were counted for %s (under 5%% of the median sample): their results are unreliable.", paste(names(tot)[tot < 0.05 * med], collapse = ", ")))
  utils::write.csv(counts, file.path(ctx$out_dir, "counts.csv"), row.names = FALSE)
  used <- utils::read.csv(file.path(tux_out, "samples_used.csv"), stringsAsFactors = FALSE)
  m <- match(normalizePath(used$file, mustWork = FALSE), normalizePath(samples$file, mustWork = FALSE))   # by FILE: two samples can share a file name, never a path
  if (anyNA(m)) m[is.na(m)] <- match(used$sample, .stem(samples$file))[is.na(m)]                           # (older runs recorded no usable path)
  ss <- data.frame(sample = used$sample, condition = samples$condition[m], time_point = samples$time_point[m],
                   replicate = samples$replicate[m], file = basename(used$file), stringsAsFactors = FALSE)
  if (all(is.na(ss$replicate))) ss$replicate <- as.integer(ave(seq_len(nrow(ss)), ss$condition, ss$time_point, FUN = seq_along))
  if (all(!nzchar(ss$time_point))) ss$time_point <- if (length(d$timepoint_names) == 1) d$timepoint_names else "t1"
  utils::write.csv(ss, file.path(ctx$out_dir, "sample_sheet.csv"), row.names = FALSE)

  list(status = "ok",
       files = c("counts.csv", "sample_sheet.csv", "tux/diff_out/gene_exp.diff", "tux/samples_used.csv",
                 if (file.exists(file.path(tux_out, "rtpstar_inputs", "rtpstar_genes_from_tux.csv"))) c("tux/rtpstar_inputs/rtpstar_genes_from_tux.csv", "tux/rtpstar_inputs/rtpstar_expression_from_tux.csv")),
       message = sprintf("%s genes × %d samples counted (Cuffdiff comparison of the groups: tux/diff_out/).%s",
                         format(nrow(counts), big.mark = ","), ncol(counts) - 1,
                         if (file.exists(file.path(tux_out, "rtpstar_inputs", "rtpstar_genes_from_tux.csv"))) " A genes file and an expression file for RTP-STAR are in tux/rtpstar_inputs/." else ""),
       data_path = file.path(ctx$out_dir, "counts.csv"), sample_sheet = file.path(ctx$out_dir, "sample_sheet.csv"))
}
