# analysis/registry.R -- the analysis steps StRym can run automatically, and the rules that decide
# when each one applies. This is the ONE place to edit if a rule changes.
#
# A "design" is a list describing the experiment (see pipeline.R: strym_default_design()):
#   data_type ("bulk" | "single_cell"), n_conditions, n_timepoints, de_methods, plots, interactive_browser ...
#
# Each step has:
#   id, label, group     how it is shown in the plan
#   dir                  where its results go, relative to the project folder (numbered, so the folder reads in order)
#   script, fn           analysis/<script> defines the function fn(ctx) that does the work
#   desc                 one-sentence explanation shown to the user
#   rule(design)         -> list(applicable, default, note)
#                         applicable: can this step make sense for this experiment?
#                         default:    should it be ticked automatically, from what the user asked for?
#                         note:       short text explaining the decision (for a step that does not apply: WHY it is greyed out)
#                         fit:        TRUE when this is an especially good choice for the design (not shown in the page any more; kept as data)
#                         why:        one sentence saying what the step is good for
#   needs_any            ids of steps of which at least one must also be selected ("*" = any other step);
#                         the step is "blocked" until that is true
#   needs_unless(design) optional: TRUE when the step does not need those steps after all (e.g. the researcher gave a gene list)

# Why a step that reads Tux (Cuffdiff) results cannot work for this design ("" = it can).
.needs_tux <- function(d) {
  if (identical(d$data_type, "single_cell")) "Not available for single cell/single nucleus data: it reads Tux (Cuffdiff) results from whole organ/tissue RNA-seq."
  else if ((!identical(d$data_kind, "fastq") || identical(d$read_platform, "long")) && !length(d$diff_path)) {              # a table: q-values are calculated with DESeq2 when it holds raw counts with replicates
    why <- strym_derive_possible(d)
    if (nzchar(why)) paste0("Needs q-values for treated vs control (a significance value for every gene). ", why, " Or start from raw reads (FASTQ) so StRym runs the Tux pipeline.")
    else if (d$n_conditions < 2) "Needs a control and a treated condition (at least 2 conditions)."
    else ""
  }
  else if (d$n_conditions < 2) "Needs a control and a treated condition (at least 2 conditions)."
  else ""
}

strym_steps <- function() list(
  list(id = "quantify", label = "Quantify reads (clean, align, count)", group = "Quantification",
       dir = "01_quantification", script = "quantify.R", fn = "strym_step_quantify",
       desc = "Turns your FASTQ reads into a counts table with the Tuxedo pipeline (fastq-mcf or fastp, HISAT2, Cufflinks, Cuffdiff).",
       needs_any = character(),
       rule = function(d) list(
         applicable = identical(d$data_kind, "fastq"),
         default = identical(d$data_kind, "fastq"),
         note = if (identical(d$data_kind, "fastq")) "Runs first: the steps below use its counts table. This can take a while (minutes to hours)."
                else "Only needed when you start from raw reads.")),
  list(id = "deseq2", label = "DESeq2 differential expression", group = "Differential expression",
       dir = "02_differential_expression/deseq2", script = "deseq2.R", fn = "strym_step_deseq2",
       desc = "Genes that change between conditions (negative-binomial model).",
       needs_any = character(),
       rule = function(d) list(
         applicable = d$n_conditions >= 2 && !identical(d$data_type, "single_cell") && !(identical(d$data_kind, "counts") && identical(strym_unit_of(d), "normalized")),
         default = "deseq2" %in% d$de_methods && !identical(d$data_type, "single_cell") && !(identical(d$data_kind, "counts") && identical(strym_unit_of(d), "normalized")),
         fit = d$n_conditions >= 2 && d$data_type == "bulk",
         why = "Compares conditions using raw read counts, so it suits whole organ/tissue experiments with replicates.",
         note = if (identical(d$data_type, "single_cell")) "Not available for single cell/single nucleus data: DESeq2 reads a counts table (genes in rows, one column per sample), not a 10x folder, .h5 file or Seurat object. The Seurat clustering step already compares conditions inside each cluster."
                else if (identical(d$data_kind, "counts") && identical(strym_unit_of(d), "normalized"))
                  { if (identical(d$expr_unit, "normalized")) "DESeq2 needs raw counts, and you chose “Already normalized” for this table. If the numbers are raw counts, choose “Raw counts”."
                    else "DESeq2 needs raw counts, but this table has decimal values, so it looks normalized (FPKM / TPM). If it holds expected counts (RSEM, Salmon), choose “Raw counts”." }
                else if (d$n_conditions < 2) "DESeq2 compares conditions, so it needs at least 2."
                else "")),
  list(id = "shannon_entropy", label = "Shannon entropy", group = "Differential expression",
       dir = "02_differential_expression/shannon_entropy", script = "shannon_entropy.R", fn = "strym_step_shannon_entropy",
       desc = "Finds genes whose expression changes over time under the treatment (Shannon entropy of the time profile) and the time points where treated differs from control.",
       needs_any = character(),
       rule = function(d) { why_not <- .needs_tux(d)
         if (!nzchar(why_not) && d$n_timepoints < 2) why_not <- sprintf("Needs at least 2 time points: it measures how a gene changes over time, and this experiment has %d.", d$n_timepoints)
         list(applicable = !nzchar(why_not),
              default = "shannon_entropy" %in% d$de_methods,
              fit = d$n_timepoints >= 3,
              why = "Scores how evenly a gene is expressed across time points; most informative for a time course (3 or more time points).",
              note = why_not) }),
  list(id = "modified_shannon_entropy", label = "Modified Shannon entropy", group = "Differential expression",
       dir = "02_differential_expression/modified_shannon_entropy", script = "modified_shannon_entropy.R",
       fn = "strym_step_modified_shannon_entropy",
       desc = "Finds genes with outlier time points under the treatment (modified Shannon entropy, MSE) that also differ from the control.",
       needs_any = character(),
       rule = function(d) { why_not <- .needs_tux(d)
         if (!nzchar(why_not) && d$n_timepoints < 3) why_not <- sprintf("Needs at least 3 time points: it looks for outlier time points, and with %d there is nothing to compare them with.", d$n_timepoints)
         list(applicable = !nzchar(why_not),
              default = "modified_shannon_entropy" %in% d$de_methods,
              fit = d$n_timepoints >= 3,
              why = "Scores how a gene’s expression changes across time points, so it fits time-course experiments.",
              note = why_not) }),
  list(id = "seurat_clustering", label = "Seurat clustering: UMAP and t-SNE", group = "Clustering and plots",
       dir = "03_clustering/seurat_umap_tsne", script = "seurat_clustering.R", fn = "strym_step_seurat_clustering",
       desc = "Seurat workflow: quality control, PCA, integration of the samples, clusters, UMAP and/or t-SNE maps, cluster markers and treated-vs-control genes per cluster.",
       needs_any = character(),
       rule = function(d) list(
         applicable = d$data_type == "single_cell",
         default = any(c("umap", "tsne") %in% d$plots),
         fit = d$data_type == "single_cell",
         why = "Groups cells and draws UMAP and/or t-SNE maps: made for single cell/single nucleus data.",
         note = if (d$data_type == "bulk") "Needs single cell/single nucleus data." else "")),
  list(id = "qc_overview", label = "Data overview", group = "Clustering and plots",
       dir = "00_data_overview", script = "qc_overview.R", fn = "strym_step_qc_overview",
       desc = "A quick look at your data: samples per condition and time point, sequencing depth, genes detected, and how alike the samples are, with any problems written out in plain words.",
       needs_any = character(),
       rule = function(d) list(applicable = identical(d$data_type, "bulk"), default = "qc" %in% d$plots && identical(d$data_type, "bulk"),
         fit = identical(d$data_type, "bulk"), why = "Shows at a glance whether the experiment is complete and the samples comparable.",
         note = if (identical(d$data_type, "single_cell")) "For single cell/single nucleus data the Seurat step draws its own quality plots." else "")),
  list(id = "pca", label = "PCA", group = "Clustering and plots",
       dir = "04_pca", script = "pca.R", fn = "strym_step_pca",
       desc = "Principal component analysis of the samples: shows how the replicates group and which look like outliers (or, in Analysis settings, one point per condition and time point averaged over the replicates).",
       needs_any = character(),
       rule = function(d) list(applicable = identical(d$data_type, "bulk"), default = "pca" %in% d$plots && identical(d$data_type, "bulk"),
         fit = identical(d$data_type, "bulk"), why = "A quick overview of how your samples group by condition and time point.",
         note = if (identical(d$data_type, "single_cell")) "For single cell/single nucleus data the PCA is part of the Seurat clustering step (its elbow plot)." else "")),
  list(id = "trajectory", label = "Trajectory plots", group = "Clustering and plots",
       dir = "05_trajectories", script = "trajectory.R", fn = "strym_step_trajectory",
       desc = "Expression over time for chosen genes, treated vs control (line = mean, points = replicates).",
       needs_any = c("shannon_entropy", "modified_shannon_entropy"),
       needs_unless = function(d) nzchar(strym_normalize_params(d$params)$trajectory_genes),           # a gene list was given
       rule = function(d) list(
         applicable = identical(d$data_type, "bulk") && d$n_timepoints >= 3,
         default = "trajectory" %in% d$plots && identical(d$data_type, "bulk"),
         fit = identical(d$data_type, "bulk") && d$n_timepoints >= 3,
         why = "Shows how expression moves across time points, for the genes an entropy step finds or for your own gene list.",
         note = if (identical(d$data_type, "single_cell")) "Trajectory plots use time-course expression from whole organ/tissue RNA-seq."
                else if (d$n_timepoints < 3) sprintf("Trajectories need at least 3 time points; this experiment has %d.", d$n_timepoints)
                else "Plots the genes found by an entropy step (tick one) or the gene list chosen in Analysis settings.")),
  list(id = "genist_network", label = "GENIST regulatory network", group = "Gene regulatory networks",
       dir = "06_regulatory_networks/GENIST_auto", script = "genist_network.R", fn = "strym_step_genist_network",     # not GENIST/: that is the GENIST tab's own results folder, and a run clears its step folder
       desc = "Builds the genes file and the time course file for GENIST from the entropy results (for every treated condition) and runs GENIST on them.",
       needs_any = c("shannon_entropy", "modified_shannon_entropy"),
       rule = function(d) { why_not <- .needs_tux(d)
         if (!nzchar(why_not) && d$n_timepoints < 3) why_not <- sprintf("Needs a time course of at least 3 time points; this experiment has %d.", d$n_timepoints)
         list(applicable = !nzchar(why_not),
              default = "genist" %in% d$grn_methods,
              fit = d$n_timepoints >= 3,
              why = "Infers a regulatory network from the genes the entropy steps find, over the time course.",
              note = if (nzchar(why_not)) why_not else "Uses the main-set genes of the entropy step(s) you tick; the run time grows steeply with the number of genes (500 at most by default).") }),
  list(id = "web_browser", label = "Interactive web browser", group = "Share",
       dir = "07_interactive_browser", script = "web_browser.R", fn = "strym_step_web_browser",
       desc = "A small web app to explore the plots and tables produced above.",
       needs_any = "*",
       rule = function(d) list(applicable = TRUE, default = isTRUE(d$interactive_browser),
                               note = "Built from the other results: tick at least one analysis step."))
)

# Network modules that already exist as their own tabs. They need input files only the researcher can
# choose (gene lists, time courses, ...), so the plan points to them instead of running them blindly.
strym_modules <- function() list(
  list(id = "tux", label = "TUX: reads → counts → differentially expressed genes", tab = "TUX · RNA-seq",
       desc = "Cleans, aligns and quantifies raw reads (HISAT2 / Cufflinks / Cuffdiff).",
       rule = function(d) list(applicable = identical(d$data_kind, "fastq"),
                               note = if (identical(d$data_kind, "fastq")) "StRym runs this for you as step 1 (“Quantify reads”); this tab offers the same pipeline with every parameter."
                                      else "Only needed if you start from raw reads.")),
  list(id = "genist", label = "GENIST: regulatory network from a time course", tab = "GENIST",
       desc = "Bayesian network inference from time-course data.",
       rule = function(d) list(applicable = d$n_timepoints >= 3, fit = d$n_timepoints >= 3,
                               why = "Builds a network from a time course (3 or more time points).",
                               note = if (d$n_timepoints >= 3) "Your time points allow a time-course network."
                                      else "Needs a time course of at least 3 points.")),
  list(id = "rtpstar", label = "RTP-STAR: directed regulatory network", tab = "RTP-STAR",
       desc = "Random-forest network with clustering and edge direction.",
       rule = function(d) list(applicable = TRUE, fit = TRUE,
                               why = "Works with replicates; a time course adds the direction of each edge.",
                               note = "Works with replicates; a time course adds edge direction."))
)
