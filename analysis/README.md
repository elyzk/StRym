# Analysis steps in StRym

Each automatic step is one file in this folder, and all of them are real. To change a step, replace the body of its function:
nothing else in StRym has to change. (`placeholder.R` is a helper that writes a labeled note and picture, for steps added later.)

| step id | file | function | status | what it does |
|---|---|---|---|---|
| `deseq2` | `deseq2.R` | `strym_step_deseq2(ctx)` | **real** | DESeq2 (Wald test) for every treated-vs-control comparison at every time point: results, DE genes, normalized counts, volcano plot, heatmap (raw counts, ≥ 2 replicates) |
| `shannon_entropy` | `shannon_entropy.R` | `strym_step_shannon_entropy(ctx)` | **real** | genes whose treated profile changes over time + the time points where treated differs from control (Tux results, ≥ 2 time points) |
| `modified_shannon_entropy` | `modified_shannon_entropy.R` | `strym_step_modified_shannon_entropy(ctx)` | **real** | MSE outlier time points per gene, with the control-drift filter (Tux results, ≥ 3 time points) |
| `seurat_clustering` | `seurat_clustering.R` | `strym_step_seurat_clustering(ctx)` | **real** | Seurat v5: QC, PCA, integration, clusters, UMAP / t-SNE, cluster markers, treated-vs-control genes per cluster |
| `qc_overview` | `qc_overview.R` | `strym_step_qc_overview(ctx)` | **real** | design grid, depth, genes detected, sample correlation and plain-language notes about the data |
| `genist_network` | `genist_network.R` | `strym_step_genist_network(ctx)` | **real** | builds the genes + time-course files for GENIST from the entropy results (for every treated condition) and runs GENIST |
| `pca` | `pca.R` | `strym_step_pca(ctx)` | **real** | replicate-level PCA (Tux FPKM, or counts per million) |
| `trajectory` | `trajectory.R` | `strym_step_trajectory(ctx)` | **real** | expression over time, treated vs control, for the entropy genes or your gene list |
| `web_browser` | `web_browser.R` | `strym_step_web_browser(ctx)` | **real** | `strym_browser.html` (a self-contained export of the plots and tables) + a live Shiny app (template) |

The real steps were generalized from the baseline scripts kept unchanged in `baseline_scripts/` (see the README there). DESeq2 follows
the structure of the lab's reference script (`easy_deseq2`, linked in the supplied document) with explicit contrasts and without its
human-only enrichment.
Shared code (reading Tux results, the treated-vs-control pairs, the entropy algorithms) is in `common.R`; the analysis settings
and their defaults (`strym_default_params()`) are in `../pipeline.R`.

## What the Tux-based steps read
StRym's Tux step runs ONE Cuffdiff over all groups, where a group is `<condition>` (one time point) or `<condition>_<time point>`.
The steps read, from `01_quantification/`: `sample_sheet.csv` (sample, condition, time_point, replicate), `tux/samples_used.csv`,
`tux/diff_out/genes.read_group_tracking` (FPKM per replicate) and `tux/diff_out/gene_exp.diff` (q-values for every pair of groups).
They compare each **treated** condition with its **control** at every time point. The control is detected from the condition names
(control, ctrl, mock, untreated, WT, ...) or set in *Analysis settings*; with several controls, each treated condition is paired with
the control that shares its name (`Pana_Heat` with `Pana_Control`), or set the pairs yourself (`Pana_Heat:Pana_Control, ...`).
A counts table has no Cuffdiff comparisons: PCA and trajectories work on it (library-size normalized), the entropy steps do not.

## The contract
Every step function receives one list, `ctx`, and returns `list(status, files, message)`:

```r
ctx$design        # the experiment: data_type ("bulk"/"single_cell"), data_kind, n_conditions, n_timepoints,
                  #   condition_names, timepoint_names, de_methods, plots, interactive_browser, params (analysis settings), ...
ctx$data_path     # the expression data the user pointed to ("" if none)
ctx$sample_sheet  # sample sheet (condition, time_point, ...) ("" if none)
ctx$out_dir       # write ALL results here (already created; numbered folder inside the project)
ctx$project_dir   # the whole project folder
ctx$prereq        # named list: ids of earlier steps in this run -> their output folders
ctx$log(msg)      # print a progress line to the run log

# return value
list(status  = "ok",                       # "ok" or "placeholder"; to fail, call stop("reason")
     files   = c("results.csv", "plot.png"),   # files written into out_dir (names only)
     message = "short note shown in the summary")
```
If a step calls `stop()`, StRym records the error, skips steps that depend on it and carries on with the rest.

## When does a step apply?
The rules (which steps are offered, ticked by default, or blocked) live in `registry.R`: one small `rule(design)`
function per step, with the reason shown when an option is grayed out. Edit them there; they drive the plan on the
*Experimental Design* page. A step can also say it does not need its prerequisites after all (`needs_unless`): trajectories
need an entropy step *unless* a gene list is given.

## Folder layout produced
```
<project>/
  00_design/            design.json, plan.csv, sample sheet
  01_quantification/    counts.csv, sample_sheet.csv, tux/ (when starting from raw reads)
  02_differential_expression/{deseq2,shannon_entropy,modified_shannon_entropy}/
  03_clustering/seurat_umap_tsne/
  04_pca/
  05_trajectories/
  06_regulatory_networks/   (suggested place for GENIST / RTP-STAR outputs)
  07_interactive_browser/
  run_summary.csv   index.html   run_log.txt
```

## Shared code added in version 7
* `derived_diff.R` -- q-values for the time-course steps from RAW counts (DESeq2 per treated/control pair and time point), cached in the project; `strym_ensure_diff()`.
* `common.R` -- `.counts_raw` (sample-sheet matching, junk cells, duplicates), `strym_load_diff` (Cuffdiff files of any packing), time-point ordering with units.
* `../longreads.R` -- the long-read quantification route (minimap2 + featureCounts).
