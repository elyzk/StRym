# StRym: notes for whoever continues this work

StRym (Streamlined transcriptomic Regulatory mapping) is an R/Shiny browser application, run in WSL, that ports the TuxNet modules (TUX = the
Tuxedo RNA-seq pipeline, GENIST, RTP-STAR) and adds an "Experimental Design" page that plans and runs analyses into a numbered project folder.
It was built for biologists: plain-language messages, checks and a reason on every grayed-out option. Developed by Jordan Reed and
Stephanie M. Stewart; maintained by the Sozzani lab at NC State University (this text is also in the footer of every page).

## Map of the folder
| Path | What it is |
|---|---|
| `app.R` | the Shiny app: navbar (About, Experimental Design, Analysis modules > TUX / GENIST / RTP-STAR, Drives), footer, the TUX tab |
| `mod_design.R` | the Experimental Design page (form, sample sheet, options, plan, confirmation window, results, clear, sample experiment) |
| `pipeline.R` | the engine: design -> plan -> run; analysis settings (`strym_default_params`), sample sheet building, run estimate, checks |
| `analysis/registry.R` | the steps and the **rules** that decide when each applies (with the reason shown when grayed out) |
| `analysis/common.R` | shared analysis code: reading Tux (Cuffdiff) output, treated-vs-control pairs, Shannon / MSE algorithms |
| `analysis/{quantify,deseq2,pca,shannon_entropy,modified_shannon_entropy,trajectory,seurat_clustering,web_browser}.R` | one file per step, `strym_step_<id>(ctx)` |
| `analysis/baseline_scripts/` | the lab's ORIGINAL scripts (unchanged) the steps were generalized from, used as the reference in tests |
| `tuxedo.R` | the Tuxedo pipeline (fastq-mcf, HISAT2, samtools, Cufflinks, Cuffmerge, Cuffdiff) incl. space-safe paths and **resume** |
| `genist.R`, `rtpstar.R`, `mod_genist.R`, `mod_rtpstar.R` | the GENIST and RTP-STAR ports (checked against GENIE3, not yet numerically against the MATLAB originals) |
| `gui_helpers.R`, `jobs.R`, `drives.R`, `io_utils.R` | path fields with drag and drop, background jobs, WSL drives, readers |
| `www/strym.css`, `www/strym.js`, `www/leaf.svg` | theme and the drag-and-drop script |
| `www/strym-logo.png`, `www/strym-logo-white.png`, `www/favicon.png`, `www/favicon-192.png` | the StRym logo (colour for white pages, white for the dark header and banner), tab icon (the "S"); the white logo is also embedded in `index.html` and `strym_browser.html` by `strym_logo_img()` |
| `run_analysis_cli.R`, `run_cli.R`, `run_genist_cli.R`, `run_rtpstar_cli.R` | command-line runners (the GUI starts these as background jobs) |
| `tests/` | `bash tests/run_tests.sh [name]`; files `test-*.R`, helpers `helper-*.R` (synthetic Tux, 10x data) |
| `setup_wsl.sh`, `setup_modules.sh`, `environment*.{yml,txt}` | install; exact versions recorded from the working machine |

## Conventions the team asked for
* **American spelling** everywhere. **NC State palette**: Genomic Green `#6F7D1C`, Bio-Indigo `#4156A1`, Innovation Blue `#427E93`,
  Carmichael Aqua `#008473`, Wolfpack neutrals; **no warm colors, except red for errors and warnings** (`#CC0000` / `#990000`).
* Options that cannot work for the design are **grayed out with the reason**; options that work are plain (no "best fit" badges).
* Read-only inputs: StRym never changes the user's input files. Results go to the output folder only.
* Every run goes through a confirmation window (disk space, time estimate + 10%, sample sheet summary).
* The sample sheet is built from the inputs and shown for confirmation; saved as `sample_sheet_<project>.csv` in the output folder.

## How the analysis steps get their data
StRym's TUX step runs ONE Cuffdiff over all groups (group = `<condition>` or `<condition>_<time point>`). Steps read
`01_quantification/{counts.csv, sample_sheet.csv, tux/samples_used.csv, tux/diff_out/{genes.read_group_tracking, gene_exp.diff}}`.
Treated vs control pairs come from the condition names (or Analysis settings). The entropy steps need Tux results (raw reads) and >= 2 / >= 3
time points; trajectories need >= 3; PCA and DESeq2 work from a counts table too. Single cell: Seurat step (10x folders or `.rds`).

## Pitfalls found the hard way
* `jsonlite::validate` masks `shiny::validate`: use `shiny::validate` / `shiny::need`.
* Calling `wsl ...` from Git Bash rewrites `/mnt/c/...` arguments: set `MSYS_NO_PATHCONV=1`, or put the commands in a script file.
* The hisat2 launcher breaks on spaces in paths: every path handed to a tool goes through `.safe()` (symlinks with plain names).
* A drive plugged in after WSL started is not mounted: the Drives tab mounts it. shinyFiles roots must be a fixed vector.
* Do not let conda "upgrade" the Tuxedo tools: the pins in `setup_*.sh` keep samtools at 1.16.1 (otherwise it downgrades to 1.5).
* Seurat's `Read10X` only reads the all-gz layout; StRym has its own reader (`.read_10x`). CCA/RPCA need `k.weight` <= the smallest sample's cells.
* Cuffdiff writes its result tables only at the very end and prints no progress in quiet mode; its testing phase is long.
* A Windows restart (the lab PC has a scheduled "weekly reboot" at 03:00 Friday) kills WSL and any running job. That is what the resume feature is for.
* Background-job state lives in `~/.tuxnet_r` (WSL): `analysis_state.json`, `analysis_last_run.log` (progress markers `[PROGRESS n]`, `[STEP]`, `[SAMPLE]`, `DONE:` / `ERROR:`).

## Real-data facts (4 single-end samples, ~80 million reads, TAIR10, 1 thread; use to recalibrate `strym_estimate_run`)
* hisat2-build ~2 min; cleaning + alignment ~40-55 min per sample (~138 s per million reads); ~97% alignment.
* **Cufflinks with the real annotation: ~26 min per sample** (31,356 loci, bias correction) -- NOT the ~3 min per sample assumed from the first test
  run (that run used a wrong annotation and found 0 loci). Cuffmerge ~5 min. Cuffdiff: preliminary estimates ~20 min, bias ~9 min, then testing
  for differential expression for **more than 1 h 40 min without finishing** (unknown total). `strym_estimate_run` still uses the old, too-small Cufflinks/Cuffdiff constants.

## State when this zip was made
* All features of the three revision batches (28 items) and the later additions (baseline analysis scripts, DESeq2, Seurat, analysis settings,
  resume of interrupted Tux runs) are in. Run the tests (`bash tests/run_tests.sh`) to see the current state.
* The real full run (4 FASTQ files from `D:\TEST FOR TUX`, conditions B and D, 1 thread, TAIR10) was cut short at Cuffdiff by the weekly reboot.
  Aligned BAMs, Cufflinks output and `merged.gtf` were intact; the resume feature is meant to continue from there (project output folder was
  `imani 2025 nat comm full run test`, not included in this zip). With 2 conditions at one time point only Quantify, DESeq2, PCA and the browser apply;
  MSE, Shannon entropy, trajectories and GENIST are correctly grayed out.

* Second real run (completed): E-MTAB-15703, 15 single-end FASTQ files (33 GB), 3 stages (earlyInit / midProto / lateSE; one condition, `-Fe`), 2 threads,
  7 GB WSL memory. 25.3 h in total: cleaning + alignment about 15 h (~60 min per sample), Cufflinks ~14.5 min per sample, Cuffmerge ~10 min,
  Cuffdiff ~6.5 h (peak ~2 GB). 33,201 genes counted; PCA separates the stages cleanly (PC1 76.6%). Output went to a folder on the D: drive (~22 GB).
  DESeq2 / MSE / Shannon are grayed out for one condition (they need control vs treated), so GENIST gets its gene list from Cuffdiff's table instead.
* GENIST tab: "Start from your Tux results" (`genist_inputs_from_tux()` in `genist.R`, tests in `test-genist-tux.R`) builds the genes file (genes with a
  significant, large change between two of the stages) and the time course file (FPKM per stage, in the design's time order, NOT Cuffdiff's order)
  from a project folder or a `diff_out` folder. On the real data: 9,379 genes pass q <= 0.05, |log2FC| >= 1, FPKM >= 1; 500 genes take ~2 min in GENIST
  (2 threads; 990 edges, FRO2 and LBD2 among the genes), but 1,500 genes had not finished after 30 minutes: the cost grows steeply with the number of
  genes and of TFs among them (candidate regulator sets of up to 3), so the tab's default limit is 500 genes. Parameter defaults on the tab now follow
  the original TuxNet GUI (fold change 1.1, time percent 0.1, 3 levels, bottom 0.05), and the TF list and gene names are prefilled as in the original.

* Single cell (from a collaborator's error report on the EBI/GEO Shahan root atlas matrix: 24,826 genes x 361,743 cells, 72 M entries, 19 sample labels in
  the barcodes, but only 3 samples with any counts, a 1e308 entry and ~777,000 fractional values): `analysis/singlecell.R` (new: `.find_10x`, barcode sample
  labels, `strym_10x_samples`) and `analysis/seurat_clustering.R` (`.read_mtx` low-memory reader, `.clean_counts`, `.sc_prepare` = clean + split by sample sheet +
  UMI filter before Seurat objects, skipped samples reported). DESeq2 is grayed out for single cell data (it reads a counts table, hence the old "cannot open
  the connection"). The design page drafts the sample sheet from the barcode labels. Tests: `test-seurat-multi.R`.
* Drives tab (`drives.R`, `app.R`): states ready / not mounted / needs reconnecting (stale: the drive was ejected and replugged, Ubuntu kept a dead connection) /
  removed; Mount drive reconnects a stale drive; Release drive unmounts before an eject (refused while any run is going); the Browse windows get a FUNCTION for
  `roots`, so a drive mounted later appears without reloading the page; written directions on the page. Tests: `test-drives.R` (Windows and the mounts are faked).

* Pitfall found on a Mac: the program folder was "StRym 3" (macOS adds " 3" to a second download). `Rscript` writes a space in its own path and
  arguments as `~+~` in `commandArgs()`, so the background runners (`run_analysis_cli.R`, `run_genist_cli.R`, `run_rtpstar_cli.R`, `run_cli.R`,
  `make_demo_data.R`, `tests/run_tests.R`) could not find `io_utils.R` ("cannot open the connection", path with `StRym~+~3`). They now decode `~+~`
  back to a space (`test-spaces-runner.R` runs them from a folder named "StRym 3"). macOS itself is still untested.

* Shahan root atlas (E-MTAB / GEO matrix, 361,743 cells, only 3 of 19 samples with counts) run end to end through `run_analysis_cli.R`: 22,729 cells, 3 samples, CCA, 30 clusters,
  UMAP, markers, interactive browser: 22 min, peak 7.3 GB of the 8 GB Ubuntu has (output in a `StRym_results\shahan_sc` folder in the researcher's Windows home folder; not included in the zip). Three memory problems had to be fixed to get there
  (each run was killed by the kernel with signal 9): (a) `JoinLayers()` on the whole object needs several GB more than it holds, so the marker tests now get a small object built from the
  normalized data and the saved object is joined only when `MemAvailable` allows (otherwise it keeps one layer per sample and the log says so); (b) `FindMarkers(<object>)` turns the whole
  normalized matrix dense (3.7 GB for 20,000 genes x 22,000 cells) whatever `max.cells.per.ident` says, so `.de_wilcox()` hands it only the compared cells (<= 1,000 per side);
  (c) the Seurat object is saved BEFORE the marker tests and released. The log now has "Memory in use ..." lines at each stage (`.mem_note`). A larger dataset than this still needs more
  memory than the default 8 GB: give Ubuntu more in `.wslconfig` (`memory=12GB`) or use RPCA / Harmony and fewer cells.
* Outside review, 11 problems (all reproduced first, then fixed; `test-review-fixes.R` has one test per item):
  (1) Seurat `.rds` without a `sample` column crashed, and conditions were assigned per cell instead of per sample; (2) Cell Ranger `Sample/outs/...` folders were
  all named "outs", and raw + filtered matrices were both read; (3) paths with `( ) + [ ]` found 0 samples (regex built from the path); (4) `write.table` files (header one
  field shorter) lost their gene ids in `read_table_any`; (5) Shannon / MSE crashed when no gene was called; (6) two pairs with the same treated condition
  overwrote each other (results now keyed `<treated>_vs_<control>` only when the treated condition repeats, so single-pair file names are unchanged);
  (7) the same FASTQ name in different condition folders gave duplicate sample names and wrong counts / conditions (names get the group in front, the
  sample sheet records the path relative to the data folder, quantify matches by file); (8) text / annotation columns (gene_name, featureCounts Chr/Start/End/Strand/Length)
  of a counts table became fake samples (`strym_count_columns`); (9) CCA / RPCA crashed on samples under 30 cells (dims limited, then Harmony, then merge only);
  (10) a symbol column in a GENIST / RTP-STAR time-course file counted as a time point (`table_to_matrix` drops non-numeric columns and says so);
  (11) an output name equal to an input (genes.csv, clusters.xlsx ...) could overwrite it (`output_overwrites_input`, checked on the page and in the runners).

## Open items / ideas
1. Recalibrate the time estimate with the real-data facts above (Cufflinks ~14.5 min/sample and Cuffdiff ~6.5 h for 15 samples at 2 threads).
2. The "select genes" step is built for GENIST only (see above). Still open: the same for RTP-STAR, TF flags and a replicate-level FPKM table (TuxOP style).
   The original MATLAB TuxOP files are in the TuxNet download (`TuxNet-MATLAB/TuxOP`).
3. GENIST / RTP-STAR numeric comparison with the MATLAB originals (MATLAB R2025a was available on the lab PC).
4. Entropy methods and trajectories on a counts table (they need Cuffdiff q-values now).
5. macOS / native Linux support is untested.

## Version 7 (overnight work, 10 November 2026): decisions worth knowing
* **Why q-values from DESeq2:** the entropy steps only needed "a q-value for treated vs control per time point". Making that source-agnostic (Tux, Cuffdiff files, or DESeq2 on counts) removed the dependence on FASTQ and made long reads / any provider's counts usable. The result is cached in `<project>/02_differential_expression/.derived_diff_cache.rds`, keyed on the counts, sample sheet, pairs and time points.
* **Expression data is a matrix:** `strym_load_counts()` returns `$wide` (genes x samples); the long `$rep` table now exists only for Tux runs. Use `strym_wide()`, `strym_long(data, genes)` and `strym_gene_ids()`.
* **GENIST step folder is `GENIST_auto`** because a run clears the folder of each step it runs and `GENIST/` is the GENIST tab's own results folder.
* **Fast GENIST inference** (`genist_infer_fast`) was written by an audit agent and verified identical on 48 random cases plus by me on real data (md5 of the Cytoscape table, 1 and 4 cores, time lapse 1 and both). It depends on R summing in long double in index order.
* **Long reads** use their own environment (`longread`) so the Tuxedo environment is never solved again.
* **Tests:** new files test-derived, test-robust, test-reads-robust, test-longreads, test-diffinput. The long-read test skips when the tools are missing; the zip test skips without a `zip` binary.
* Open ideas: strand-aware counting for direct RNA; merging lane-split files automatically (it only warns now); caching parsed Cuffdiff files between steps; a cluster-scheduler front end (the CLI `run_analysis_cli.R` is the building block).