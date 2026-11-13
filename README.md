<img src="www/strym-logo.png" alt="StRym" width="260">

# StRym — Streamlined transcriptomic Regulatory mapping

One browser app that goes from RNA-seq reads to regulatory networks. StRym is an R re-implementation of the
TuxNet modules (https://github.com/rspurney/TuxNet), running on Windows through WSL or directly on macOS:

| Module | What it does | Original |
|---|---|---|
| **TUX** | RNA-seq reads → differentially expressed genes (fastq-mcf, HISAT2, SAMtools, Cufflinks, Cuffmerge, Cuffdiff) | TUX tab |
| **GENIST** | gene regulatory networks from time-course data (Bayesian scoring, optional spatial clustering) | GENIST tab |
| **RTP-STAR** | directed regulatory networks with random forests (GENIE3), clustering and time-course direction | RTP-STAR tab |

Not (yet) ported: TuxNet's **TuxOP** step, which turns Cuffdiff output into the DEG gene lists and expression tables
that GENIST and RTP-STAR take as input. For now, supply those files yourself (Excel/CSV).

## What is new in version 8 (fixes found by running 40 test experiments)
* Excel tables with a blank row or several title lines above the header are read correctly.
* DESeq2 no longer stops on a very small gene panel (falls back to one common dispersion value).
* Very large tables (tens of thousands of genes, 100+ samples): the entropy steps use fewer parallel workers and repeat any comparison that lost its worker.
* Long reads: a cut-off `.gz` file is named before alignment, as for short reads.
* A run started without the page (command line, scheduler) refuses a drive root, a file, or the input folder as the output folder.
* When quantification fails, the later steps are skipped with a plain message instead of confusing errors.
* The Python example-data scripts are replaced by R (`make_synthetic.R`, `make_longread_data.R`).

## What is new in version 7
* **Entropy steps from any counts table with replicates** (q-values calculated with DESeq2, cached in the project), from an FPKM table plus Cuffdiff files (any packing: files, folders, zips inside folders inside zips; labels matched by name, folder or expression values), or from long-read counts.
* **Long reads** (Nanopore, PacBio, Plasmidsaurus FASTQ): *Sequencing platform: Long reads* uses minimap2 + featureCounts (`longreads.R`, `setup_modules.sh` installs the tools into their own environment).
* **Data overview** step (design grid, depth, genes detected, sample correlation, plain-language notes) and an automatic **GENIST** step (genes from the entropy results, time course from your data, optional TF list; results in `06_regulatory_networks/GENIST_auto`).
* **GENIST inference about 18x faster**, byte-identical output (`options(strym.genist_orig = TRUE)` runs the original loop).
* **Robust to real data**: samples not in the sheet, junk cells, Inf/negatives, duplicate or empty gene ids, empty libraries, missing time points, semicolon/decimal-comma tables, Excel title rows, time labels with units; read checks (damaged `.gz`, mates that do not belong together, lanes, long reads on the short-read route, Phred offset, empty or tiny samples); a quantification that yields an all-zero sample is an error, not a quiet success.
* **Lighter and faster**: expression values held as one genes x samples matrix (about half the memory), `--no-update-check` for Cufflinks/Cuffdiff (about 17 s saved per call), capped sorter memory, non-blocking drive list, cached table checks.
* **Safer page**: guarded checks, a second run cannot be started while one is going, a crashed child shows as failed (not cancelled), clearing the form resets every field.
* **Command line:** `bash run_analysis.sh config.json` runs the design page's analysis without a browser (it sets the R path and calls `run_analysis_cli.R`; `config.json` = `{"design": {...00_design/design.json...}, "steps": ["deseq2", "pca"]}`). Prints `DONE: <project>/index.html` (exit 0) or `ERROR: ...` (exit 1).
* Full protocol with screenshots: `StRym_Protocol.docx`. Tests: `tests/testthat` now 1,100+ expectations.

## Quick start
1. One-time setup:
   * **Windows 10/11:** run `wsl --install -d Ubuntu`, then, from this folder, `wsl -d Ubuntu -u root -e bash setup_wsl.sh`.
   * **macOS:** install Apple's command-line developer tools with `xcode-select --install`. On Apple Silicon, also install Rosetta 2 with
     `softwareupdate --install-rosetta`. Then run `bash setup_wsl.sh` from this folder. The script name is retained for compatibility,
     but it detects macOS and creates the correct environments.
   * Already set up before GENIST/RTP-STAR existed? Run `setup_modules.sh` in the same way to add the analysis packages.
   * **Docker on an ARM host:** the Ubuntu container must be created with `--platform linux/amd64`; the legacy
     Cufflinks/ea-utils toolchain has no Linux ARM packages. See `MOVING_TO_ANOTHER_MACHINE.md` for a complete command.
2. Start StRym: double-click **`start_gui.bat`** on Windows, or run `bash start_gui.sh` on macOS/Linux, then open
   http://localhost:3838. The **About** page opens first.
3. **Experimental Design** page: describe your experiment and what you want out of it; StRym plans the analysis and
   saves everything in an organized project folder (see below). **Download the sample experiment and fill everything in**
   saves a small made-up study in a `sample_experiment` folder next to StRym and fills in every field for you.
4. Prefer to do things yourself? Every module is under **Analysis modules** with all its parameters: work left to
   right through *Inputs → Parameters → Run → Results*. GENIST and RTP-STAR have a **Load the Arabidopsis example data** button.
5. On Windows/WSL, use the **Drives** tab to mount or release external drives. macOS volumes are available normally under `/Volumes`.
6. Runs happen in the background: close the page and reopen it later to see progress; *Cancel* stops everything.

## Experimental Design (automatic analysis)
You answer: type of data (whole organ/tissue RNA-seq, or single cell/single nucleus RNA-seq), number of conditions, number of time points, differential
expression method(s) (DESeq2, Shannon entropy, modified Shannon entropy), plots (PCA, Seurat UMAP and/or t-SNE, trajectories),
gene regulatory network (GRN) inference (GENIST, RTP-STAR) and whether you want an interactive web browser — plus where your
input data are and which output folder StRym should save the results in. In **What would you like?** every option that cannot
work for your design is grayed out with the reason (e.g. Seurat clustering needs single cell/single nucleus data; modified
Shannon entropy needs ≥ 2 time points); everything that works is simply available. The rules live in `analysis/registry.R`.
**Clear all fields / start a new project** empties the form (nothing on disk is deleted). StRym then:
* shows the **checks** (what is missing or inconsistent, e.g. a sample sheet that disagrees with the number of conditions);
* builds a **plan**: each step is ticked or not from your answers and explains why (e.g. trajectories need ≥ 3 time points
  and build on an entropy or clustering result). Tick/untick anything; **Run only this** runs a single step. GRN methods you
  tick appear in the plan with an **Open the … tab** button, because they need the gene list and expression files you choose there;
* runs the steps and writes an organized **project folder**:
```
<project>/
  index.html  run_summary.csv                      summary page and table
  00_design/            design.json, plan.csv, sample sheet
  01_quantification/    only for raw reads: counts.csv, sample_sheet.csv, tux/ (the Tuxedo run)
  02_differential_expression/{deseq2, shannon_entropy, modified_shannon_entropy}/
  03_clustering/seurat_umap_tsne/     04_pca/     05_trajectories/
  06_regulatory_networks/             (suggested place for GENIST / RTP-STAR outputs)
  07_interactive_browser/             strym_browser.html: ONE self-contained file (plots + tables inside) to open, export and share;
                                      plus app.R, a live Shiny version (template) you can launch from the page
```
**Starting from raw reads:** choose *Raw reads (FASTQ files)* under "What do you have?", then under *Input data* select all
your FASTQ files **at once** (Ctrl- or Shift-click in the file window, or drag and drop them onto the field) — or choose the
folder that holds them (only the files you select are used, even if the folder holds others). Give the annotation (GTF/GFF) and
genome (FASTA): there is a one-click link for the Arabidopsis reference. StRym adds **Quantify reads** as step 1: it runs the same
Tuxedo pipeline as the TUX tab, writes `01_quantification/counts.csv` and a derived `sample_sheet.csv`, and every later step
then works from those counts (later single-step runs reuse them). Quantification runs in the background and can take minutes to
hours depending on the reads.

**Sample sheet.** Section 3 builds the sample sheet for you from your inputs (the real FASTQ names, or the columns of your counts
table; conditions and time points are matched using the names you typed). It is shown right away so you can confirm that StRym
understood the experiment, and any cell can be corrected by double-clicking it. **Save to the project folder** writes it as
`sample_sheet_<project name>.csv` in the output folder (it is also saved there automatically when you run), and **Download a copy**
gives it the same name. Already have a sheet (`condition`, `time_point`, and for raw reads `file` columns)? Use it instead.

**Before a run starts**, a window summarizes the experiment, the steps, the disk space needed against what is free, the estimated
time (and the estimate + 10%), and the sample sheet, and asks you to confirm. The time estimate is calibrated on a real run on this
computer (4 samples, 79 million reads: 3 h 17 min on 1 thread) and scales with the read count (estimated from file sizes) and threads.

**Drag and drop** works on every file and folder field. A browser only reveals a dropped file's name and size, so StRym looks for the
file on this computer (the folder you are already using, Downloads/Desktop/Documents, your profile, every mounted drive). If it
cannot find a small file it takes a copy; large files (FASTQ) are never copied, so use *Browse* for those if they are not found.

**The analysis steps.** DESeq2, PCA, Shannon entropy, modified Shannon entropy, trajectory plots, Seurat clustering and the
interactive browser are all real; they were generalized from the lab's baseline scripts (kept unchanged in `analysis/baseline_scripts/`)
and are checked against them (`tests/testthat/test-baselines.R`).
* **DESeq2** compares every treated condition with its control at every time point (explicit contrasts, Wald test) on RAW counts (the Tux
  counts, or your table; it refuses normalized values) and needs replicates. Results per comparison, the DE genes (padj < 0.05, |log2 FC| ≥ 1
  by default), normalized counts, a volcano plot and a heatmap.
* **Entropy methods and trajectories need q-values**: from a Tux run (raw reads), calculated by StRym from a raw counts table or long-read counts with replicates (DESeq2, one model per comparison), or read from Cuffdiff `gene_exp.diff` files you give next to an FPKM table. They compare each *treated*
  condition with its *control* (found from names like control / mock / WT, or set in **Analysis settings**) at every time point, using
  Cuffdiff's q-values. Defaults are those of the original scripts: mean FPKM ≥ 2, |log2 FC| ≥ 1, q < 0.05, main set = changed at ≥ 2
  time points, first time point = baseline. Shannon entropy needs ≥ 2 time points, modified Shannon entropy ≥ 3.
* **PCA** works on Tux FPKM per replicate, or on a counts table (counts per million). *Samples to leave out* removes outliers from every step.
* **Trajectories** plot the genes found by an entropy step, or your own gene list (ids or names), treated vs control over time.
* **Single cell / single nucleus** data go through the Seurat step: 10x folders (one per sample, or ONE folder whose cell barcodes carry a
  sample label such as `SAMN15315800-AAAC…`, `Ctrl_AAAC…` or the cellranger-aggr `-1/-2`), a `.h5` file, or a Seurat `.rds`; PCs from the elbow rule,
  integration across samples (CCA, RPCA, Harmony or none), UMAP and/or t-SNE as ticked, cluster markers and treated-vs-control genes per cluster.
  A matrix holding several samples is split by the sample sheet (its *sample* column lists the labels; the design page drafts it from the barcodes;
  labels it does not list are left out, and without a sheet all cells are analyzed as one sample). Impossible counts (missing, infinite, negative, above
  10 million) are set to 0 and logged; cells below the UMI filter are dropped before any Seurat object is built; a sample with no cells left is skipped, and the
  log and `cells_per_sample.csv` say why. Big matrices are read with a low-memory reader (data.table). DESeq2 is not offered for single cell data: it reads a counts table.
* GENIST and RTP-STAR appear in the plan as hand-offs to their own tabs; the TUX module remains available under *Analysis modules*.
  Details of every step and of the rules that decide when it applies: `analysis/README.md`.
To add or change a step, edit its file in `analysis/`; the contract and the rules (`analysis/registry.R`) are described in `analysis/README.md`.
The new packages (ggplot2, Seurat, DESeq2, Harmony ...) are in `setup_wsl.sh`; on an existing install run `setup_modules.sh`
(it pins samtools 1.16.1 so the Tuxedo tools are not downgraded).

Paths can be typed (`D:\data\file.xlsx` or `/mnt/d/data/file.xlsx`) or browsed. Input data are only ever read.

## TUX (RNA-seq → DEGs)
* **Inputs**: a folder of FASTQ files (flat — groups guessed from file names, editable — or one sub-folder per condition),
  a GTF/GFF annotation, an uncompressed genome FASTA, and an output folder (any drive). Single- or paired-end
  (`A_/B_`, `_R1/_R2`, `_1/_2`); `.gz` is fine. Arabidopsis files: `reference/genes.gtf`, `reference/genome.fa`.
* Disk-friendly: HISAT2 is piped straight into `samtools sort` (BAM, no raw SAM), cleaned reads are deleted after
  each sample, free space is checked before starting. Spaces in paths are handled.
* **Output**: `<output>/diff_out/gene_exp.diff` and friends; the Results tab shows the table and a volcano plot.

## GENIST (time course → network)
* **Inputs**: genes file (list of genes), time-course file (≥ 3 time points), optional clustering data, TF list,
  gene-name file. All Excel/CSV, first row = header, gene ids in column 1.
* **Start from your Tux results**: choose a StRym project folder (or a Tux `diff_out` folder) and press *Build the genes file and the
  time course file*. StRym picks the genes that Tux finds significantly different between two time points (q ≤ 0.05, |log2 fold change| ≥ 1,
  FPKM ≥ 1, at most 500 genes by default) and writes the FPKM of every gene at each time point, in the time order of your design, to
  `<project>/06_regulatory_networks/GENIST/inputs`. The run time grows steeply with the gene count (about 2 min for 500 genes, over
  30 min for 1,500, on 2 threads). Nothing in the Tux folder is changed.
* **Parameters** (the defaults shown on the original TuxNet GENIST tab): time lapse 1 / 0 / both; regulation fold-change 1.1; regulation
  time percent 0.1; discretization levels 3; bottom fraction of regulations dropped 0.05. (The R function `genist_run()` itself keeps the
  defaults of the MATLAB code, 1.3 / 0.3 / 2 / 0.2.) The TF list and gene names are filled in at the start, as in the original.
* **Output**: `cytoscape_table.txt` (source, target, width, sign — import into Cytoscape) and `clusters.xlsx`;
  the Results tab draws the network (green = activates, red = inhibits) and lists the edges and clusters.

## RTP-STAR (expression → directed network)
* **Inputs**: genes file, expression file (plus a TF list or a `TF` column), optional gene names, clustering data,
  time course (edge direction), or your own clusters file.
* **Parameters** (TuxNet defaults): iterations 1, edge proportion 0.33 (iterations > 1), spatial/temporal clustering,
  connect cluster hubs, clustering seed, 10000 trees per forest.
* **Output**: `biograph.txt` (lines like `GENE_A activates GENE_B`; with several iterations also the proportion of runs
  in which each edge appeared), cluster files and per-iteration files.

## How the R ports relate to the MATLAB originals
The algorithms, defaults and tie-breaking follow the MATLAB code (GENIST2_0, GRN_preprocessing2_0, GRN_inference2_0,
RTPSTAR_MAIN, regression_tree_pipeline, run_regressiontree, clustering). They were checked against
the MATLAB formulas on small hand-worked cases and against the published GENIE3 package; **MATLAB itself was not
available, so results were not compared number-for-number with TuxNet's.** Deliberate differences:
* genes are matched across files by **name** (MATLAB matched rows by position in places and fails or mislabels genes
  when files are ordered differently or contain different genes);
* GENIE3 is re-implemented with `ranger` (same method as the published GENIE3) instead of TuxNet's GENIE3 MATLAB code;
* name files are read without assuming a header row (the original silently dropped their first row);
* cluster numbers may differ (R's `cutree` numbering) and k-means uses several random starts.

## Tests
`bash tests/run_tests.sh quick` runs the quick suite in the installed environment on WSL, Linux, or macOS.
On Windows, use `wsl -d Ubuntu -u root -e bash tests/run_tests.sh quick`.

## Files
| file | purpose |
|---|---|
| `app.R`, `about.R`, `mod_design.R`, `mod_genist.R`, `mod_rtpstar.R`, `gui_helpers.R`, `jobs.R`, `drives.R`, `www/` | the GUI (theme in `www/strym.css`) |
| `pipeline.R`, `analysis/`, `run_analysis_cli.R`, `make_demo_data.R` | Experimental Design engine, analysis step scripts, baseline scripts, runner, demo data |
| `tuxedo.R`, `run_cli.R` | TUX pipeline and its command-line runner |
| `genist.R`, `rtpstar.R`, `io_utils.R`, `run_genist_cli.R`, `run_rtpstar_cli.R` | GENIST / RTP-STAR and their runners |
| `start_gui.bat`, `start_gui.sh`, `stop_gui.sh` | launch / stop on Windows/WSL, Linux, or macOS |
| `setup_wsl.sh`, `setup_modules.sh` | platform-detecting one-time install / add packages to an existing install |
| `reference/` | Arabidopsis annotation, genome, TF list, gene names (from TuxNet) |
| `example_data/` | TuxNet's Arabidopsis example data for GENIST and RTP-STAR |
| `make_synthetic.R`, `make_longread_data.R`, `run_test.sh`, `verify.sh`, `all_tests.sh`, `tests/` | tests |

## Notes
* WSL runs as root here and uses `/root/mamba`. macOS uses `~/micromamba`; Apple Silicon runs the Intel-only Bioconda
  packages through Rosetta and keeps the legacy Tuxedo tools separate from current R/Bioconductor packages.
* Linux containers on Apple Silicon must use the `linux/amd64` platform. `setup_wsl.sh` rejects native Linux ARM early
  instead of allowing micromamba to fail after dependency resolution.
* The default WSL memory limit (half of RAM) may be too small for large genomes; raise it in `%UserProfile%\.wslconfig`.
* Credits: GENIST — M. Angels de Luis Balaguer; RTP-STAR — Natalie M. Clark (MIT license); GENIE3 — Huynh-Thu, Irrthum,
  Wehenkel & Geurts (PLoS ONE 2010); TuxNet — https://github.com/rspurney/TuxNet. The TuxNet repository has no top-level
  license file: confirm licensing and attribution with its authors before distributing StRym.

## Using your own Cuffdiff results with an expression table

Shannon entropy and modified Shannon entropy need q-values (Cuffdiff). If you did not start from raw reads, give StRym the expression table (FPKM / TPM or counts) **and** your Cuffdiff `gene_exp.diff` files under *Differential-expression results from another run*: one file, several (e.g. one per genotype and time point), a folder, or a `.zip`. StRym matches each comparison to your sample sheet by name (word order does not matter, and a bare "Control" is completed from the folder name); if the names do not match at all, it matches them by their expression values. Genes are matched by `test_id`, `gene_id` or gene name, whichever fits your table. A `gene_short_name` column in the table is used to name genes in the results. *What are the numbers in the table?* tells StRym whether the values are raw counts (turned into counts per million) or already normalized (used as given); "work it out" decides from whether there are decimals. If no time points are typed, they are read from the end of the sample names (`_T3`, `_5d`, `_day2`).

## Long reads (Oxford Nanopore, PacBio; e.g. Plasmidsaurus)

Choose **Raw reads (FASTQ)** and, under *Sequencing platform*, **Long reads**. StRym then aligns each file with minimap2 (splice-aware) and counts the reads per gene with featureCounts (long-read mode) instead of using the short-read Tuxedo route, and writes the same `counts.csv` the other steps read. Because there is no Cuffdiff run, the q-values for the entropy steps are calculated from the counts with DESeq2 (needs replicates). Install the two tools once with `setup_modules.sh` (they go into their own `longread` environment). If your provider already supplies a counts table (Plasmidsaurus does), you can simply use *A counts table*, which is quicker. A test set is in `example_data/protocol/exp6_longread_tiny` (made by `make_longread_data.R`).
