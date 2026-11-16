# StRym tests

## The test suite (testthat)
Run inside WSL (Windows) or on Linux, from the StRym folder:

    bash tests/run_tests.sh            # everything (20 to 30 minutes; prints a short summary and "PASSED: n  FAILED: 0")
    bash tests/run_tests.sh v10        # one file, here tests/testthat/test-v10.R
    Rscript tests/run_tests.R          # everything, with the full testthat report

What the files cover:

| file | covers |
|---|---|
| `test-genist*.R`, `test-rtpstar*.R` | the GENIST and RTP-STAR algorithms, their checks, Tux-result inputs |
| `test-pipeline.R`, `test-steps.R`, `test-derived.R`, `test-diffinput.R`, `test-robust.R`, `test-edgecases.R` | the Experimental Design engine, every analysis step, messy input |
| `test-quantify.R`, `test-reads-robust.R`, `test-longreads.R`, `test-resume.R` | the Tux route from FASTQ (needs the Tux tools), long reads, resuming a run |
| `test-seurat*.R`, `test-singlecell-light.R` | single cell clustering |
| `test-drives.R`, `test-v10.R`, `test-v11.R` | the drives code on Windows/WSL, Mac and Linux; the Browse window (listing, places), the drive list that never waits, releasing a busy drive (uses a test file system when run as root in WSL) |
| `test-v9.R`, `test-rtpstar-tux.R` | PCA levels, entropy heatmaps, network names and pictures, results filter, RTP-STAR inputs from Tux |
| `test-baselines.R` | the analysis steps against the lab's original scripts (`analysis/baseline_scripts/`) |

A few tests skip themselves when a tool is missing (for example the zip program); the summary line counts them as SKIPPED.

## Test cases for people (made-up data)
`Rscript tests/make_test_sets.R <output folder> <StRym folder>` builds five folders, each with a `README.txt` that says what to do and what a correct run shows:

1. `V10_1_drive_reads_tux_to_rtpstar`: reads and a GFF3 annotation on a drive; Tux; RTP-STAR from the Tux results; unplug test (needs a USB stick or external drive)
2. `V10_2_drive_timecourse_clutter`: a time course on a drive; an output folder that already holds other files; sample sheet follows the page; PCA levels; heatmaps; GENIST with and without a TF list (needs a drive)
3. `V10_3_heatmap_pca_edge_cases`: tiny, 20,000-gene, many-condition and one-responder tables
4. `V10_4_cuffdiff_single_file_modes`: one `gene_exp.diff` file; the GENIST and RTP-STAR tabs started from Tux results; your own clusters file
5. `V10_5_network_tabs_browser`: network file names and pictures; the TF list as advice only; extra result folders in the interactive browser

`make_synthetic.R` (reads) and `make_longread_data.R` (long reads) in the StRym folder make smaller sets for the TUX route.