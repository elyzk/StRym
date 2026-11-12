#!/usr/bin/env bash
# Run an analysis from the command line (no browser): bash run_analysis.sh /path/to/config.json
# config.json = {"design": { ...contents of 00_design/design.json... }, "steps": ["deseq2", "pca", "qc_overview"]}
# Prints  DONE: <project>/index.html  (exit 0)  or  ERROR: ...  (exit 1).
export PATH=/root/mamba/envs/tuxnet/bin:$PATH
cd "$(dirname "$0")" || exit 1
exec Rscript run_analysis_cli.R "$1"