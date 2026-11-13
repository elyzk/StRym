#!/usr/bin/env bash
# Run an analysis from the command line (no browser): bash run_analysis.sh /path/to/config.json
# config.json = {"design": { ...contents of 00_design/design.json... }, "steps": ["deseq2", "pca", "qc_overview"]}
# Prints  DONE: <project>/index.html  (exit 0)  or  ERROR: ...  (exit 1).
if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
  [[ "$(uname -s)" == "Darwin" ]] && MAMBA_ROOT_PREFIX="${HOME}/micromamba" || MAMBA_ROOT_PREFIX="${HOME}/mamba"
fi
export MAMBA_ROOT_PREFIX
export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet/bin:${MAMBA_ROOT_PREFIX}/envs/tuxnet-tools/bin:$PATH"
cd "$(dirname "$0")" || exit 1
exec Rscript run_analysis_cli.R "$1"