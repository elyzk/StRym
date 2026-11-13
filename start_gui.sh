#!/usr/bin/env bash
# Start the StRym GUI at http://localhost:3838.
if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
  [[ "$(uname -s)" == "Darwin" ]] && MAMBA_ROOT_PREFIX="${HOME}/micromamba" || MAMBA_ROOT_PREFIX="${HOME}/mamba"
fi
export MAMBA_ROOT_PREFIX
export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet/bin:${MAMBA_ROOT_PREFIX}/envs/tuxnet-tools/bin:$PATH"
cd "$(dirname "$0")"
bash ./stop_gui.sh >/dev/null 2>&1; sleep 1
exec Rscript -e 'app_dir_override <- getwd(); shiny::runApp("app.R", port = 3838, host = "0.0.0.0", launch.browser = FALSE)'
