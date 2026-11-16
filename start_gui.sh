#!/usr/bin/env bash
# Start the StRym GUI. STRYM_PORT defaults to 3838.
set -euo pipefail

STRYM_PORT="${STRYM_PORT:-3838}"
if [[ ! "$STRYM_PORT" =~ ^[0-9]+$ ]] || (( ${#STRYM_PORT} > 5 )); then
  echo "Invalid STRYM_PORT: $STRYM_PORT. Use a number from 1 to 65535." >&2
  exit 2
fi
STRYM_PORT=$((10#$STRYM_PORT))
if (( STRYM_PORT < 1 || STRYM_PORT > 65535 )); then
  echo "Invalid STRYM_PORT: $STRYM_PORT. Use a number from 1 to 65535." >&2
  exit 2
fi
export STRYM_PORT

if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
  [[ "$(uname -s)" == "Darwin" ]] && MAMBA_ROOT_PREFIX="${HOME}/micromamba" || MAMBA_ROOT_PREFIX="${HOME}/mamba"
fi
export MAMBA_ROOT_PREFIX
export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet/bin:${MAMBA_ROOT_PREFIX}/envs/tuxnet-tools/bin:$PATH"
cd "$(dirname "$0")"
bash ./stop_gui.sh >/dev/null 2>&1; sleep 1
echo "Starting StRym at http://localhost:${STRYM_PORT}"
exec Rscript -e 'app_dir_override <- getwd(); shiny::runApp("app.R", port = as.integer(Sys.getenv("STRYM_PORT")), host = "0.0.0.0", launch.browser = FALSE)'
