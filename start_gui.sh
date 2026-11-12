#!/usr/bin/env bash
# Start the StRym GUI; open http://localhost:3838 in a Windows browser.
export PATH=/root/mamba/envs/tuxnet/bin:$PATH
cd "$(dirname "$0")"
bash ./stop_gui.sh >/dev/null 2>&1; sleep 1      # never fail on a stale server holding the port
exec Rscript -e 'app_dir_override <- getwd(); shiny::runApp("app.R", port = 3838, host = "0.0.0.0", launch.browser = FALSE)'
