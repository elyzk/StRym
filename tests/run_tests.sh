#!/usr/bin/env bash
# Run the tests inside WSL and show just the summary (the real tools print a lot).
#   bash tests/run_tests.sh            all tests (30-60 minutes: the real Tuxedo tools and the single cell clustering runs are slow)
#   bash tests/run_tests.sh quick      everything EXCEPT the slow groups (resume, quantify, seurat, seurat-multi): about 6 minutes
#   bash tests/run_tests.sh quantify   only tests/testthat/test-quantify.R
if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
  [[ "$(uname -s)" == "Darwin" ]] && MAMBA_ROOT_PREFIX="${HOME}/micromamba" || MAMBA_ROOT_PREFIX="${HOME}/mamba"
fi
export MAMBA_ROOT_PREFIX
export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet/bin:${MAMBA_ROOT_PREFIX}/envs/tuxnet-tools/bin:$PATH"
cd "$(dirname "$0")/.." || exit 1
log=$(mktemp)
SUMMARY="cat('\nPASSED:', sum(r\$passed), ' FAILED:', sum(r\$failed), ' ERRORS:', sum(r\$error), ' SKIPPED:', sum(r\$skipped), '\n')"
if [ "$1" = "quick" ]; then
  Rscript -e "library(testthat); setwd('tests/testthat'); r <- as.data.frame(test_dir('.', filter = 'baselines|batch3|deseq2|drives|edgecases|genist|pipeline|review-fixes|rtpstar|rtpstar-tux|singlecell-light|spaces-runner|steps|v9|v10|v11', reporter = 'summary', stop_on_failure = FALSE)); $SUMMARY" > "$log" 2>&1
elif [ -n "$1" ]; then
  Rscript -e "library(testthat); setwd('tests/testthat'); r <- as.data.frame(test_file('test-$1.R', reporter = 'summary', stop_on_failure = FALSE)); $SUMMARY" > "$log" 2>&1
else
  Rscript tests/run_tests.R > "$log" 2>&1
fi
awk '/^══ Failed/{show=1} show{print} /^(baselines|batch3|deseq2|drives|edgecases|genist|genist-tux|pipeline|quantify|resume|review-fixes|rtpstar|rtpstar-tux|seurat|seurat-multi|singlecell-light|spaces-runner|steps|v9|v10|v11):/{print substr($0, 1, 90)} /^PASSED:/{print}' "$log" | grep -v "built under R version"
rm -f "$log"
