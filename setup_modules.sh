#!/usr/bin/env bash
# Adds GENIST/RTP-STAR and analysis packages to an existing StRym installation.
set -euo pipefail

OS="$(uname -s)"
if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
  if [[ "$OS" == "Darwin" ]]; then
    MAMBA_ROOT_PREFIX="${HOME}/micromamba"
  else
    MAMBA_ROOT_PREFIX="${HOME}/mamba"
  fi
fi
export MAMBA_ROOT_PREFIX
export PATH="${HOME}/.local/bin:$PATH"
MAMBA="$(command -v micromamba)"
CHANNELS=(--override-channels -c conda-forge -c bioconda)
PLATFORM_ARGS=()
TOOL_PINS=("samtools=1.16.1" "htslib=1.16" "libzlib=1.2.13" "zlib=1.2.13")
if [[ "$OS" == "Darwin" ]]; then
  PLATFORM_ARGS=(--platform osx-64)
  TOOL_PINS=()
fi

"$MAMBA" install -y -n tuxnet "${PLATFORM_ARGS[@]}" "${CHANNELS[@]}" "${TOOL_PINS[@]}" \
  r-readxl r-writexl r-igraph r-ranger r-visnetwork r-cluster r-htmlwidgets bioconductor-genie3 \
  r-ggplot2 r-viridis r-rcolorbrewer r-data.table r-patchwork r-tidyr r-seurat bioconductor-deseq2 \
  r-harmony r-rtsne r-uwot r-ggrepel r-pheatmap r-hdf5r

export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet/bin:$PATH"
if [[ "$OS" == "Darwin" ]]; then
  export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet-tools/bin:$PATH"
fi
samtools --version | /usr/bin/head -1
Rscript -e 'for(p in c("readxl","writexl","igraph","ranger","visNetwork","cluster","GENIE3","htmlwidgets","ggplot2","Seurat","DESeq2","harmony","Rtsne","uwot","patchwork","viridis","data.table")) cat(p, requireNamespace(p,quietly=TRUE), as.character(tryCatch(packageVersion(p),error=function(e)NA)), "\n")'
echo MODULES_SETUP_DONE

# Long-read tools stay in their own environment so the pinned Tuxedo toolchain is never changed.
LONG_BIN="${MAMBA_ROOT_PREFIX}/envs/longread/bin"
if [[ ! -x "${LONG_BIN}/minimap2" ]]; then
  "$MAMBA" create -y -n longread "${PLATFORM_ARGS[@]}" "${CHANNELS[@]}" minimap2 subread samtools
fi
"${LONG_BIN}/minimap2" --version
"${LONG_BIN}/featureCounts" -v 2>&1 | /usr/bin/head -2
echo LONGREAD_SETUP_DONE