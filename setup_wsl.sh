#!/usr/bin/env bash
# Installs micromamba + an isolated "tuxnet" env with the Tuxedo tools, R, Shiny and the packages for GENIST and RTP-STAR.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq && apt-get install -y -qq curl bzip2 ca-certificates >/dev/null
if [ ! -x /root/bin/micromamba ]; then
  mkdir -p /root/bin
  curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /root bin/micromamba
fi
export MAMBA_ROOT_PREFIX=/root/mamba
# samtools / htslib / zlib are pinned: the analysis packages below would otherwise make conda pick samtools 1.5.
/root/bin/micromamba create -y -n tuxnet -c conda-forge -c bioconda \
  "samtools=1.16.1" "htslib=1.16" "libzlib=1.2.13" "zlib=1.2.13" \
  hisat2 cufflinks ea-utils r-base r-shiny r-dt r-shinyfiles r-jsonlite r-processx r-testthat \
  r-readxl r-writexl r-igraph r-ranger r-visnetwork r-cluster r-htmlwidgets bioconductor-genie3 \
  r-ggplot2 r-viridis r-rcolorbrewer r-data.table r-patchwork r-tidyr r-seurat bioconductor-deseq2 r-harmony r-rtsne r-uwot r-ggrepel r-pheatmap r-hdf5r
echo SETUP_DONE
