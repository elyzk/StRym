#!/usr/bin/env bash
# Adds the packages for the GENIST / RTP-STAR modules and for the analysis steps (PCA, entropy, trajectories, Seurat, DESeq2)
# to the existing "tuxnet" env. The Tuxedo tools are pinned to the versions StRym was built with (samtools 1.16.1 ...): without the
# pins, conda would downgrade samtools to 1.5 to fit the new packages.
set -euo pipefail
export MAMBA_ROOT_PREFIX=/root/mamba
/root/bin/micromamba install -y -n tuxnet -c conda-forge -c bioconda \
  "samtools=1.16.1" "htslib=1.16" "libzlib=1.2.13" "zlib=1.2.13" \
  r-readxl r-writexl r-igraph r-ranger r-visnetwork r-cluster r-htmlwidgets bioconductor-genie3 \
  r-ggplot2 r-viridis r-rcolorbrewer r-data.table r-patchwork r-tidyr r-seurat bioconductor-deseq2 r-harmony r-rtsne r-uwot r-ggrepel r-pheatmap r-hdf5r
export PATH=/root/mamba/envs/tuxnet/bin:$PATH
samtools --version | head -1
Rscript -e 'for(p in c("readxl","writexl","igraph","ranger","visNetwork","cluster","GENIE3","htmlwidgets","ggplot2","Seurat","DESeq2","harmony","Rtsne","uwot","patchwork","viridis","data.table")) cat(p, requireNamespace(p,quietly=TRUE), as.character(tryCatch(packageVersion(p),error=function(e)NA)), "\n")'
echo MODULES_SETUP_DONE

# ---- Long reads (Oxford Nanopore, PacBio; e.g. Plasmidsaurus): minimap2 + featureCounts in their OWN environment, so the Tuxedo
# environment above is never changed. Skip this block if you only have short reads or a counts table.
if [ ! -x /root/mamba/envs/longread/bin/minimap2 ]; then
  /root/bin/micromamba create -y -n longread -c conda-forge -c bioconda minimap2 subread samtools
fi
/root/mamba/envs/longread/bin/minimap2 --version
/root/mamba/envs/longread/bin/featureCounts -v 2>&1 | head -2
echo LONGREAD_SETUP_DONE