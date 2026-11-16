#!/usr/bin/env bash
# Installs micromamba + isolated environments with StRym's R packages and command-line tools.
# Linux/WSL uses one environment. macOS uses Rosetta for packages that Bioconda does not publish for Apple Silicon.
set -euo pipefail

OS="$(uname -s)"
ARCH="$(uname -m)"
BIN_DIR="${HOME}/.local/bin"
if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
  if [[ "$OS" == "Linux" ]]; then
    MAMBA_ROOT_PREFIX="${HOME}/mamba"
  else
    MAMBA_ROOT_PREFIX="${HOME}/micromamba"
  fi
fi
export MAMBA_ROOT_PREFIX

case "$OS" in
  Linux)  PLATFORM="linux" ;;
  Darwin) PLATFORM="osx" ;;
  *)      echo "Unsupported OS: $OS" >&2; exit 1 ;;
esac

case "$ARCH" in
  x86_64|amd64)  ARCH_TAG="64" ;;
  arm64|aarch64) ARCH_TAG="aarch64" ;;
  *)             echo "Unsupported architecture: $ARCH" >&2; exit 1 ;;
esac

# The following packages aren't provided for ARM systems by Bioconda. Came across this during testing. This problem doesn't occur on Apple silicon because the packages run through Rosetta 2
if [[ "$OS" == "Linux" && "$ARCH_TAG" == "aarch64" ]]; then
  echo "Unsupported Linux architecture: $ARCH." >&2
  echo "Cufflinks 2.2.1, ea-utils, and the pinned SAMtools/HTSlib builds require linux/amd64." >&2
  echo "For Docker on an ARM host, recreate the container with: --platform linux/amd64" >&2
  exit 1
fi

if [[ "$PLATFORM" == "osx" && "$ARCH_TAG" == "aarch64" ]]; then
  ARCH_TAG="arm64"
fi

mkdir -p "$BIN_DIR"
export PATH="$BIN_DIR:$PATH"
export DEBIAN_FRONTEND=noninteractive


# This was needed to make this setup script work on the docker linux/amd64 container I was using for testing. It should not be required for a native Linux x86 system and could be removed.
NEED_BOOTSTRAP_DEPS=false
if ! command -v curl >/dev/null 2>&1 || ! command -v bzip2 >/dev/null 2>&1; then
  NEED_BOOTSTRAP_DEPS=true
fi
if [[ "$OS" == "Linux" ]] && ! command -v bsdtar >/dev/null 2>&1; then
  NEED_BOOTSTRAP_DEPS=true
fi
if [[ "$NEED_BOOTSTRAP_DEPS" == true ]]; then
  if command -v apt-get >/dev/null 2>&1; then
    SUDO=""
    if [[ "$(id -u)" -ne 0 ]] && command -v sudo >/dev/null 2>&1; then
      SUDO="sudo"
    fi
    PACKAGES=(curl bzip2 ca-certificates)
    [[ "$OS" == "Linux" ]] && PACKAGES+=(libarchive-tools)
      echo "Need curl and bzip2 to be installed. Install them manually or enable sudo access to install them automatically."
    $SUDO apt-get update -qq
    $SUDO apt-get install -y -qq "${PACKAGES[@]}" >/dev/null
  elif ! command -v curl >/dev/null 2>&1 || ! command -v bzip2 >/dev/null 2>&1; then
    echo "Error: curl and bzip2 are required. Install them first." >&2
    exit 1
  fi
fi

RELEASE_URL="https://micro.mamba.pm/api/micromamba/${PLATFORM}-${ARCH_TAG}/latest"
if ! command -v micromamba >/dev/null 2>&1; then
  TMP_DIR="$(mktemp -d)"
  trap 'rm -rf "${TMP_DIR:-}" "${EA_TMP_DIR:-}"' EXIT
  MAMBA_ARCHIVE="${TMP_DIR}/micromamba.tar.bz2"
  curl -fLsS --retry 5 --retry-delay 2 --retry-max-time 120 "$RELEASE_URL" -o "$MAMBA_ARCHIVE"
  if command -v bsdtar >/dev/null 2>&1; then
    bsdtar -xjf "$MAMBA_ARCHIVE" -C "$TMP_DIR" bin/micromamba
  else
    tar -xjf "$MAMBA_ARCHIVE" -C "$TMP_DIR" bin/micromamba
  fi
  install -m 0755 "$TMP_DIR/bin/micromamba" "$BIN_DIR/micromamba"
fi
MAMBA="$(command -v micromamba)"
CHANNELS=(--override-channels -c conda-forge -c bioconda)
R_PACKAGES=(
  r-base r-shiny r-dt r-shinyfiles r-jsonlite r-processx r-testthat
  r-readxl r-writexl r-igraph r-ranger r-visnetwork r-cluster r-htmlwidgets bioconductor-genie3
  r-ggplot2 r-viridis r-rcolorbrewer r-data.table r-patchwork r-tidyr r-seurat bioconductor-deseq2
  r-harmony r-rtsne r-uwot r-ggrepel r-pheatmap r-hdf5r
)

if [[ "$OS" == "Darwin" ]]; then
  # Check Rosetta 2 installation if on Apple Silicon
  if [[ "$ARCH" == "arm64" ]] && ! /usr/bin/arch -x86_64 /usr/bin/true >/dev/null 2>&1; then
    echo "Rosetta 2 is required for Cufflinks and GENIE3." >&2
    echo "Install it with: softwareupdate --install-rosetta" >&2
    exit 1
  fi
  # Check for Xcode
  for tool in g++ make shasum; do
    if ! command -v "$tool" >/dev/null 2>&1; then
      echo "Apple command-line developer tools are required. Run: xcode-select --install" >&2
      exit 1
    fi
  done

  # Create the conda environment
  "$MAMBA" create -y -n tuxnet --platform osx-64 "${CHANNELS[@]}" "${R_PACKAGES[@]}"

  # Need to separate packages with different dependency versions
  echo "Creating the macOS Tuxedo tool environment..."
  "$MAMBA" create -y -n tuxnet-tools --platform osx-64 "${CHANNELS[@]}" \
    "samtools=1.16.1" "htslib=1.16" "libzlib=1.2.13" "zlib=1.2.13" \
    "hisat2=2.2.3" "cufflinks=2.2.1"

  TOOLS_BIN="${MAMBA_ROOT_PREFIX}/envs/tuxnet-tools/bin"

  # Need to install fastq-mcf manually because Bioconda doesn't publish ea-utils for macOS
  if [[ ! -x "${TOOLS_BIN}/fastq-mcf" ]]; then
    echo "Installing ea-utils and setting up fastq-mcf..."
    EA_VERSION="1.04.807"
    EA_SHA256="aa09d25e6aa7ae71d2ce4198a98e58d563f151f8ff248e4602fa437f12b8d05f"
    EA_TMP_DIR="$(mktemp -d)"
    trap 'rm -rf "${TMP_DIR:-}" "${EA_TMP_DIR:-}"' EXIT
    curl -Ls "https://github.com/ExpressionAnalysis/ea-utils/archive/refs/tags/${EA_VERSION}.tar.gz" \
      -o "${EA_TMP_DIR}/ea-utils.tar.gz"
    printf '%s  %s\n' "$EA_SHA256" "${EA_TMP_DIR}/ea-utils.tar.gz" | shasum -a 256 -c -
    tar -xzf "${EA_TMP_DIR}/ea-utils.tar.gz" -C "$EA_TMP_DIR"
    EA_SRC="${EA_TMP_DIR}/ea-utils-${EA_VERSION}/clipper"
    (
      cd "${EA_SRC}/sparsehash-2.0.2"
      ./configure >/dev/null
      make src/sparsehash/internal/sparseconfig.h >/dev/null
    )
    mkdir "${EA_SRC}/sparsehash"
    cp -R "${EA_SRC}/sparsehash-2.0.2/src/sparsehash/." "${EA_SRC}/sparsehash/"
    g++ -O3 -I"$EA_SRC" "${EA_SRC}/fastq-mcf.cpp" "${EA_SRC}/fastq-lib.cpp" \
      -o "${EA_TMP_DIR}/fastq-mcf"
    install -m 0755 "${EA_TMP_DIR}/fastq-mcf" "${TOOLS_BIN}/fastq-mcf"
  fi

  export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet/bin:${TOOLS_BIN}:$PATH"
else
  "$MAMBA" create -y -n tuxnet "${CHANNELS[@]}" \
    "samtools=1.16.1" "htslib=1.16" "libzlib=1.2.13" "zlib=1.2.13" \
    hisat2 cufflinks ea-utils "${R_PACKAGES[@]}"
  export PATH="${MAMBA_ROOT_PREFIX}/envs/tuxnet/bin:$PATH"
fi

Rscript -e "library(shiny); library(Seurat); library(DESeq2); cat('R OK\n')"
