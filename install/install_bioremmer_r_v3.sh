#!/bin/bash
# =============================================================================
#  BioRemmer v3.0 — install_bioremmer_r_v3.sh
#  Creates or recreates the bioremmer_r conda environment.
#  Includes all R packages needed for the v3 report.
#
#  Uso: bash install/install_bioremmer_r_v3.sh
#  Ejecutar desde la raíz de BioRemmer.
# =============================================================================
set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; NC='\033[0m'
info()    { echo -e "${BLUE}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIOREM_DIR="$(dirname "$SCRIPT_DIR")"
ENV_NAME="bioremmer_r"
YML_FILE="${SCRIPT_DIR}/environment_bioremmer_r_v3.yml"
LOG="${BIOREM_DIR}/Results/logs/install_bioremmer_r_v3.log"
mkdir -p "$(dirname "$LOG")"
exec > >(tee -a "$LOG") 2>&1

info "BioRemmer root: $BIOREM_DIR"
info "YML file: $YML_FILE"

command -v conda >/dev/null 2>&1 || error "conda not found. Install Miniconda first."
source "${HOME}/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true

# Install mamba if not available
if ! conda run -n base mamba --version >/dev/null 2>&1; then
  info "Installing mamba in base..."
  conda install -n base -c conda-forge mamba -y
fi

# Remove existing env if present
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  warn "Removing existing env: $ENV_NAME"
  conda env remove -n "$ENV_NAME" -y || true
fi

# Create environment from yml
info "Creating R environment: $ENV_NAME"
[ -f "$YML_FILE" ] || error "YML not found: $YML_FILE"
mamba env create -f "$YML_FILE" \
  --override-channels -c conda-forge -c bioconda --strict-channel-priority
success "R environment created"

# Verify key tools
info "Verifying R and Rscript..."
for cmd in R Rscript; do
  if conda run -n "$ENV_NAME" bash -lc "command -v $cmd" >/dev/null 2>&1; then
    success "  $cmd found"
  else
    warn "  $cmd NOT found"
  fi
done

# Verify all required R packages
info "Verifying R packages..."
conda run -n "$ENV_NAME" Rscript -q -e "
pkgs <- c(
  'rmarkdown', 'ggplot2', 'data.table', 'dplyr', 'tidyr',
  'viridis', 'huxtable', 'knitr', 'kableExtra', 'RColorBrewer',
  'jsonlite', 'ape', 'phangorn', 'ggtree', 'treeio'
)
ok  <- sapply(pkgs, requireNamespace, quietly = TRUE)
cat('Installed in R env:\n')
print(ok)
cat('\nTotal:', sum(ok), '/', length(ok), '\n')
if (any(!ok)) cat('MISSING:', paste(names(ok)[!ok], collapse = ', '), '\n')
"

echo ""
echo "============================================"
echo "  R environment ready — BioRemmer v3.0"
echo "  Log: $LOG"
echo "============================================"
echo ""
echo "  Quick checks:"
echo "    conda activate bioremmer_r"
echo "    R --version"
echo "    Rscript scripts/report.R /path/to/BioRemmer"
echo ""
