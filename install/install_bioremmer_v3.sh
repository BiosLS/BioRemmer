#!/bin/bash
# =============================================================================
#  BioRemmer v3.0 — install_bioremmer_v3.sh
#  Instala fastp y dependencias Python (matplotlib, pandas, numpy)
#  en el entorno bioremmer_core existente.
#
#  Uso: bash install/install_bioremmer_v3.sh
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
# install/ is inside BioRemmer root, so go up one level
BIOREM_DIR="$(dirname "$SCRIPT_DIR")"
LOG="${BIOREM_DIR}/Results/logs/install_bioremmer_v3.log"
mkdir -p "$(dirname "$LOG")"
exec > >(tee -a "$LOG") 2>&1

echo ""
echo "============================================"
echo "  BioRemmer v3.0 — Installation"
echo "  $(date)"
echo "  Root: $BIOREM_DIR"
echo "============================================"
echo ""

command -v conda >/dev/null 2>&1 || error "conda not found. Install Miniconda first."
source "${HOME}/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true

# ── 1. fastp → bioremmer_core ─────────────────────────────────────────────────
info "Installing fastp in bioremmer_core..."
conda install -n bioremmer_core -c bioconda -c conda-forge fastp -y
conda run -n bioremmer_core fastp --version 2>&1 | head -1 && \
  success "fastp installed" || error "fastp failed"

# ── 2. Python dependencies for heatmap → bioremmer_core ──────────────────────
info "Installing matplotlib, pandas, numpy in bioremmer_core..."
conda install -n bioremmer_core -c conda-forge matplotlib pandas numpy -y
conda run -n bioremmer_core python3 -c \
  "import matplotlib, pandas, numpy; print('OK — matplotlib', matplotlib.__version__)" && \
  success "Python heatmap dependencies installed" || \
  error "Python dependencies failed"

# ── 3. R packages → bioremmer_r ───────────────────────────────────────────────
info "Installing new R packages in bioremmer_r..."
conda install -n bioremmer_r -c conda-forge \
  r-jsonlite r-kableextra r-rcolorbrewer -y
conda run -n bioremmer_r Rscript -e \
  "pkgs <- c('jsonlite','kableExtra','RColorBrewer')
   ok <- sapply(pkgs, requireNamespace, quietly=TRUE)
   cat('R packages:\n'); print(ok)" && \
  success "R packages installed" || \
  warn "Some R packages may need manual install"

# ── 4. Copy scripts if not already in place ───────────────────────────────────
info "Checking scripts..."
for script in plastic_heatmap.py fix_meg.py parse_cog_rpsblast.py; do
  if [ -f "${BIOREM_DIR}/scripts/${script}" ]; then
    success "  ${script} already in scripts/"
  else
    warn "  ${script} not found in scripts/ — copy it manually"
  fi
done

# ── 5. Final verification ─────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Final verification"
echo "============================================"
conda run -n bioremmer_core bash -c "
  echo -n '  fastp:      '; fastp --version 2>&1 | head -1
  echo -n '  python:     '; python3 -c 'import matplotlib, pandas, numpy; print(\"OK\")'
  echo -n '  matplotlib: '; python3 -c 'import matplotlib; print(matplotlib.__version__)'
"
conda run -n bioremmer_r Rscript -e "
  pkgs <- c('rmarkdown','ggplot2','viridis','huxtable','ape','ggtree',
            'knitr','kableExtra','RColorBrewer','jsonlite','phangorn')
  ok <- sapply(pkgs, requireNamespace, quietly=TRUE)
  cat('R packages OK:', sum(ok), '/', length(ok), '\n')
  if (any(!ok)) cat('Missing:', paste(names(ok)[!ok], collapse=', '), '\n')
"

echo ""
echo "============================================"
echo "  Installation complete — BioRemmer v3.0"
echo "  Log: $LOG"
echo "============================================"
echo ""
echo "  Test the heatmap:"
echo "    conda run -n bioremmer_core python3 scripts/plastic_heatmap.py \\"
echo "        --results Results --output Results/Plots"
echo ""
