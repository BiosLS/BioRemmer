#!/bin/bash
# =============================================================================
#  BioRemmer v3.0 — setup_git.sh
#  Prepara el repositorio para subir a GitHub.
#  Corre desde la raíz de BioRemmer UNA SOLA VEZ.
# =============================================================================
set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()    { echo -e "${BLUE}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo ""
echo "============================================"
echo "  BioRemmer — Git setup"
echo "  Root: $ROOT"
echo "============================================"
echo ""

# ── 1. Limpiar Zone.Identifier (Windows) ─────────────────────────────────────
info "Removing Windows Zone.Identifier files..."
find . -name "*:Zone.Identifier" -delete 2>/dev/null && success "Cleaned"

# ── 2. Eliminar duplicado en raíz ────────────────────────────────────────────
if [ -f "install_bioremmer_v3.sh" ]; then
    rm install_bioremmer_v3.sh
    success "Removed duplicate install_bioremmer_v3.sh from root"
fi

# ── 3. Crear carpetas placeholder para gdrive (para que git las ignore) ───────
info "Creating .gitkeep placeholders for gdrive folders..."
mkdir -p databases/COG_LE
mkdir -p bin/vamphyre/genomes/genome_DB
mkdir -p test

for dir in databases/COG_LE bin/vamphyre/genomes/genome_DB test; do
    if [ ! -f "$dir/.gitkeep" ]; then
        touch "$dir/.gitkeep"
        success "  $dir/.gitkeep created"
    fi
done

# ── 4. Inicializar git si no existe ──────────────────────────────────────────
if [ ! -d ".git" ]; then
    info "Initializing git repository..."
    git init
    success "Git initialized"
else
    warn "Git already initialized — skipping git init"
fi

# ── 5. Agregar .gitignore ─────────────────────────────────────────────────────
info "Adding .gitignore..."
git add .gitignore
success ".gitignore staged"

# ── 6. Primer commit con todos los scripts ────────────────────────────────────
info "Staging all tracked files..."
git add \
    README.md \
    biorem_pipeline_v3.sh \
    config.sh \
    scripts/ \
    install/ \
    assets/Biorem_logo.png \
    bin/vamphyre/VPS/ \
    bin/vamphyre/bin/ \
    bin/vamphyre/scripts/ \
    bin/vamphyre/genomes/.gitkeep \
    databases/BioRemDB.csv \
    databases/pfam/ \
    databases/COG_LE/.gitkeep \
    test/README.md \
    test/.gitkeep \
    2>/dev/null || true

git status --short

echo ""
echo "============================================"
echo "  Ready for first commit!"
echo "  Run:"
echo ""
echo '    git commit -m "BioRemmer v3.0 — initial commit"'
echo ""
echo "  Then connect to GitHub:"
echo '    git remote add origin https://github.com/YOUR_USER/BioRemmer.git'
echo '    git branch -M main'
echo '    git push -u origin main'
echo "============================================"
echo ""
