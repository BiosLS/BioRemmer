#!/bin/bash
# =============================================================================
#  BioRemmer v1.0 — config.sh
#  Sourced by biorem_pipeline_v2.sh at runtime.
#  All paths are resolved relative to the BioRemmer root (BASE_DIR).
# =============================================================================

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Working directory (BioRemmer root) ───────────────────────────────────────
DIR="$BASE_DIR"

# ── Conda environments ────────────────────────────────────────────────────────
CORE_ENV="bioremmer_core"   # bioinformatics tools (steps 1–9)
R_ENV="bioremmer_r"         # R + report generation (step 6 tree + step 10)

# ── Databases ─────────────────────────────────────────────────────────────────
COG_database="${BASE_DIR}/databases/COG_LE/Cog"
pfam="${BASE_DIR}/databases/pfam"

# ── VAMPhyRE ──────────────────────────────────────────────────────────────────
vamphyre="${BASE_DIR}/bin/vamphyre/bin"
probe="${BASE_DIR}/bin/vamphyre/VPS/vps13.txt"

# genome: list of reference genome .fna paths (used by VH5cmdl and VFAT)
genome="${BASE_DIR}/bin/vamphyre/genomes/genome_list.txt"

# genomeDB: pre-computed VGF fingerprints for the reference DB (used by MergeVGF)
# NOTE: this is the fingerprint file (.txt), NOT the distance matrix (.meg)
genomeDB="${BASE_DIR}/bin/vamphyre/genomes/VGF_13mer_genomeDB.txt"

# ── Scripts ───────────────────────────────────────────────────────────────────
cog_counter="${BASE_DIR}/scripts"

# ── Output paths ──────────────────────────────────────────────────────────────
RESULTS="${BASE_DIR}/Results"
PLOTS="${RESULTS}/Plots"
LOGS="${RESULTS}/logs"
