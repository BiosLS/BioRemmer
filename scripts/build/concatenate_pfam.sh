#!/bin/bash
# Concatenates individual Pfam HMM profiles per plastic type.
# Run from the repo root: bash scripts/build/concatenate_pfam.sh
set -euo pipefail
PFAM_DIR="databases/pfam"
for plastic_dir in "$PFAM_DIR"/*/; do
    plastic=$(basename "$plastic_dir")
    out="${plastic_dir}conc_${plastic}.hmm"
    mapfile -t hmms < <(find "$plastic_dir" -maxdepth 1 -name "PF*.hmm" | sort)
    if [ ${#hmms[@]} -eq 0 ]; then
        echo "  No individual HMMs for $plastic — skipping"
        continue
    fi
    cat "${hmms[@]}" > "$out"
    echo "  ✓ ${plastic}: $(grep -c '^NAME' "$out") profiles → $out"
done
