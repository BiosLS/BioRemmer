#!/usr/bin/env python3
"""
BioRemmer — fix_meg.py
Corrects the .meg format output by VFAT/VAMPhyRE so it is compatible
with ape::read.dist() and build_tree.R.

Fixes:
  1. Collapses multi-line !Description into one line
  2. Removes column index header row  [    1    2    3 ...]
  3. Strips [N] row/taxon index prefixes from all lines

Usage:
    python3 fix_meg.py <input.meg> [output.meg]
    (if output is omitted, overwrites input in place)
"""
import re
import sys
from pathlib import Path

def fix_meg(input_path: str, output_path: str = None):
    if output_path is None:
        output_path = input_path

    with open(input_path, 'r') as f:
        lines = f.readlines()

    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.rstrip()

        # Fix multi-line !Description → one line with ; at end
        if stripped == '!Description':
            desc_parts = []
            i += 1
            while i < len(lines) and lines[i].strip() != ';':
                desc_parts.append(lines[i].strip())
                i += 1
            out.append('!Description ' + ' '.join(desc_parts) + ';\n')
            i += 1  # skip the ';' line
            continue

        # Remove column index header line e.g. "[    1    2    3 ...]"
        if re.match(r'^\[\s*\d+\s+\d+', stripped):
            i += 1
            continue

        # Taxa lines: "[ 1] #Taxon_name" → "#Taxon_name"
        m_taxa = re.match(r'^\[\s*\d+\]\s+(#\S+)', stripped)
        if m_taxa:
            out.append(m_taxa.group(1) + '\n')
            i += 1
            continue

        # Data lines: "[ 2]    0.719040 ..." → "0.719040 ..."
        m_data = re.match(r'^\[\s*\d+\]\s*(.*)', stripped)
        if m_data:
            rest = m_data.group(1).strip()
            out.append((rest + '\n') if rest else '\n')
            i += 1
            continue

        out.append(line)
        i += 1

    with open(output_path, 'w') as f:
        f.writelines(out)

    taxa = [l for l in out if l.startswith('#') and not l.startswith('#mega')]
    data = [l for l in out if re.match(r'^\s*\d', l)]
    print(f"[fix_meg.py] Fixed: {Path(output_path).name} — "
          f"{len(taxa)} taxa, {len(data)} data rows")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python3 fix_meg.py <input.meg> [output.meg]")
    fix_meg(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
