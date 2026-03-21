#!/usr/bin/env python3
"""
BioRemmer — parse_cog_rpsblast.py
Parses rpsblast text output and returns best COG hit per protein.
Used by plastic_heatmap.py and the R report via CSV output.

Usage:
    python3 parse_cog_rpsblast.py --blast <COGs_results.faa> \
                                   --output <cog_hits.csv>
"""
import re
import sys
import argparse
import csv
from pathlib import Path

COG_CATEGORIES = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking and secretion",
    "O": "Posttranslational modification, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis",
    "R": "General function prediction only",
    "S": "Function unknown",
}

# COG categories linked to plastic degradation
PLASTIC_RELATED_COGS = {"I", "Q", "C", "E", "G", "P"}

def parse_rpsblast(blast_file: str) -> dict:
    """
    Parse rpsblast text output.
    Returns dict: protein_id -> {cog_id, cog_name, cog_cat, evalue, description}
    """
    results = {}
    current_query = None
    best_hit = None
    best_evalue = float('inf')

    cdd_pattern   = re.compile(r'CDD:\d+\s+(COG\d+),\s*([^,]+),\s*(.+?)\s+[\d.e+-]+\s+([\d.e+-]+)\s*$')
    cdd_pattern2  = re.compile(r'(COG\d+),\s*[^,]*,\s*(.+?)\[(.+?)\]')
    query_pattern = re.compile(r'^Query=\s+(\S+)')
    evalue_pattern = re.compile(r'>\s*CDD:\d+\s+(COG\d+),\s*[^,]+,\s*(.+)')
    score_pattern  = re.compile(r'Expect\s*=\s*([\d.e+-]+)')

    with open(blast_file) as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]

        # New query
        m = query_pattern.match(line)
        if m:
            if current_query and best_hit:
                results[current_query] = best_hit
            current_query = m.group(1)
            best_hit = None
            best_evalue = float('inf')
            i += 1
            continue

        # Hit line in alignment summary table
        # e.g.: "CDD:223489 COG0412, ... 64.7  1e-12"
        m = cdd_pattern.match(line.strip())
        if m and current_query:
            cog_id   = m.group(1)
            evalue_s = m.group(4)
            try:
                evalue = float(evalue_s)
            except ValueError:
                evalue = float('inf')
            if evalue < best_evalue:
                best_evalue = evalue
                best_hit = {
                    'cog_id':   cog_id,
                    'evalue':   evalue,
                    'cog_name': m.group(2).strip(),
                    'cog_desc': m.group(3).strip(),
                    'cog_cat':  '',
                }
            i += 1
            continue

        # Full alignment header: ">CDD:XXXXX COG####, name, desc [category]"
        if line.startswith('>CDD:') and current_query:
            m2 = evalue_pattern.match(line.strip())
            if m2:
                cog_id   = m2.group(1)
                full_desc = m2.group(2).strip()
                # Extract category from brackets if present
                cat_match = re.search(r'\[(.+?)\]', full_desc)
                cat = cat_match.group(1) if cat_match else ''
                # Look ahead for Expect value
                for j in range(i+1, min(i+10, len(lines))):
                    em = score_pattern.search(lines[j])
                    if em:
                        try:
                            evalue = float(em.group(1))
                        except ValueError:
                            evalue = float('inf')
                        if best_hit and best_hit['cog_id'] == cog_id and not best_hit['cog_cat']:
                            best_hit['cog_cat'] = cat
                        break
            i += 1
            continue

        i += 1

    # Save last query
    if current_query and best_hit:
        results[current_query] = best_hit

    # Infer COG single-letter category from known mappings
    for pid, hit in results.items():
        if not hit['cog_cat']:
            # Try to infer from description
            desc_lower = hit.get('cog_desc', '').lower() + hit.get('cog_name', '').lower()
            if any(k in desc_lower for k in ['lipid', 'esterase', 'lipase', 'fatty']):
                hit['cog_cat'] = 'I'
            elif any(k in desc_lower for k in ['secondary', 'oxygenase', 'monooxygenase']):
                hit['cog_cat'] = 'Q'
            elif any(k in desc_lower for k in ['energy', 'oxidoreductase', 'dehydrogenase']):
                hit['cog_cat'] = 'C'

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast',  required=True, help='rpsblast output file')
    parser.add_argument('--output', required=True, help='Output CSV file')
    args = parser.parse_args()

    print(f"[parse_cog_rpsblast.py] Parsing: {args.blast}")
    results = parse_rpsblast(args.blast)
    print(f"[parse_cog_rpsblast.py] Proteins with COG hit: {len(results):,}")

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            'protein', 'cog_id', 'cog_name', 'cog_desc', 'cog_cat', 'evalue',
            'plastic_related'
        ])
        writer.writeheader()
        for pid, hit in sorted(results.items()):
            writer.writerow({
                'protein':         pid,
                'cog_id':          hit['cog_id'],
                'cog_name':        hit.get('cog_name', ''),
                'cog_desc':        hit.get('cog_desc', ''),
                'cog_cat':         hit.get('cog_cat', ''),
                'evalue':          hit['evalue'],
                'plastic_related': hit.get('cog_cat', '') in PLASTIC_RELATED_COGS,
            })
    print(f"[parse_cog_rpsblast.py] Written: {args.output}")


if __name__ == '__main__':
    main()
