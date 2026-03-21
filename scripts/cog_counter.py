#!/usr/bin/env python3
"""
BioRemmer v1.0 — cog_counter.py
Parses RPS-BLAST tabular output and counts COG functional categories.
Replaces COGcounter_v2.pl.

Usage:
    python cog_counter.py --blast <blast_results.faa> \
                          --cognames <cognames2003-2014.tab> \
                          --output <cog_frequencies.csv>

    # Or using a file list (same interface as the Perl version):
    python cog_counter.py --list <cogs_list.txt> \
                          --cognames <cognames2003-2014.tab> \
                          --output <cog_frequencies.csv>
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# ── COG categories (same order as original Perl script) ──────────────────────
COG_CATEGORIES = [
    "J", "A", "K", "L", "B", "D", "Y", "V", "T", "M",
    "N", "Z", "W", "U", "O", "X", "C", "G", "E", "F",
    "H", "I", "P", "Q", "R", "S"
]

COG_DESCRIPTIONS = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


def load_cognames(cognames_path: str) -> dict:
    """
    Load cognames2003-2014.tab into a dict: COG_ID -> functional_letter(s).
    Format: COG0001  J  Translation...
    """
    cog_map = {}
    with open(cognames_path, encoding='latin-1') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                cog_id = parts[0]
                func   = parts[1]
                cog_map[cog_id] = func
    return cog_map


def best_hit_from_blast(blast_file: str) -> dict:
    """
    Parse RPS-BLAST output (default BLAST text format, not tabular).
    Returns dict: query_id -> best_COG_id (lowest e-value).

    Handles the text-based BLAST output produced by:
        rpsblast -db Cog -query *.faa -out results.faa
    """
    best_hits = {}        # query -> (evalue_exp, evalue_coef, cog_id)
    current_query = None
    best_evalue  = (0, 1000)    # (exponent, coefficient) — smaller = better
    best_cog     = None

    def _parse_evalue(evalue_str: str):
        """Parse e-value like '1e-12' or '2.5e-34' into (exp, coef) tuple."""
        evalue_str = evalue_str.strip()
        if "e" in evalue_str.lower():
            parts = re.split(r"[eE]", evalue_str)
            coef = float(parts[0]) if parts[0] not in ("-", "") else 1.0
            exp  = int(parts[1])
            return (abs(exp), coef)   # higher abs(exp) = smaller e-value
        try:
            val = float(evalue_str)
            return (0, val)
        except ValueError:
            return (0, 1000)

    def _save_hit(query, cog, evalue_tuple):
        if query is None or cog is None:
            return
        prev = best_hits.get(query)
        if prev is None:
            best_hits[query] = (evalue_tuple, cog)
        else:
            prev_eval, _ = prev
            # Compare: larger exponent wins; on tie, smaller coefficient wins
            if (evalue_tuple[0] > prev_eval[0] or
                    (evalue_tuple[0] == prev_eval[0] and
                     evalue_tuple[1] <= prev_eval[1])):
                best_hits[query] = (evalue_tuple, cog)

    with open(blast_file) as fh:
        for line in fh:
            line = line.rstrip()

            # New query
            if line.startswith("Query="):
                if current_query is not None:
                    _save_hit(current_query, best_cog, best_evalue)
                current_query = line.split("=", 1)[1].strip().split()[0]
                best_evalue = (0, 1000)
                best_cog    = None

            # Hit line with CDD accession and e-value
            elif "CDD" in line and not line.startswith(">"):
                parts = line.split()
                if len(parts) >= 2:
                    cog_id    = parts[1].rstrip(",")
                    evalue_s  = parts[-1]
                    ev_tuple  = _parse_evalue(evalue_s)
                    if best_cog is None or (
                        ev_tuple[0] > best_evalue[0] or
                        (ev_tuple[0] == best_evalue[0] and
                         ev_tuple[1] <= best_evalue[1])
                    ):
                        best_evalue = ev_tuple
                        best_cog    = cog_id

            # Reset on "No hits found"
            elif "found" in line.lower():
                best_evalue = (0, 1000)
                best_cog    = None

    # Save last query
    if current_query is not None:
        _save_hit(current_query, best_cog, best_evalue)

    # Return only COG IDs
    return {q: info[1] for q, info in best_hits.items() if info[1]}


def count_categories(best_hits: dict, cog_map: dict) -> dict:
    """Count COG category occurrences from best-hit assignments."""
    counts = defaultdict(int)
    for cog_id in best_hits.values():
        func_letters = cog_map.get(cog_id, "")
        for letter in func_letters:
            if letter in COG_CATEGORIES:
                counts[letter] += 1
    return counts


def write_csv(counts_per_file: list, output_path: str):
    """
    Write one row per input file with counts for each COG category.
    Compatible with the original Perl output format consumed by Biorem_report.Rmd.
    """
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COG_CATEGORIES)      # header
        for counts in counts_per_file:
            writer.writerow([counts.get(cat, 0) for cat in COG_CATEGORIES])
    print(f"[cog_counter.py] Results written to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Count COG functional categories from RPS-BLAST results."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--blast",  help="Single RPS-BLAST output file")
    group.add_argument("--list",   help="Text file with one BLAST result path per line")

    parser.add_argument("--cognames", required=True,
                        help="cognames2003-2014.tab file")
    parser.add_argument("--output",   required=True,
                        help="Output CSV file (cog_frequencies.csv)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print per-file category totals")
    args = parser.parse_args()

    # Collect input files
    if args.blast:
        blast_files = [args.blast]
    else:
        list_path = Path(args.list)
        if not list_path.exists():
            sys.exit(f"ERROR: List file not found: {args.list}")
        blast_files = [l.strip() for l in list_path.read_text().splitlines()
                       if l.strip()]

    if not Path(args.cognames).exists():
        sys.exit(f"ERROR: COG names file not found: {args.cognames}")

    print(f"[cog_counter.py] Loading COG names from: {args.cognames}")
    cog_map = load_cognames(args.cognames)
    print(f"[cog_counter.py] Loaded {len(cog_map):,} COG entries")

    counts_per_file = []
    for bf in blast_files:
        if not Path(bf).exists():
            print(f"  WARNING: File not found, skipping: {bf}", file=sys.stderr)
            counts_per_file.append({})
            continue

        print(f"[cog_counter.py] Processing: {bf}")
        best_hits = best_hit_from_blast(bf)
        counts    = count_categories(best_hits, cog_map)
        counts_per_file.append(counts)

        if args.verbose:
            total = sum(counts.values())
            print(f"  Queries with COG hit: {len(best_hits):,} | "
                  f"Categorized: {total:,}")
            for cat in COG_CATEGORIES:
                if counts.get(cat, 0) > 0:
                    print(f"    {cat} ({COG_DESCRIPTIONS[cat][:40]}): "
                          f"{counts[cat]}")

    write_csv(counts_per_file, args.output)


if __name__ == "__main__":
    main()
