#!/usr/bin/env python3
"""
BioRemmer — plastic_heatmap_curated.py

Curated version for BioRemmer reports.

Main idea:
  - HMMER still detects all significant hits.
  - The main report/heatmap only shows curated candidates.
  - Promiscuous PFAMs and biologically non-specific annotations are excluded
    from the main report and saved as supplementary tables.

Outputs:
  Results/HMMER/plastic_heatmap_summary_curated.csv
  Results/HMMER/plastic_heatmap_summary_all.csv
  Results/HMMER/plastic_heatmap_summary_excluded.csv
  Results/HMMER/plastic_heatmap_summary.csv   # legacy name, same as curated
  Results/Plots/plastic_degrading_heatmap.png
  Results/Plots/plastic_degrading_heatmap_all_hits.png

Usage:
  python3 scripts/plastic_heatmap_curated.py --results Results --output Results/Plots
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Patch

    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "legend.title_fontsize": 9,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })
except ImportError as e:
    sys.exit(
        f"ERROR: Missing dependency — {e}\n"
        "Install: conda run -n bioremmer_core pip install pandas matplotlib numpy"
    )


PLASTICS = ["PS", "PET", "PE", "PUR", "IP", "PA", "PBAT", "PHB", "PLA"]

PLASTIC_FULL = {
    "PS": "Polystyrene",
    "PET": "Polyethylene terephthalate",
    "PE": "Polyethylene",
    "PUR": "Polyurethane",
    "IP": "Isoprene",
    "PA": "Polyamide",
    "PBAT": "Polybutylene adipate terephthalate",
    "PHB": "Polyhydroxybutyrate",
    "PLA": "Polylactic acid",
}

COG_LEVELS = ["J","A","K","L","B","D","Y","V","T","M","N","Z",
              "W","U","O","X","C","G","E","F","H","I","P","Q","R","S"]

COG_DESCRIPTIONS = {
    "J": "Translation, ribosomal structure",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, division",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction",
    "M": "Cell wall/membrane biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking",
    "O": "Posttranslational modification",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate metabolism",
    "E": "Amino acid metabolism",
    "F": "Nucleotide metabolism",
    "H": "Coenzyme metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport",
    "Q": "Secondary metabolites biosynthesis",
    "R": "General function prediction",
    "S": "Function unknown",
}

PLASTIC_COGS = {"I", "Q", "C", "E", "G", "P"}

# Always excluded from the main report because they are too broad/promiscuous.
# They are still preserved in plastic_heatmap_summary_all.csv and excluded.csv.
PROMISCUOUS_PFAMS = {
    "PF00135",
    "PF00144",
    "PF01425",
    "PF03576",
    "PF13472",
    "PF00561",
    "PF12146",
    "PF12697",
}

# Always excluded from the main report when present in annotation/description.
EXCLUDE_KEYWORDS = [
    "penicillin-binding",
    "beta-lactamase",
    "lactamase",
    "amidotransferase",
    "aminopeptidase",
    "d-alanyl-d-alanine",
    "housekeeping",
    "ribosomal",
    "trna",
    "glutamyl-trna",
    "hypothetical protein",
]

# Positive-support terms. These do not rescue excluded PFAMs for the main table,
# but help assign confidence among retained candidates.
SUPPORTIVE_KEYWORDS = [
    "cutinase",
    "petase",
    "mhetase",
    "polyesterase",
    "depolymerase",
    "carboxylesterase",
    "esterase",
    "lipase",
    "alkane hydroxylase",
    "alkb",
    "monooxygenase",
    "dioxygenase",
    "laccase",
    "multicopper oxidase",
    "peroxidase",
    "haloalkane dehalogenase",
    "epoxide hydrolase",
    "polyhydroxybutyrate depolymerase",
    "phb depolymerase",
]


def pfam_id(hmm_name: str) -> str:
    m = re.search(r"(PF\d+)", str(hmm_name))
    return m.group(1) if m else str(hmm_name)


def contains_any(text: str, terms) -> bool:
    text = str(text).lower()
    return any(term.lower() in text for term in terms)


def parse_hmmer(hmmer_dir: Path) -> "pd.DataFrame":
    rows = []
    for plastic in PLASTICS:
        f = hmmer_dir / f"{plastic}_search.txt"
        if not f.exists():
            continue
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 9:
                    continue

                desc = " ".join(parts[18:]) if len(parts) > 18 else ""
                try:
                    evalue = float(parts[4])
                    score = float(parts[5])
                    dom_evalue = float(parts[7])
                    dom_score = float(parts[8])
                except ValueError:
                    continue

                rows.append({
                    "protein": parts[0],
                    "plastic": plastic,
                    "evalue": evalue,
                    "score": score,
                    "dom_evalue": dom_evalue,
                    "dom_score": dom_score,
                    "hmm": parts[2],
                    "pfam": pfam_id(parts[2]),
                    "desc": desc,
                })
    return pd.DataFrame(rows)


def parse_taxonomy(tax_file: Path, protein_ids: set) -> dict:
    tax_map = {}
    id_pat = re.compile(r"gnl\|X\|(\S+)")
    with open(tax_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            m = id_pat.match(parts[0])
            if not m:
                continue
            pid = m.group(1)
            if pid not in protein_ids:
                continue
            taxonomy = parts[1].strip() if len(parts) > 1 else ""
            if not taxonomy:
                tax_map[pid] = "Unclassified"
                continue
            levels = [t.strip() for t in taxonomy.rstrip(";").split(";") if t.strip()]
            phylum = levels[1] if len(levels) >= 2 else (levels[0] if levels else "Unclassified")
            tax_map[pid] = phylum
    return tax_map


def parse_prokka_tsv(prokka_tsv: Path, protein_ids: set) -> dict:
    annot = {}
    if not prokka_tsv.exists():
        return annot
    with open(prokka_tsv) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            locus, ftype = parts[0], parts[1]
            if ftype == "CDS" and locus in protein_ids:
                annot[locus] = {
                    "product": (parts[6] or "hypothetical protein").strip(),
                    "cog": (parts[5] if len(parts) > 5 else "").strip(),
                }
    return annot


def parse_cog_blast(blast_file: Path, protein_ids: set) -> dict:
    cog_hits = {}
    current_query = None
    cdd_re = re.compile(r'CDD:\d+\s+(COG\d+),\s*([^,]+),\s*(.+?)\s+([\d.e+-]+)\s+([\d.e+-]+)\s*$')

    if not blast_file.exists():
        return cog_hits

    with open(blast_file) as fh:
        for line in fh:
            m_q = re.match(r'^Query=\s+(\S+)', line)
            if m_q:
                current_query = m_q.group(1)
                continue
            if current_query not in protein_ids:
                continue
            m_c = cdd_re.match(line.strip())
            if m_c and current_query not in cog_hits:
                try:
                    evalue = float(m_c.group(5))
                except ValueError:
                    evalue = 1.0
                cog_hits[current_query] = {
                    "cog_id": m_c.group(1),
                    "cog_name": m_c.group(2).strip(),
                    "cog_desc": m_c.group(3).strip(),
                    "evalue": evalue,
                }
    return cog_hits


def parse_cog_frequencies(cog_csv: Path) -> dict:
    if not cog_csv.exists():
        return {cat: 0 for cat in COG_LEVELS}
    with open(cog_csv) as fh:
        reader = csv.DictReader(fh)
        row = next(reader, None)
        if row is None:
            return {cat: 0 for cat in COG_LEVELS}
        return {cat: int(float(row.get(cat, 0) or 0)) for cat in COG_LEVELS}


def classify_hits(df: "pd.DataFrame", annot_map: dict) -> "pd.DataFrame":
    rows = []

    for _, r in df.iterrows():
        pid = r["protein"]
        annotation = annot_map.get(pid, {}).get("product", r.get("desc", "hypothetical protein"))
        text_blob = f"{annotation} {r.get('desc', '')} {r.get('hmm', '')}".lower()

        pfam = r["pfam"]
        is_promiscuous = pfam in PROMISCUOUS_PFAMS
        has_exclusion = contains_any(text_blob, EXCLUDE_KEYWORDS)
        has_supportive = contains_any(text_blob, SUPPORTIVE_KEYWORDS)

        strong_hmmer = (r["evalue"] <= 1e-20 and r["dom_evalue"] <= 1e-12 and r["score"] >= 70)
        very_strong_hmmer = (r["evalue"] <= 1e-50 and r["dom_evalue"] <= 1e-20 and r["score"] >= 100)

        exclusion_reasons = []
        if is_promiscuous:
            exclusion_reasons.append("promiscuous PFAM excluded from main report")
        if has_exclusion:
            exclusion_reasons.append("excluded annotation keyword")

        keep_curated = len(exclusion_reasons) == 0

        if keep_curated:
            if very_strong_hmmer and has_supportive:
                confidence = "High confidence"
                reason = "Very strong HMMER support and compatible functional annotation."
            elif strong_hmmer and has_supportive:
                confidence = "Medium confidence"
                reason = "Strong HMMER support and compatible functional annotation."
            elif strong_hmmer:
                confidence = "Medium confidence"
                reason = "Strong HMMER support; functional annotation requires manual validation."
            else:
                confidence = "Low confidence"
                reason = "Weak HMMER support; retained only if not excluded, but should be interpreted cautiously."
        else:
            confidence = "Excluded from main report"
            reason = "; ".join(exclusion_reasons)

        rr = r.to_dict()
        rr.update({
            "annotation_for_filter": annotation,
            "Promiscuous PFAM": bool(is_promiscuous),
            "Excluded keyword": bool(has_exclusion),
            "Supportive annotation": bool(has_supportive),
            "Keep curated": bool(keep_curated),
            "Confidence": confidence,
            "Classification reason": reason,
        })
        rows.append(rr)

    return pd.DataFrame(rows)


def build_heatmap(df, tax_map, annot_map, cog_blast, output_path, title_suffix=""):
    if df.empty:
        print(f"[plastic_heatmap_curated.py] No rows available for heatmap: {output_path}")
        return False

    df = df.copy()
    df["-log10_evalue"] = df["evalue"].apply(lambda e: min(-math.log10(e), 60) if e > 0 else 60)
    pivot = df.pivot_table(
        index="protein", columns="plastic",
        values="-log10_evalue", aggfunc="max"
    ).reindex(columns=PLASTICS)

    pivot["n_plastics"] = pivot.notna().sum(axis=1)
    pivot = pivot.sort_values(["n_plastics", "protein"], ascending=[False, True])
    pivot = pivot.drop(columns="n_plastics")

    proteins = pivot.index.tolist()
    n_prot = len(proteins)

    phyla = [tax_map.get(p, "No taxonomy data") for p in proteins]
    unique_phyla = sorted(set(phyla))
    cmap_tax = matplotlib.colormaps.get_cmap("Set2").resampled(max(len(unique_phyla), 3))
    phylum_color = {ph: cmap_tax(i) for i, ph in enumerate(unique_phyla)}
    tax_colors = [phylum_color[ph] for ph in phyla]

    ylabels = []
    for p in proteins:
        ann = annot_map.get(p, {})
        prod = ann.get("product", "hypothetical protein")
        cb = cog_blast.get(p, {})
        cog_id = cb.get("cog_id", "")
        cog_name = cb.get("cog_name", "")
        cog_str = f" [{cog_id}: {cog_name[:28]}]" if cog_id and cog_name else (f" [{cog_id}]" if cog_id else "")
        short = prod[:32] + ("…" if len(prod) > 32 else "")
        ylabels.append(f"{p}  {short}{cog_str}" if short else f"{p}{cog_str}")

    fig_h = max(6, n_prot * 0.45 + 4.5)
    fig, ax = plt.subplots(figsize=(14.5, fig_h))

    mat = pivot.values.copy()
    masked = np.ma.masked_invalid(mat)
    cmap_heat = plt.cm.YlOrRd.copy()
    cmap_heat.set_bad(color="#f0f4f8")
    im = ax.imshow(masked, aspect="auto", cmap=cmap_heat, vmin=0, vmax=60, interpolation="nearest")

    for i in range(n_prot):
        for j, pl in enumerate(PLASTICS):
            val = mat[i, j]
            if not np.isnan(val):
                orig_e = df[(df["protein"] == proteins[i]) & (df["plastic"] == pl)]["evalue"].min()
                txt = f"{orig_e:.0e}".replace("e-0", "e-").replace("e+0", "e+")
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=6.5, color="black" if val < 40 else "white",
                        fontweight="bold")

    ax.set_xticks(range(len(PLASTICS)))
    ax.set_xticklabels([PLASTIC_FULL[p] for p in PLASTICS],
                       fontsize=8.5, fontweight="bold", rotation=45, ha="right", va="top")
    ax.tick_params(axis="x", pad=4)
    ax.xaxis.set_ticks_position("bottom")
    ax.xaxis.set_label_position("bottom")

    ax.set_yticks(range(n_prot))
    ax.set_yticklabels(ylabels, fontsize=7.5)
    ax.tick_params(axis="y", pad=22)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("-log10(E-value)\n(higher = more significant)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    strip_x = -0.04
    strip_width = 0.025
    for i, col in enumerate(tax_colors):
        ax.add_patch(plt.Rectangle((strip_x, i - 0.48), strip_width, 0.96,
                                   color=col, transform=ax.get_yaxis_transform(),
                                   clip_on=False, zorder=3))

    legend_handles = [Patch(facecolor=phylum_color[ph], label=ph) for ph in unique_phyla]
    ax.legend(handles=legend_handles, title="Phylum (Metaxa2)",
              loc="upper center", bbox_to_anchor=(0.5, -0.26),
              ncol=min(5, len(unique_phyla)), fontsize=7.5, title_fontsize=8.5,
              framealpha=0.9, edgecolor="#cccccc")

    ax.set_xticks(np.arange(-0.5, len(PLASTICS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_prot, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    plt.title(f"Potential Plastic-Degrading Enzymes — Taxonomic Origin{title_suffix}",
              fontsize=12, fontweight="bold", pad=10)
    fig.subplots_adjust(left=0.38, right=0.96, top=0.94, bottom=0.22)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"[plastic_heatmap_curated.py] Heatmap saved: {output_path}")
    return True


def build_cog_highlighted(cog_freq, hmmer_cogs, output_path):
    cats = COG_LEVELS
    counts = [cog_freq.get(c, 0) for c in cats]
    colors = []
    for cat in cats:
        if cat in hmmer_cogs and cat in PLASTIC_COGS:
            colors.append("#d62728")
        elif cat in PLASTIC_COGS:
            colors.append("#2c7bb6")
        else:
            colors.append("#b0b8c1")

    ymax = max(counts) if counts else 1
    ymax = ymax if ymax > 0 else 1

    fig, ax = plt.subplots(figsize=(13.5, 7.8))
    ax.bar(range(len(cats)), counts, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_ylim(0, ymax * 1.20)

    for i, cat in enumerate(cats):
        if cat in hmmer_cogs:
            ax.text(i, counts[i] + ymax * 0.03, "★",
                    ha="center", va="bottom", fontsize=10, color="#d62728")

    ax.set_xlabel("COG Category", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_ylabel("Frequency", fontsize=11, fontweight="bold")
    ax.set_title("Functional Profiling — Plastic Degradation Context",
                 fontsize=12, fontweight="bold", pad=18)
    ax.spines[["top", "right"]].set_visible(False)

    legend_elements = [
        Patch(facecolor="#d62728", label="Plastic-related COG with enzyme hits (★)"),
        Patch(facecolor="#2c7bb6", label="Plastic-related COG (no hits in this sample)"),
        Patch(facecolor="#b0b8c1", label="Other COG categories"),
    ]
    ax.legend(handles=legend_elements, fontsize=8.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.62), ncol=1, framealpha=0.9,
              edgecolor="#cccccc", title="Legend", title_fontsize=9,
              bbox_transform=ax.transAxes)

    desc_labels = [COG_DESCRIPTIONS.get(c, "") for c in cats]
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(desc_labels, fontsize=7, rotation=45, ha="right", va="top")
    ax.tick_params(axis="x", pad=4)

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(range(len(cats)))
    ax_top.set_xticklabels(cats, fontsize=9)
    ax_top.xaxis.set_ticks_position("top")
    ax_top.tick_params(axis="x", pad=2)
    ax_top.set_xlabel("COG category", fontsize=10, fontweight="bold", labelpad=6)

    fig.subplots_adjust(left=0.07, right=0.97, top=0.88, bottom=0.44)
    fig.savefig(output_path, dpi=400, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"[plastic_heatmap_curated.py] COG figure saved: {output_path}")


def write_summary_csv(df, tax_map, annot_map, cog_blast_map, out_csv):
    rows = []
    for _, row in df.iterrows():
        pid = row["protein"]
        ann = annot_map.get(pid, {})
        cb = cog_blast_map.get(pid, {})
        tax = tax_map.get(pid, "Unclassified")
        rows.append({
            "Protein ID": pid,
            "Plastic type": row["plastic"],
            "HMM profile": row["hmm"],
            "PFAM": row["pfam"],
            "E-value": row["evalue"],
            "Score": row["score"],
            "Domain E-value": row["dom_evalue"],
            "Domain score": row["dom_score"],
            "Annotation": ann.get("product", "hypothetical protein"),
            "COG ID": cb.get("cog_id", "—"),
            "COG function": cb.get("cog_name", "—"),
            "Phylum": tax,
            "Confidence": row.get("Confidence", "Not classified"),
            "Classification reason": row.get("Classification reason", ""),
            "Keep curated": row.get("Keep curated", False),
            "Promiscuous PFAM": row.get("Promiscuous PFAM", False),
            "Excluded keyword": row.get("Excluded keyword", False),
            "Supportive annotation": row.get("Supportive annotation", False),
        })

    with open(out_csv, "w", newline="") as fh:
        if rows:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        else:
            writer = csv.writer(fh)
            writer.writerow([
                "Protein ID", "Plastic type", "HMM profile", "PFAM", "E-value",
                "Annotation", "Confidence", "Classification reason"
            ])
    print(f"[plastic_heatmap_curated.py] Summary CSV: {out_csv}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", default=None,
                        help="Output directory for plots. Default: <results>/Plots")
    parser.add_argument("--use-all-heatmap", action="store_true",
                        help="Use all hits in heatmap. Default uses curated hits only.")
    args = parser.parse_args()

    results = Path(args.results)
    outdir = Path(args.output) if args.output else (results / "Plots")
    outdir.mkdir(parents=True, exist_ok=True)

    hmmer_dir = results / "HMMER"
    tax_file = results / "Metaxa2_results" / "metaxa2.taxonomy.txt"
    prokka_tsv = results / "Prokka_results" / "prokka_results.tsv"
    cog_blast = results / "COG" / "COGs_results.faa"
    cog_csv = results / "COG" / "cog_frequencies.csv"

    print("[plastic_heatmap_curated.py] Parsing HMMER results...")
    df = parse_hmmer(hmmer_dir)
    if df.empty:
        sys.exit("ERROR: No HMMER hits found.")

    protein_ids = set(df["protein"].unique())
    print(f"  Raw HMMER: {len(df)} hits | {len(protein_ids)} proteins | {df['plastic'].nunique()} plastic types")

    print("[plastic_heatmap_curated.py] Parsing Metaxa2 taxonomy...")
    tax_map = parse_taxonomy(tax_file, protein_ids) if tax_file.exists() else {}

    print("[plastic_heatmap_curated.py] Parsing Prokka annotations...")
    annot_map = parse_prokka_tsv(prokka_tsv, protein_ids)

    print("[plastic_heatmap_curated.py] Parsing COG blast results...")
    cog_blast_map = parse_cog_blast(cog_blast, protein_ids) if cog_blast.exists() else {}

    print("[plastic_heatmap_curated.py] Classifying and filtering hits...")
    classified = classify_hits(df, annot_map)
    curated = classified[classified["Keep curated"] == True].copy()
    excluded = classified[classified["Keep curated"] == False].copy()

    print(f"  All classified hits: {len(classified)}")
    print(f"  Curated main-report hits: {len(curated)}")
    print(f"  Excluded supplementary hits: {len(excluded)}")

    all_csv = hmmer_dir / "plastic_heatmap_summary_all.csv"
    curated_csv = hmmer_dir / "plastic_heatmap_summary_curated.csv"
    excluded_csv = hmmer_dir / "plastic_heatmap_summary_excluded.csv"
    legacy_csv = hmmer_dir / "plastic_heatmap_summary.csv"

    write_summary_csv(classified, tax_map, annot_map, cog_blast_map, all_csv)
    write_summary_csv(curated, tax_map, annot_map, cog_blast_map, curated_csv)
    write_summary_csv(excluded, tax_map, annot_map, cog_blast_map, excluded_csv)
    write_summary_csv(curated, tax_map, annot_map, cog_blast_map, legacy_csv)

    heatmap_df = classified if args.use_all_heatmap else curated

    if heatmap_df.empty:
        print("[plastic_heatmap_curated.py] WARNING: no curated hits. Main heatmap will not be generated.")
    else:
        build_heatmap(
            heatmap_df,
            tax_map,
            annot_map,
            cog_blast_map,
            outdir / "plastic_degrading_heatmap.png",
            title_suffix=" (curated candidates)" if not args.use_all_heatmap else " (all hits)",
        )

    build_heatmap(
        classified,
        tax_map,
        annot_map,
        cog_blast_map,
        outdir / "plastic_degrading_heatmap_all_hits.png",
        title_suffix=" (all HMMER hits)",
    )

    cog_freq = parse_cog_frequencies(cog_csv)
    hmmer_cogs = set()
    for cat in PLASTIC_COGS:
        if cog_freq.get(cat, 0) > 0:
            hmmer_cogs.add(cat)

    build_cog_highlighted(cog_freq, hmmer_cogs, outdir / "cog_plastic_highlighted.png")

    print("[plastic_heatmap_curated.py] Done.")
    print("Main-report table: Results/HMMER/plastic_heatmap_summary_curated.csv")
    print("Supplementary all-hits table: Results/HMMER/plastic_heatmap_summary_all.csv")
    print("Supplementary excluded table: Results/HMMER/plastic_heatmap_summary_excluded.csv")


if __name__ == "__main__":
    main()
