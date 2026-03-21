#!/usr/bin/env python3
"""
BioRemmer v3.0 — plastic_heatmap.py
Generates:
  1. Heatmap: plastic-degrading enzymes × plastic types (annotated with taxonomy)
  2. COG figure: all COG categories with plastic-related ones highlighted
  3. Summary CSVs for the R report

Usage:
    python3 plastic_heatmap.py --results <Results/> --output <Results/Plots/>
"""
import argparse, os, re, math, sys, csv, textwrap
from pathlib import Path

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import numpy as np
    from matplotlib.patches import Patch
    # Global font settings — uniform across all figures
    matplotlib.rcParams.update({
        "font.family":          "sans-serif",
        "font.sans-serif":      ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":            10,
        "axes.titlesize":       12,
        "axes.titleweight":     "bold",
        "axes.labelsize":       11,
        "axes.labelweight":     "normal",
        "xtick.labelsize":      9,
        "ytick.labelsize":      9,
        "legend.fontsize":      8.5,
        "legend.title_fontsize":9,
        "figure.dpi":           150,
        "savefig.dpi":          300,
        "savefig.bbox":         "tight",
        "savefig.facecolor":    "white",
    })
except ImportError as e:
    sys.exit(f"ERROR: Missing dependency — {e}\n"
             f"Install: conda run -n bioremmer_core pip install pandas matplotlib numpy")

# ── Constants ─────────────────────────────────────────────────────────────────
PLASTICS = ["PS", "PET", "PE", "PUR", "IP", "PA", "PBAT", "PHB", "PLA"]

PLASTIC_FULL = {
    "PS":   "Polystyrene",
    "PET":  "Polyethylene\nterephthalate",
    "PE":   "Polyethylene",
    "PUR":  "Polyurethane",
    "IP":   "Isoprene",
    "PA":   "Polyamide",
    "PBAT": "Polybutylene adipate\nterephthalate",
    "PHB":  "Polyhydroxybutyrate",
    "PLA":  "Polylactic acid",
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

# COG categories functionally linked to plastic degradation
PLASTIC_COGS = {"I", "Q", "C", "E", "G", "P"}

# ── Parsers ───────────────────────────────────────────────────────────────────
def parse_hmmer(hmmer_dir: Path) -> pd.DataFrame:
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
                if len(parts) < 6:
                    continue
                rows.append({
                    "protein": parts[0],
                    "plastic": plastic,
                    "evalue":  float(parts[4]),
                    "hmm":     parts[2],
                    "desc":    " ".join(parts[18:]) if len(parts) > 18 else "",
                })
    return pd.DataFrame(rows)


def parse_taxonomy(tax_file: Path, protein_ids: set) -> dict:
    tax_map = {}
    id_pat = re.compile(r"gnl\|X\|(\S+)")
    with open(tax_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            m = id_pat.match(parts[0])
            if not m:
                continue
            pid = m.group(1)
            if pid not in protein_ids:
                continue
            taxonomy = parts[1].strip() if len(parts) > 1 else ""
            if not taxonomy:
                tax_map[pid] = "Unclassified (no assignment)"
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
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            locus, ftype = parts[0], parts[1]
            product = parts[6] if len(parts) > 6 else ""
            cog_col = parts[5] if len(parts) > 5 else ""
            if ftype == "CDS" and locus in protein_ids:
                annot[locus] = {
                    "product": product.strip() or "hypothetical protein",
                    "cog":     cog_col.strip(),
                }
    return annot


def parse_cog_blast(blast_file: Path, protein_ids: set) -> dict:
    """Extract best COG hit per protein from rpsblast output."""
    cog_hits = {}
    current_query = None
    cdd_re = re.compile(r'CDD:\d+\s+(COG\d+),\s*([^,]+),\s*(.+?)\s+([\d.e+-]+)\s+([\d.e+-]+)\s*$')

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
                    "cog_id":   m_c.group(1),
                    "cog_name": m_c.group(2).strip(),
                    "cog_desc": m_c.group(3).strip(),
                    "evalue":   evalue,
                }
    return cog_hits


def parse_cog_frequencies(cog_csv: Path) -> dict:
    """Load cog_frequencies.csv → dict category_letter -> count."""
    if not cog_csv.exists():
        return {cat: 0 for cat in COG_LEVELS}
    import csv as csvmod
    with open(cog_csv) as fh:
        reader = csvmod.DictReader(fh)
        row = next(reader, None)
        if row is None:
            return {cat: 0 for cat in COG_LEVELS}
        return {cat: int(float(row.get(cat, 0) or 0)) for cat in COG_LEVELS}


# ── Figure 1: Heatmap ─────────────────────────────────────────────────────────
def build_heatmap(df, tax_map, annot_map, cog_blast, output_path):
    df["-log10_evalue"] = df["evalue"].apply(
        lambda e: min(-math.log10(e), 60) if e > 0 else 60)
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
        ann  = annot_map.get(p, {})
        prod = ann.get("product", "hypothetical protein")
        cb   = cog_blast.get(p, {})
        cog_id   = cb.get("cog_id", "")
        cog_name = cb.get("cog_name", "")
        if cog_id and cog_name:
            cog_str = f" [{cog_id}: {cog_name[:28]}]"
        elif cog_id:
            cog_str = f" [{cog_id}]"
        else:
            cog_str = ""
        short = prod[:28] + ("…" if len(prod) > 28 else "")
        if short == "hypothetical protein":
            ylabels.append(f"{p}{cog_str}")
        else:
            ylabels.append(f"{p}  {short}{cog_str}")

    fig_h = max(6, n_prot * 0.45 + 4.5)
    fig, ax = plt.subplots(figsize=(14.5, fig_h))

    mat = pivot.values.copy()
    masked = np.ma.masked_invalid(mat)
    cmap_heat = plt.cm.YlOrRd.copy()
    cmap_heat.set_bad(color="#f0f4f8")
    im = ax.imshow(masked, aspect="auto", cmap=cmap_heat,
                   vmin=0, vmax=60, interpolation="nearest")

    for i in range(n_prot):
        for j, pl in enumerate(PLASTICS):
            val = mat[i, j]
            if not np.isnan(val):
                orig_e = df[(df["protein"] == proteins[i]) &
                            (df["plastic"] == pl)]["evalue"].min()
                txt = f"{orig_e:.0e}".replace("e-0", "e-").replace("e+0", "e+")
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=6.5, color="black" if val < 40 else "white",
                        fontweight="bold")

    ax.set_xticks(range(len(PLASTICS)))
    ax.set_xticklabels(
        [PLASTIC_FULL[p].replace("\n", " ") for p in PLASTICS],
        fontsize=8.5, fontweight="bold",
        rotation=45, ha="right", va="top"
    )
    ax.tick_params(axis="x", pad=4)
    ax.xaxis.set_ticks_position("bottom")
    ax.xaxis.set_label_position("bottom")

    ax.set_yticks(range(n_prot))
    ax.set_yticklabels(ylabels, fontsize=7.5)
    ax.tick_params(axis="y", pad=22)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("-log10(E-value)\n(higher = more significant)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # Taxonomy strip — narrow band between y-tick labels and the axis line
    # get_yaxis_transform: x in axes coords (0=axis line), y in data coords
    strip_x     = -0.04   # just left of the axis line
    strip_width =  0.025  # narrow strip
    for i, col in enumerate(tax_colors):
        ax.add_patch(plt.Rectangle(
            (strip_x, i - 0.48), strip_width, 0.96,
            color=col,
            transform=ax.get_yaxis_transform(),
            clip_on=False,
            zorder=3
        ))

    legend_handles = [Patch(facecolor=phylum_color[ph], label=ph) for ph in unique_phyla]
    ax.legend(handles=legend_handles, title="Phylum (Metaxa2)",
              loc="upper center", bbox_to_anchor=(0.5, -0.26),
              ncol=min(5, len(unique_phyla)), fontsize=7.5, title_fontsize=8.5,
              framealpha=0.9, edgecolor="#cccccc")

    ax.set_xticks(np.arange(-0.5, len(PLASTICS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_prot, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    plt.title("Potential Plastic-Degrading Enzymes — Taxonomic Origin",
              fontsize=12, fontweight="bold", pad=10)
    fig.subplots_adjust(left=0.38, right=0.96, top=0.94, bottom=0.22)
    fig.savefig(output_path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"[plastic_heatmap.py] Heatmap saved: {output_path}")


# ── Figure 2: COG highlighted ─────────────────────────────────────────────────
def build_cog_highlighted(cog_freq, hmmer_cogs, output_path):
    """
    Bar chart of all COG categories.
    Bottom axis: wrapped COG descriptions
    Top axis: single-letter COG categories
    Bars with HMMER hits get a star annotation.
    """
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

    ymax = max(counts) if counts else 0
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
    ax.legend(handles=legend_elements, fontsize=8.5,
              loc="upper center",
              bbox_to_anchor=(0.5, -0.62),
              ncol=1, framealpha=0.9, edgecolor="#cccccc",
              title="Legend", title_fontsize=9,
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
    fig.savefig(output_path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"[plastic_heatmap.py] COG figure saved: {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", default=None,
                        help="Output directory for plots (default: <results>/Plots)")
    args = parser.parse_args()

    results = Path(args.results)
    outdir  = Path(args.output) if args.output else (results / "Plots")
    outdir.mkdir(parents=True, exist_ok=True)

    hmmer_dir  = results / "HMMER"
    tax_file   = results / "Metaxa2_results" / "metaxa2.taxonomy.txt"
    prokka_tsv = results / "Prokka_results"  / "prokka_results.tsv"
    cog_blast  = results / "COG"             / "COGs_results.faa"
    cog_csv    = results / "COG"             / "cog_frequencies.csv"

    print("[plastic_heatmap.py] Parsing HMMER results...")
    df = parse_hmmer(hmmer_dir)
    if df.empty:
        sys.exit("ERROR: No HMMER hits found.")
    protein_ids = set(df["protein"].unique())
    print(f"  {len(df)} hits | {len(protein_ids)} proteins | "
          f"{df['plastic'].nunique()} plastic types")

    print("[plastic_heatmap.py] Parsing Metaxa2 taxonomy...")
    tax_map = parse_taxonomy(tax_file, protein_ids) if tax_file.exists() else {}

    print("[plastic_heatmap.py] Parsing Prokka annotations...")
    annot_map = parse_prokka_tsv(prokka_tsv, protein_ids)

    print("[plastic_heatmap.py] Parsing COG blast results...")
    cog_blast_map = parse_cog_blast(cog_blast, protein_ids) if cog_blast.exists() else {}

    # ── Figure 1: Heatmap ─────────────────────────────────────────────────────
    build_heatmap(df, tax_map, annot_map, cog_blast_map,
                  outdir / "plastic_degrading_heatmap.png")

    # ── Figure 2: COG highlighted ─────────────────────────────────────────────
    cog_freq = parse_cog_frequencies(cog_csv)
    # Which plastic-related COG categories have HMMER hits?
    hmmer_cogs = set()
    for pid in protein_ids:
        ann = annot_map.get(pid, {})
        cog_tsv = ann.get("cog", "")
        if cog_tsv:
            # Extract single-letter category from COG description
            # Prokka TSV col 6 has e.g. "COG0749"
            pass
        cb = cog_blast_map.get(pid, {})
        cat = cb.get("cog_cat", "")
        if cat:
            hmmer_cogs.add(cat)
    # Also flag I, Q, C if they have any count (known plastic-related)
    for cat in PLASTIC_COGS:
        if cog_freq.get(cat, 0) > 0:
            hmmer_cogs.add(cat)

    build_cog_highlighted(cog_freq, hmmer_cogs,
                          outdir / "cog_plastic_highlighted.png")

    # ── CSV summary for R report ───────────────────────────────────────────────
    rows = []
    for _, row in df.iterrows():
        pid  = row["protein"]
        ann  = annot_map.get(pid, {})
        cb   = cog_blast_map.get(pid, {})
        tax  = tax_map.get(pid, "Unclassified")
        rows.append({
            "Protein ID":   pid,
            "Plastic type": row["plastic"],
            "HMM profile":  row["hmm"],
            "E-value":      row["evalue"],
            "Annotation":   ann.get("product", "hypothetical protein"),
            "COG ID":       cb.get("cog_id", "—"),
            "COG function": cb.get("cog_name", "—"),
            "Phylum":       tax,
        })

    summary_csv = results / "HMMER" / "plastic_heatmap_summary.csv"
    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[plastic_heatmap.py] Summary CSV: {summary_csv}")
    print(f"[plastic_heatmap.py] Done — {len(rows)} records")


if __name__ == "__main__":
    main()
