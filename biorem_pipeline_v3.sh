#!/bin/bash
# =============================================================================
#  BioRemmer v3.0 — Main pipeline (with resume support)
#  Usage: ./biorem_pipeline_v2.sh <R1.fastq.gz> <R2.fastq.gz> <sample> <threads>
#
#  Resume: if a step's output already exists, it is skipped automatically.
#  To rerun a specific step, delete its output folder:
#    rm -rf Results/Trimmomatic      → reruns step 2 (fastp)
#    rm -rf Results/SPAdes_results   → reruns step 3
#    rm -rf Results/Prokka_results   → reruns step 4
#    rm -rf Results/MetaBat2         → reruns step 5
#    rm -rf Results/VAMPhyRE         → reruns step 6
#    rm -rf Results/HMMER            → reruns step 7
#    rm -rf Results/COG              → reruns step 8
#    rm -rf Results/Metaxa2_results  → reruns step 9
#    rm -f  Results/Biorem_report.html → reruns step 10
#    rm -f  Results/Plots/plastic_degrading_heatmap.png → reruns step 11
# =============================================================================

if [ $# -ne 4 ]; then
    echo "
  BioRemmer v3.0
  Usage:
    ./biorem_pipeline_v2.sh <R1> <R2> <sample_name> <threads>

  Example:
    ./biorem_pipeline_v2.sh test/Test_short_1.fastq.gz \\
                            test/Test_short_2.fastq.gz \\
                            test_short 4

  Resume: the pipeline skips steps whose output already exists.
  Delete a step's Results/ subfolder to force it to rerun.
    "
    exit 1
fi

# ── Load config ───────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ── Helpers ───────────────────────────────────────────────────────────────────
run_in() {
    local env="$1"; shift
    conda run --no-capture-output -n "$env" "$@"
}

# skip_if_done <description> <file_or_dir_to_check>
# Returns 0 (skip) if the output exists, 1 (run) if it doesn't
skip_if_done() {
    local description="$1"
    local output="$2"
    if [ -e "$output" ]; then
        printf "  [SKIP] %s — output already exists: %s\n" "$description" "$output"
        return 0
    fi
    return 1
}

# ── Pipeline ──────────────────────────────────────────────────────────────────
(
set -euo pipefail

start=$(date +%s)
R1="$1"
R2="$2"
samplename="$3"
threads="$4"
EXT=".fastq.gz"
RESULTS="${DIR}/Results"
    PLOTS="${RESULTS}/Plots"
    LOGS="${RESULTS}/logs"

printf "\n============================================\n"
printf "  BioRemmer v1.0\n"
printf "  Sample  : %s\n" "$samplename"
printf "  Threads : %s\n" "$threads"
printf "  Core env: %s\n" "$CORE_ENV"
printf "  R env   : %s\n" "$R_ENV"
printf "  Output  : %s\n" "$RESULTS"
printf "  Resume  : steps with existing output are skipped\n"
printf "============================================\n\n"

# ── Step 1: Validate inputs ───────────────────────────────────────────────────
printf "\nStep 1: Validating input files\n"
[ -f "$R1" ] || { echo "ERROR: R1 not found: $R1"; exit 1; }
[ -f "$R2" ] || { echo "ERROR: R2 not found: $R2"; exit 1; }
printf "  R1: %s\n  R2: %s\n" "$R1" "$R2"

# ── Step 2: Quality control — fastp ──────────────────────────────────────────
printf "\nStep 2: Quality control — fastp\n"
mkdir -p "${RESULTS}/Trimmomatic" "${RESULTS}/QC"

if ! skip_if_done "fastp" "${RESULTS}/Trimmomatic/fastp_summary.json"; then
    run_in "$CORE_ENV" fastp \
        -i "$R1" -I "$R2" \
        -o "${RESULTS}/Trimmomatic/${samplename}.t_R1${EXT}" \
        -O "${RESULTS}/Trimmomatic/${samplename}.t_R2${EXT}" \
        --unpaired1 "${RESULTS}/Trimmomatic/${samplename}.u_R1${EXT}" \
        --unpaired2 "${RESULTS}/Trimmomatic/${samplename}.u_R2${EXT}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 15 \
        --length_required 36 \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 15 \
        --thread "$threads" \
        --json "${RESULTS}/Trimmomatic/fastp_summary.json" \
        --html "${RESULTS}/QC/fastp_report.html" \
        --report_title "${samplename} — fastp QC report"
fi

# ── Step 3: Metagenome assembly — metaSPAdes ──────────────────────────────────
printf "\nStep 3: Metagenome assembly — metaSPAdes\n"
if ! skip_if_done "metaSPAdes" "${RESULTS}/SPAdes_results/contigs.fasta"; then
    run_in "$CORE_ENV" metaspades.py \
        -t "$threads" \
        -1 "${RESULTS}/Trimmomatic/${samplename}.t_R1${EXT}" \
        -2 "${RESULTS}/Trimmomatic/${samplename}.t_R2${EXT}" \
        -o "${RESULTS}/SPAdes_results"
fi

# ── Step 4: Functional annotation — Prokka ────────────────────────────────────
printf "\nStep 4: Functional annotation — Prokka\n"
if ! skip_if_done "Prokka" "${RESULTS}/Prokka_results/prokka_results.faa"; then
    run_in "$CORE_ENV" prokka \
        --outdir "${RESULTS}/Prokka_results" \
        "${RESULTS}/SPAdes_results/contigs.fasta" \
        --kingdom Bacteria \
        --cpus "$threads" \
        --metagenome \
        --prefix "prokka_results" \
        --fast
fi

# ── Step 5: MAG binning — MetaBat2 ───────────────────────────────────────────
printf "\nStep 5: MAG binning — MetaBat2\n"
mkdir -p "${RESULTS}/MetaBat2"
if ! skip_if_done "MetaBat2" "${RESULTS}/MetaBat2/genome_list.txt"; then
    run_in "$CORE_ENV" metabat2 \
        -i "${RESULTS}/SPAdes_results/contigs.fasta" \
        -o "${RESULTS}/MetaBat2/Metabat2_results" \
        -m 1500 -t "$threads"

    ls -1 "${RESULTS}/MetaBat2/Metabat2_results".*.fa \
        > "${RESULTS}/MetaBat2/genome_list.txt" 2>/dev/null || \
        touch "${RESULTS}/MetaBat2/genome_list.txt"
fi

# ── Step 6: Phylogenomics — VAMPhyRE + build_tree.R ──────────────────────────
printf "\nStep 6: Phylogenomics — VAMPhyRE + NJ tree\n"
mkdir -p "${RESULTS}/VAMPhyRE"
chmod +x "${vamphyre}/VH5cmdl" \
         "${vamphyre}/MergeVGF" \
         "${vamphyre}/VFAT" 2>/dev/null || true

if ! skip_if_done "VAMPhyRE + NJ tree" "${RESULTS}/VAMPhyRE/VGF_tree.nwk"; then

    # Generar lista de genomas de referencia
    ls -1 "${DIR}/bin/vamphyre/genomes/genome_DB/"*.fna \
        > "${RESULTS}/VAMPhyRE/genome_list.txt"
    n_genomes=$(wc -l < "${RESULTS}/VAMPhyRE/genome_list.txt")
    printf "  Reference genomes found: %s\n" "$n_genomes"
    if [ "$n_genomes" -lt 3 ]; then
        printf "  ERROR: Less than 3 reference genomes in genome_DB\n"
        exit 1
    fi

    if [ -s "${RESULTS}/MetaBat2/genome_list.txt" ]; then
        printf "  MAG bins found — running full VAMPhyRE merge with reference DB\n"

        # 6a: VH5cmdl sobre MAGs
        if ! skip_if_done "VH5cmdl (MAGs)" "${RESULTS}/VAMPhyRE/VGF_13mer.txt"; then
            "${vamphyre}/VH5cmdl" \
                -PROBEFILE  "$probe" \
                -TARGETLIST "${RESULTS}/MetaBat2/genome_list.txt" \
                -OUTFILE    "${RESULTS}/VAMPhyRE/VGF_13mer.txt" \
                -MISMATCHES 1 -STRAND both
        fi

        # 6b: MergeVGF — combina MAGs + DB de referencia
        if ! skip_if_done "MergeVGF" "${RESULTS}/VAMPhyRE/VGF_Merged.txt"; then
            "${vamphyre}/MergeVGF" \
                "${RESULTS}/VAMPhyRE/VGF_13mer.txt" \
                "$genomeDB" \
                "${RESULTS}/VAMPhyRE/VGF_Merged.txt"
        fi

        # 6c: Concatenar listas de genomas
        cat "${RESULTS}/MetaBat2/genome_list.txt" "${RESULTS}/VAMPhyRE/genome_list.txt" \
            > "${RESULTS}/VAMPhyRE/genome_list_conc.txt"

        # 6d: VFAT — matriz de distancias con MAGs + referencia
        if ! skip_if_done "VFAT (merged)" "${RESULTS}/VAMPhyRE/VGFdistances.meg"; then
            "${vamphyre}/VFAT" \
                -VHFILE     "${RESULTS}/VAMPhyRE/VGF_Merged.txt" \
                -TARGETLIST "${RESULTS}/VAMPhyRE/genome_list_conc.txt" \
                -OUTFILE    "${RESULTS}/VAMPhyRE/VGFdistances" \
                -LEFTEXT 6 -RIGHTEXT 6 -THRESHOLD 23 \
                -MODE DISTANCE -TRACKING YES -FORMAT MEGA
            [ -f track.txt ] && mv track.txt "${RESULTS}/VAMPhyRE/"
            run_in "$CORE_ENV" python3 "${cog_counter}/fix_meg.py" \
                "${RESULTS}/VAMPhyRE/VGFdistances.meg"
        fi

        # 6e: Árbol NJ desde distancias fusionadas
        run_in "$R_ENV" Rscript "${cog_counter}/build_tree.R" \
            "${RESULTS}/VAMPhyRE/VGFdistances.meg" \
            "${RESULTS}/VAMPhyRE/VGF_tree.nwk"

    else
        printf "  WARNING: No MAG bins — building reference-only tree\n"

        # 6a: VH5cmdl sobre genomas de referencia
        if ! skip_if_done "VH5cmdl (reference DB)" "${RESULTS}/VAMPhyRE/VGF_13mer_genomeDB.txt"; then
            "${vamphyre}/VH5cmdl" \
                -PROBEFILE  "$probe" \
                -TARGETLIST "${RESULTS}/VAMPhyRE/genome_list.txt" \
                -OUTFILE    "${RESULTS}/VAMPhyRE/VGF_13mer_genomeDB.txt" \
                -MISMATCHES 1 -STRAND both
        fi

        # 6b: VFAT — matriz de distancias solo referencia
        if ! skip_if_done "VFAT (reference DB)" "${RESULTS}/VAMPhyRE/VGFdistances_BioRemDB.meg"; then
            "${vamphyre}/VFAT" \
                -VHFILE     "${RESULTS}/VAMPhyRE/VGF_13mer_genomeDB.txt" \
                -TARGETLIST "${RESULTS}/VAMPhyRE/genome_list.txt" \
                -OUTFILE    "${RESULTS}/VAMPhyRE/VGFdistances_BioRemDB" \
                -LEFTEXT 6 -RIGHTEXT 6 -THRESHOLD 23 \
                -MODE DISTANCE -FORMAT MEGA
            run_in "$CORE_ENV" python3 "${cog_counter}/fix_meg.py" \
                "${RESULTS}/VAMPhyRE/VGFdistances_BioRemDB.meg"
        fi

        # 6c: Árbol NJ desde distancias de referencia
        run_in "$R_ENV" Rscript "${cog_counter}/build_tree.R" \
            "${RESULTS}/VAMPhyRE/VGFdistances_BioRemDB.meg" \
            "${RESULTS}/VAMPhyRE/VGF_tree.nwk"
    fi
fi

# ── Step 7: HMM enzyme search — HMMER ────────────────────────────────────────
printf "\nStep 7: HMM-based enzyme search — HMMER\n"
mkdir -p "${RESULTS}/HMMER"

for plastic in PET PUR PE PS IP PA PBAT PHB PLA; do
    out="${RESULTS}/HMMER/${plastic}_search.txt"
    if skip_if_done "$plastic HMM" "$out"; then continue; fi

    hmm_file=$(ls "${pfam}/${plastic}/conc_${plastic}.hmm" \
                  "${pfam}/${plastic}/${plastic}.hmm" 2>/dev/null | head -1 || true)
    if [ -z "$hmm_file" ]; then
        printf "  WARNING: No HMM profile for %s — skipping\n" "$plastic"
        continue
    fi
    printf "  Searching %s...\n" "$plastic"
    run_in "$CORE_ENV" hmmsearch \
        --cpu "$threads" -E 1e-12 \
        --tblout "$out" \
        "$hmm_file" \
        "${RESULTS}/Prokka_results/"*.faa
done

# ── Step 8: COG classification — RPS-BLAST + cog_counter.py ──────────────────
printf "\nStep 8: COG classification — RPS-BLAST\n"
mkdir -p "${RESULTS}/COG"

if ! skip_if_done "RPS-BLAST" "${RESULTS}/COG/COGs_results.faa"; then
    run_in "$CORE_ENV" rpsblast \
        -db "$COG_database" \
        -query "${RESULTS}/Prokka_results/"*.faa \
        -evalue 0.0001 \
        -out "${RESULTS}/COG/COGs_results.faa" \
        -num_threads "$threads"
fi

if ! skip_if_done "COG counter" "${RESULTS}/COG/cog_frequencies.csv"; then
    ls -1 "${RESULTS}/COG/COGs_results.faa" \
        > "${RESULTS}/COG/cogs_list.txt"

    run_in "$CORE_ENV" python3 "${cog_counter}/cog_counter.py" \
        --list     "${RESULTS}/COG/cogs_list.txt" \
        --cognames "${cog_counter}/cognames2003-2014.tab" \
        --output   "${RESULTS}/COG/cog_frequencies.csv"
fi

# ── Step 9: 16S rRNA taxonomy — Metaxa2 ──────────────────────────────────────
printf "\nStep 9: Taxonomic assignment — Metaxa2\n"
mkdir -p "${RESULTS}/Metaxa2_results"

if ! skip_if_done "Metaxa2" "${RESULTS}/Metaxa2_results/rarefaction_out.level_6.txt"; then
    run_in "$CORE_ENV" metaxa2 \
        -i "${RESULTS}/Prokka_results/"*.fna \
        --mode m --reltax T --table T -f f -c T --plus T -t b \
        -o "${RESULTS}/Metaxa2_results/metaxa2" --cpu "$threads"

    run_in "$CORE_ENV" metaxa2_ttt \
        -i "${RESULTS}/Metaxa2_results/"*.taxonomy.txt \
        -o "${RESULTS}/Metaxa2_results/taxonomic" -m 7 -t b -n 5

    run_in "$CORE_ENV" metaxa2_dc \
        -i "${RESULTS}/Metaxa2_results/"*level_6.txt \
        -o "${RESULTS}/Metaxa2_results/taxa_combined.txt" -p ".*"

    run_in "$CORE_ENV" metaxa2_si \
        -i "${RESULTS}/Metaxa2_results/"*.taxonomy.txt \
        -o "${RESULTS}/Metaxa2_results/species_inferred.txt" \
        -l 6 --multiple keep --low_identity keep

    run_in "$CORE_ENV" metaxa2_rf \
        -i "${RESULTS}/Metaxa2_results/"*species_inferred.txt \
        -o "${RESULTS}/Metaxa2_results/rarefaction_out" -t b -n 5
fi

# ── Step 10: Report — R Markdown ─────────────────────────────────────────────
printf "\nStep 10: Generating HTML report — R Markdown\n"
mkdir -p "${PLOTS}" "${LOGS}"

if ! skip_if_done "R Markdown report" "${RESULTS}/Biorem_report.html"; then
    run_in "$R_ENV" Rscript "${cog_counter}/report.R" "${DIR}"
fi

# ── Step 11: Plastic degradation heatmap ─────────────────────────────────────
printf "\nStep 11: Plastic degradation heatmap + COG figure\n"

if ! skip_if_done "Heatmap" "${PLOTS}/plastic_degrading_heatmap.png"; then
    run_in "$CORE_ENV" python3 "${cog_counter}/plastic_heatmap.py" \
        --results "${RESULTS}" \
        --output  "${PLOTS}"
fi

# ── Done ──────────────────────────────────────────────────────────────────────
end=$(date +%s)
elapsed=$(( (end - start) / 60 ))
printf "\n============================================\n"
printf "  BioRemmer — Finished!\n"
printf "  Time   : %d minutes\n" "$elapsed"
printf "  Report : %s/Results/Biorem_report.html\n" "$DIR"
printf "============================================\n"

) 2>&1 | tee -a "${LOGS}/bioremmer_${3}.log"
