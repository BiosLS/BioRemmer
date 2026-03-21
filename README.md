# BioRemmer v1.0

**Pipeline for the identification of a functional profile for plastic microbial biodegradation**

BioRemmer is a command-line bioinformatics pipeline that takes raw paired-end metagenomic sequencing data and produces a complete functional profile of plastic-degrading potential, including taxonomic classification, functional annotation, phylogenomics, and identification of putative plastic-degrading enzyme homologs.

---

## Repository structure

```
BioRemmer/
├── biorem_pipeline_v2.sh   # Main pipeline script
├── config.sh               # Tool and database paths (auto-configured)
├── environment.yml         # Conda environment specification
├── install_bioremmer.sh    # Full installation script
├── scripts/
│   ├── report.R            # Report generator (called by pipeline)
│   ├── Biorem_report.Rmd   # R Markdown report template
│   └── COGcounter_v2.pl    # COG frequency counter
├── databases/              # Created by install_bioremmer.sh
│   ├── COG_LE/             # NCBI COG database (little endian)
│   ├── pfam/               # Pfam-A + per-plastic HMM profiles
│   └── vamphyre/           # VAMPhyRE reference genomes and VGF database
├── bin/
│   └── vamphyre/           # VAMPhyRE binaries (manual install)
└── assets/
    └── logo_carta.png      # Logo for HTML report
```

---

## Requirements

- **OS**: Ubuntu 22.04 LTS (recommended) or Ubuntu 24.04
- **RAM**: ≥ 16 GB (32 GB recommended for large metagenomes)
- **Storage**: ≥ 50 GB free (databases + results)
- **CPU**: ≥ 4 threads
- **Conda/Mamba**: installed by `install_bioremmer.sh` if not present

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/BioRemmer/BioRemmer.git
cd BioRemmer
```

### 2. Run the installer

```bash
chmod +x install_bioremmer.sh
./install_bioremmer.sh
```

The installer will:
- Install Miniconda (if not already present)
- Install mamba (fast conda solver)
- Create the `bioremmer` conda environment with all dependencies
- Install MEGA 11 (`.deb`, system-wide)
- Install VAMPhyRE binaries
- Download the NCBI COG database
- Download and build Pfam HMM profiles for each plastic type

> **MEGA 11 and VAMPhyRE** cannot be distributed through conda.  
> If automatic download fails, see the [manual installation](#manual-installation) section below.

### 3. Activate the environment

```bash
conda activate bioremmer
```

---

## Usage

```bash
conda activate bioremmer

./biorem_pipeline_v2.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name> <threads>
```

### Example with demo data

```bash
./biorem_pipeline_v2.sh \
    demo/demo_R1.fastq.gz \
    demo/demo_R2.fastq.gz \
    demo_sample \
    4
```

### Output

All results are written to `Results/` in the working directory:

```
Results/
├── Trimmomatic/
│   ├── trimmomatic_summary.txt
│   └── *.fastq.gz               # Trimmed reads
├── QC/
│   └── *_fastqc.html            # FastQC reports
├── SPAdes_results/
│   └── contigs.fasta
├── Prokka_results/
│   ├── prokka_results.faa       # Predicted proteins
│   └── prokka_results.fna       # Genomic sequences
├── MetaBat2/
│   └── Metabat2_results.*.fa    # MAG bins
├── VAMPhyRE/
│   └── VGF_tree.nwk             # Phylogenomic tree
├── Metaxa2_results/
│   └── rarefaction_out.*        # Taxonomic profiles
├── COG/
│   └── cog_frequencies.csv      # COG category frequencies
├── HMMER/
│   └── *_search.txt             # Enzyme homolog hits per plastic
└── Biorem_report.html           # Complete HTML report
```

---

## Pipeline stages

| Step | Tool | Version | Function |
|------|------|---------|----------|
| 1 | Input validation | — | Checks that R1/R2 files exist |
| 2 | Trimmomatic + FastQC | 0.39 / 0.11.9 | Quality filtering and QC |
| 3 | metaSPAdes | 3.15.4 | Metagenome assembly |
| 4 | Prokka | 1.14.6 | Functional annotation |
| 5 | MetaBat2 | 2.15 | MAG binning |
| 6 | VAMPhyRE + MEGA 11 | 1.0 / 11.0.11 | Phylogenomic tree |
| 7 | HMMER3 | 3.3.2 | HMM-based enzyme search |
| 8 | RPS-BLAST + COGcounter | 2.12 | COG functional classification |
| 9 | Metaxa2 | 2.2.3 | 16S rRNA taxonomic assignment |
| 10 | R Markdown | — | HTML report generation |

---

## BRMR Database

BioRemmer includes the **BRMR (BioRemmer Reference) database**, a curated collection of plastic-degrading enzymes compiled from:

- [PMBD](https://doi.org/10.1093/database/baz119) — Plastics Microbial Biodegradation Database
- [PlasticDB](https://doi.org/10.1093/database/baac008) — Plastic-degrading organisms and proteins
- [PAZy](https://doi.org/10.1002/prot.26325) — Plastics-Active Enzymes Database
- Recent primary literature

Each enzyme record includes: plastic type, synonyms, molecular formula, SMILES, reaction, enzyme name, source organism, gene, UniProt accession, EC number, Pfam family, KEGG/GO annotations, and amino acid sequence.

**Plastic types covered:** PET, PLA, PHB, PUR, PE, PA, PBAT, PS, IP

> **Important note on PS and PE:** The enzyme homologs identified for polystyrene (PS)
> and polyethylene (PE) are based on proteins with reported activity on monomers or
> structurally related substrates *in vitro*. No functionally verified enzymes capable
> of depolymerizing PS or PE polymers under environmental conditions are currently known.
> Results for these plastic types should be interpreted as putative homologs only.

---

## Manual installation

### MEGA 11

```bash
# Download from https://www.megasoftware.net/
wget https://www.megasoftware.net/do_force_download/megacc_11.0.13-1_amd64.deb
sudo dpkg -i megacc_11.0.13-1_amd64.deb
```

### VAMPhyRE

```bash
# Download from the VAMPhyRE project page and extract to bin/vamphyre/
mkdir -p bin/vamphyre
tar -xzf VAMPhyRE_linux_x86_64.tar.gz -C bin/vamphyre --strip-components=1
chmod +x bin/vamphyre/VH5cmdl bin/vamphyre/MergeVGF bin/vamphyre/VFAT
```

---

## Citation

If you use BioRemmer in your research, please cite:

> Cano-Sánchez J, Maldonado-Rodríguez R, Méndez-Tenorio A, Díaz-Ocampo E,
> Larios-Serrato V. (2025). BioRemmer: pipeline for the identification of a
> functional profile for plastic microbial biodegradation.
> *Journal of Bioengineering and Biomedicine Research.*

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Violeta Larios-Serrato — viosdatafactory@gmail.com  
Laboratorio de Biotecnología y Bioinformática Genómica  
Escuela Nacional de Ciencias Biológicas, IPN — Ciudad de México, México
