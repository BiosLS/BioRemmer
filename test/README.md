# BioRemmer — Test datasets

Three paired-end FASTQ datasets for pipeline validation.
Files are tracked with **Git LFS** — run `git lfs pull` after cloning.

## Datasets

| Files | Purpose |
|-------|---------|
| `Test_short_1/2.fastq.gz`  | Quick smoke test (~2 min, 4 threads) |
| `Test_medium_1/2.fastq.gz` | Intermediate validation (~15 min) |
| `Test_1/2.fastq.gz`        | Full test — same data as the article |

Source: Biofilm metagenomes from plastic particles, North Pacific Gyre
(NCBI SRA: SRP151008)

## Usage

```bash
conda activate bioremmer

# Quick test
./biorem_pipeline_v2.sh test/Test_short_1.fastq.gz test/Test_short_2.fastq.gz test_short 4

# Full test
./biorem_pipeline_v2.sh test/Test_1.fastq.gz test/Test_2.fastq.gz test_full 4
```

Results → `Results/Biorem_report.html`
