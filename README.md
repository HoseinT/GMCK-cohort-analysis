# GMCK Cohort Analysis Pipeline

Snakemake pipeline for somatic variant calling and copy number analysis in paired primary/relapse tumour samples. Each patient contributes three sample types: germline (lymph node), primary tumour, and relapse tumour, sequenced with a targeted panel using paired-end reads with inline UMIs.

## Pipeline overview

```
FASTQ (per sample)
  └─ UMI extraction (fgbio FastqToBam)
  └─ Adapter trimming (Cutadapt)
  └─ Alignment (BWA mem → GRCh37/hg19)
  └─ Merge aligned + unaligned BAM (GATK MergeBamAlignment)
  └─ UMI-aware duplicate marking (GATK UmiAwareMarkDuplicatesWithMateCigar)
  └─ Base quality score recalibration (GATK BQSR)
  └─ Somatic variant calling (GATK Mutect2, tumour vs germline)
  └─ Annotation (VEP v108 + CADD pre-scored)
  └─ Filtering (depth ≥ 50, VAF > 5%, consequence, CADD ≥ 0.20)
  └─ Copy number analysis (snp-pileup + FACETS)

Cohort (all patients combined)
  └─ Tumour mutational burden (SNVs / 2.4 Mb panel)
  └─ Mutual exclusivity (DISCOVER)
  └─ Ploidy-adjusted CNV burden + recurrence
  └─ Paired Wilcoxon signed-rank test (primary vs relapse, BH-corrected)
  └─ Plotly visualisations
```

## Software versions

### Bioinformatics tools

| Tool | Version | Notes |
|------|---------|-------|
| BWA | 0.7.17 | Alignment to GRCh37/hg19 |
| GATK | 4.2.6.1 | UMI dedup, BQSR, Mutect2, FilterMutectCalls |
| Cutadapt | 4.1 | Adapter and quality trimming |
| fgbio | — | UMI extraction (FastqToBam, MergeShards); verify version with `fgbio --version` |
| samtools | 1.6+ | BAM sorting and indexing |
| VEP | 108 | Variant annotation, GRCh37 cache |
| snp-pileup | — | Allele counts for FACETS; verify version |
| FACETS | 0.6.2 | Copy number analysis (R package) |

### Python packages

| Package | Version used in study |
|---------|----------------------|
| Python | 3.8.12 |
| pandas | 1.5.3 |
| numpy | 1.26.4 |
| scipy | 1.14.0 |
| statsmodels | 0.14.2 |
| plotly | 5.12 |
| pysam | — | Required for CADD lookup; verify version |
| pyranges | — | Required for CNV segment→gene mapping; verify version |
| discover | 0.9.5 | Mutual exclusivity (DISCOVER Python package) |

### R packages

| Package | Version |
|---------|---------|
| R | 3.6.3+ |
| facets | 0.6.2 |
| optparse | — |
| data.table | — |

## Installation

### Python environment

```bash
conda env create -f envs/python.yaml
conda activate gmck-pipeline
```

### Bioinformatics tools

```bash
conda env create -f envs/tools.yaml
conda activate gmck-tools
```

Or load via your HPC module system (e.g. `module load bwa/0.7.17 gatk/4.2.6.1`).

### R packages

```r
install.packages("optparse")
install.packages("data.table")
# FACETS from GitHub:
devtools::install_github("mskcc/facets", ref = "v0.6.2")
```

## Usage

1. Copy `sample_sheet_template.tsv` to `sample_sheet.tsv` and fill in paths.
2. Edit `config/config.yaml` — fill in all `/path/to/...` entries.
3. Dry-run to verify the DAG:
   ```bash
   cd workflow
   snakemake -n --cores 1
   ```
4. Full run:
   ```bash
   snakemake --cores 16
   ```

## Variant filtering criteria

- PASS calls from FilterMutectCalls only
- Read depth ≥ 50 in the tumour sample
- VAF > 5% in either primary **or** relapse (union)
- Consequences excluded: intronic, downstream gene, upstream gene, 3′/5′ UTR, non-coding transcript exon
- CADD PHRED ≥ 0.20
- Reference genome: GRCh37/hg19

## CNV class definitions (FACETS)

| Class | Total copy number |
|-------|------------------|
| Deep loss | TCN = 0 |
| Loss | TCN < 2 |
| CNLOH | TCN = 2, LCN = 0 |
| Neutral | TCN = 2, LCN = 1 |
| Gain | TCN 3–5 |
| Amplification | TCN ≥ 6 |

FACETS run with `--max-depth 20000` to accommodate deep panel sequencing.
