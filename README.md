# GWAS Analysis Pipeline: NHGRI-EBI Population & Gene-Disease Architecture

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A complete, production-grade, reproducible bioinformatics pipeline in Python for auditing, filtering, analyzing, and visualizing human genome-wide association study (GWAS) data from the **NHGRI-EBI GWAS Catalog**.

---

## Project Objectives

1. **Complete Data Ingestion & Profiling**: Stream and audit the full NHGRI-EBI GWAS Catalog v1.0.2 release (`1,192,032` records, 38 columns) and ancestry metadata (`229,942` studies).
2. **Analysis 1 (Low-Frequency Allele Filtering)**: Isolate variants and genes with Risk Allele Frequency (RAF) between **0.5% and 2.0% inclusive** for **European (EUR)** and **South Asian (SAS)** cohorts.
3. **Analysis 2 (Gene-Disease-Population Landscape)**: Map the unique set of genes associated with diseases/traits, distinguishing **European-only (10,637)**, **South Asian-only (124)**, and **shared (19,293)** disease genes.
4. **Publication-Quality Visualizations**: Generate 300 DPI publication plots (PNG & vector PDF) exploring allele-frequency distributions, trait discovery breakdowns, and ancestral gene overlap.
5. **Reproducibility**: Modular, fully documented, tested, and executable with a single command.

---

## Key Empirical Results

| Analysis Stage | Dimension | Exact Result |
| :--- | :--- | :---: |
| **Catalog Dimensions** | Total Associations Scanned | **1,192,032** |
| | Total Studies Mapped with Ancestry | **229,942** |
| | Total Unique Diseases / Traits | **38,572** |
| **Data Quality Audit** | Numeric Allele Frequency Records | **575,750** (48.3%) |
| | Unreported Allele Frequency (`NR`) | **616,154** (51.7%) |
| **Analysis 1: Low Frequency (0.5%–2.0% RAF)** | Filtered Association Records | **27,872** |
| | European-Associated Records | **27,782** |
| | South Asian-Associated Records | **1,718** |
| | Shared Records (Both EUR & SAS) | **1,628** |
| | Unique Low-Frequency Variants (`rsIDs`) | **19,808** |
| | Unique Low-Frequency Genes | **11,730** |
| **Analysis 2: Gene-Disease Landscape** | Total Unique Disease Genes | **30,054** |
| | European-Only Genes | **10,637** (35.4%) |
| | South Asian-Only Genes | **124** (0.4%) |
| | Shared Genes (Both EUR & SAS) | **19,293** (64.2%) |

---

## Publication Figures

### Figure 1: Allele Frequency Distribution & Cohort Comparison
`figures/fig1_allele_frequency_distribution.png`  
Dual-panel comparison of Risk Allele Frequency across European and South Asian cohorts in the 0.5%–2.0% window via Kernel Density Estimates (KDE) and Box/Stripplots.

### Figure 2: Gene-Disease Discovery Landscape
`figures/fig2_gene_disease_population_landscape.png`  
Stacked horizontal bar chart of top 20 complex traits and proportion heatmap illustrating gene discovery by population context.

### Figure 3: Ancestral Specificity & Overlap
`figures/fig3_population_gene_overlap.png`  
Proportional breakdown of the 30,054 disease-linked genes across European, South Asian, and shared sets.

---

## Repository Structure

```
gwas-analysis-pipeline/
├── data/
│   ├── raw/                             # Downloaded raw datasets (associations & ancestry)
│   └── processed/                       # Processed result tables
│       ├── low_freq_variants_eur_sas.tsv.gz          # Filtered records (0.5%-2.0% RAF)
│       ├── low_freq_genes_and_variants_summary.tsv   # Variant-gene summary pairs
│       ├── gene_disease_population_mapping.tsv.gz    # Gene-disease-population triplets
│       ├── population_gene_classification.tsv        # 30,054 genes categorized
│       └── gwas_executive_summary.csv                # Core summary table for reporting
├── figures/                             # 300 DPI figures (PNG & vector PDF)
│   ├── fig1_allele_frequency_distribution.png / .pdf
│   ├── fig2_gene_disease_population_landscape.png / .pdf
│   └── fig3_population_gene_overlap.png / .pdf
├── gwas_pipeline/                       # Modular Python package
│   ├── __init__.py
│   ├── downloader.py                    # Download manager with progress tracking
│   ├── inspector.py                     # 38-column profiling & missingness auditor
│   ├── filter_variants.py               # Analysis 1 filtering engine
│   ├── gene_disease_mapper.py           # Analysis 2 mapping engine
│   └── visualizer.py                    # Matplotlib/Seaborn visualization suite
├── docs/
│   ├── data_dictionary.md               # 38-column data dictionary
│   └── methodology.md                   # Formal methodology documentation
├── tests/
│   └── test_pipeline.py                 # Automated unit and integration tests
├── GWAS_ANALYSIS_REPORT_FOR_CHATGPT.md  # Comprehensive standalone report for LLMs/ChatGPT
├── requirements.txt                     # Pinned dependencies
├── run_pipeline.py                      # Master CLI pipeline orchestrator
└── README.md                            # Repository guide
```

---


---

## Nextflow Workflow & Shell Script Runner

The entire analysis pipeline has been implemented as a production **Nextflow DSL2** pipeline (`main.nf`), orchestrated via `nextflow.config`, with support for local, Conda, Docker, and HPC/SLURM cluster profiles.

### 1. Unified Shell Script Runner (`run_pipeline.sh`)
The repository includes an executable POSIX shell script (`run_pipeline.sh`) that automatically detects the runtime environment and executes the pipeline:

```bash
# Make script executable
chmod +x run_pipeline.sh

# Display help and options
./run_pipeline.sh --help

# Run the pipeline (auto-detects Nextflow or runs Python)
./run_pipeline.sh

# Run specifically via Nextflow with Conda or Docker
./run_pipeline.sh --nextflow --profile conda
./run_pipeline.sh --nextflow --profile docker

# Run directly via Python / UV
./run_pipeline.sh --python
```

### 2. Running Directly with Nextflow
```bash
# Local standard run
nextflow run main.nf -profile standard --outdir results

# Conda environment run
nextflow run main.nf -profile conda --outdir results

# Docker container run
nextflow run main.nf -profile docker --outdir results

# HPC Cluster (SLURM) run
nextflow run main.nf -profile slurm --outdir results
```

### Nextflow Pipeline Processes:
- **`DOWNLOAD_DATA`**: Fetches the official NHGRI-EBI GWAS Catalog associations and ancestry datasets if absent.
- **`INSPECT_CATALOG`**: Generates a 38-column quality and missingness report (`catalog_inspection_report.txt`).
- **`FILTER_LOW_FREQUENCY`**: Filters for 0.5%–2.0% RAF across European and South Asian cohorts, emitting filtered tables.
- **`MAP_GENE_DISEASE`**: Extracts 30,054 disease-associated genes, mapping them to European-only, South Asian-only, and shared categories.
- **`GENERATE_VISUALIZATIONS`**: Produces 300 DPI publication plots (PNG and PDF) for allele frequency distributions and trait landscapes.
- **`GENERATE_REPORT`**: Synthesizes the executive CSV and workflow completion summaries.

## Quickstart & Reproduction

### Prerequisites
- Python 3.10+
- `pip` or `uv`

### Installation
```bash
git clone https://github.com/your-username/gwas-analysis-pipeline.git
cd gwas-analysis-pipeline
pip install -r requirements.txt
```

### Running the Complete Pipeline
```bash
# Execute all steps end-to-end:
python run_pipeline.py --all
```

Or run specific stages individually:
```bash
python run_pipeline.py --download      # Step 1: Download raw data
python run_pipeline.py --inspect       # Step 2: Inspect schema & missingness
python run_pipeline.py --filter        # Step 3: Analysis 1 low-frequency filtering
python run_pipeline.py --map-genes     # Step 4: Analysis 2 gene-disease mapping
python run_pipeline.py --plot          # Step 5: Render publication figures
```

### Running Unit Tests
```bash
pytest tests/test_pipeline.py -v
```

---

## Citation & Attribution
Data sourced from the **NHGRI-EBI GWAS Catalog**:  
Sollis E, et al. *The NHGRI-EBI GWAS Catalog: knowledgebase and summary statistics.* Nucleic Acids Res. 2023 Jan 6;51(D1):D977-D985.
