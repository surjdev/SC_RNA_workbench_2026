# Downstream Single-Cell Transcriptomics Analysis Workbench (v2.0)

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)

A reproducible, production-grade downstream analysis workbench for single-cell transcriptomics data (e.g. SMART-seq2 full-length and droplet protocols) using a **Hybrid R/Python Architecture**. Built around **Scanpy**, **AnnData**, reference **R DESeq2** (via an isolated `rpy2` interop layer), and **PyDESeq2** as a pure-Python alternative.

Complies with the technical and scientific specifications in [scRNAseq_Workbench_Requirements_v2.md](file:///home/surj/Workspace/SC_RNA_workbench_2026/scRNAseq_Workbench_Requirements_v2.md).

---

## 1. Architectural Philosophy & Language Selection Principle (NFR-5, FR-11, FR-12)

The workbench follows the **Language Selection Principle**: use the best-in-class tool for each stage rather than forcing an all-Python or all-R pipeline, while maintaining a single, unified Python environment for analysts.

| Analysis Stage | Primary Tool | Runtime Language | Rationale |
|---|---|---|---|
| **AnnData I/O & Slicing** | AnnData / Scanpy | Python | Fast in-memory array manipulation, native `.h5ad` support |
| **Cell QC & Filtering** | Scanpy / NumPy | Python | Vectorized filtering, ERCC spike-in & mito computation |
| **Plate Layout QC** | seaborn / matplotlib | Python | Interactive 96/384-well plate spatial heatmaps |
| **Doublet Detection** | Scrublet | Python | Standard simulation-based doublet scoring |
| **Normalization & HVG** | Scanpy | Python | Native CPM/TPM library size scaling, log1p transformation |
| **Embeddings & Clustering**| Scanpy (PCA, UMAP, Leiden) | Python | Fast C/C++ graph community detection |
| **Pseudobulk DE (Default)**| **DESeq2** | **R (via rpy2)** | **Gold-standard reference method** for Negative Binomial GLMs |
| **Pseudobulk DE (Alt)**    | **limma-voom** | **R (via rpy2)** | Gold-standard linear modeling for large sample cohorts |
| **Pseudobulk DE (Python)** | **PyDESeq2** | **Python** | Pure-Python alternative when R is not installed (FR-11) |
| **Visualization & Reports**| matplotlib / Papermill | Python | Nature/Cell styles, HTML report rendering |

### Key Guarantees:
1. **Zero Raw R for Analysts:** Analysts never need to write R scripts, manage R sessions, or manually marshal data. The `workbench_utils.de` interop layer converts Pandas count matrices and metadata to R objects, runs the statistical pipeline, and returns clean Pandas DataFrames.
2. **No Silent Fallback (FR-11):** If `deseq2_r` is configured and R or required Bioconductor packages are unavailable, the workbench raises `RDependencyError` with actionable commands. It will never silently substitute PyDESeq2 for R DESeq2.
3. **Optional Dependency Isolation (FR-12):** The `r-stats` feature (`rpy2`, Bioconductor) is isolated. Core Python processing (Scanpy, QC, clustering) will never crash if R is absent.
4. **No Monolithic Wrapper (NFR-5):** Tools sit *around* standard libraries. Analysts call native Scanpy functions directly.

---

## 2. Directory Structure (FR-1)

```text
SC_RNA_workbench_2026/
├── configs/                          # YAML analysis configurations (FR-5, NFR-1)
│   ├── default_analysis.yaml         # Standard parameters (CPM, R DESeq2 default)
│   ├── strict_qc_analysis.yaml       # High-stringency QC thresholds
│   ├── pseudobulk_de_analysis.yaml   # Pseudobulk DE with contrast configuration
│   └── smartseq2_plate_analysis.yaml # SMART-seq2 multi-plate analysis with ERCC QC
├── data/
│   ├── raw/                          # Upstream count matrices & plate metadata (Read-only)
│   └── processed/                    # Clustered .h5ad files (Git LFS tracked)
├── docs/                             # User guide, versioning & architectural documentation
│   ├── user_guide.md
│   └── data_versioning.md
├── notebooks/                        # Paired template notebooks (Jupytext .ipynb + .py:percent)
│   ├── 01_qc_and_filtering.ipynb
│   ├── 02_clustering_and_annotation.ipynb
│   └── 03_differential_expression.ipynb  # Dual-engine pseudobulk DE (R DESeq2 & PyDESeq2)
├── reports/                          # Executed HTML reports, CSVs & figures (FR-8)
│   ├── differential_expression_results.csv
│   ├── significant_differential_expression_results.csv
│   └── figures/
├── scripts/                          # Automation, benchmarking and scaffolding
│   ├── init_workbench.sh
│   ├── run_notebook.sh
│   └── run_comprehensive_benchmark.py
├── src/
│   └── workbench_utils/              # Shared utility package (FR-6, NFR-4)
│       ├── __init__.py
│       ├── config.py                 # YAML config loader & validation
│       ├── de.py                     # Hybrid R/Python DE interop layer & PyDESeq2
│       ├── io.py                     # Matrix ingestion & Seurat export
│       ├── plotting.py               # Nature/Cell themes, violins, volcano & plate heatmaps
│       └── qc.py                     # ERCC spike-ins, mito QC, plate parsing & filtering
├── tests/                            # Unit & reproducibility test suite (pytest - 42 tests)
├── .github/workflows/ci.yml          # GitHub Actions CI/CD pipeline
├── .gitattributes                    # Git LFS data versioning tracking (FR-4)
├── .pre-commit-config.yaml           # Git hygiene hooks: nbstripout + jupytext + ruff (FR-3)
├── environment.yml                   # Conda / Mamba environment definition (FR-2)
├── pixi.toml / pixi.lock             # Pixi reproducible binary lockfile
├── pyproject.toml                    # Package metadata & optional dependency extras
└── Makefile                          # CLI automation targets (FR-7, FR-8)
```

---

## 3. Quickstart & Installation (NFR-2, FR-12)

### Option A: Using Pixi (Recommended)
```bash
# Pure-Python default environment:
pixi install
make setup

# With R & Bioconductor statistical packages (r-stats feature):
pixi run -e r-stats pytest
```

### Option B: Using Mamba / Conda
```bash
mamba env create -f environment.yml
conda activate sc_workbench
pip install -e .
```

If R or Bioconductor packages are missing on your host machine, install them using:
```bash
# In R:
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "limma", "edgeR"))
```

---

## 4. End-to-End Analysis Workflow (FR-10)

The workbench provides 3 standard template notebooks:

### Template 1: Ingestion, ERCC Spike-ins & Plate QC
- **Notebook:** [notebooks/01_qc_and_filtering.ipynb](file:///home/surj/Workspace/SC_RNA_workbench_2026/notebooks/01_qc_and_filtering.ipynb)
- **Workflow:** Ingests count matrices and plate coordinates, computes Ambion ERCC spike-in & mito percentages, displays 96/384-well plate spatial heatmaps, filters failed wells, and exports clean AnnData.
- **Run & Export Report:**
  ```bash
  make report NB=01_qc_and_filtering CONFIG=configs/default_analysis.yaml
  ```

### Template 2: Normalization, Clustering & UMAP
- **Notebook:** [notebooks/02_clustering_and_annotation.ipynb](file:///home/surj/Workspace/SC_RNA_workbench_2026/notebooks/02_clustering_and_annotation.ipynb)
- **Workflow:** Excludes synthetic ERCC controls from library scaling, performs CPM/log1p normalization, HVG selection, PCA, Harmony/plate batch correction, UMAP projection, and Leiden community detection.
- **Run & Export Report:**
  ```bash
  make report NB=02_clustering_and_annotation CONFIG=configs/default_analysis.yaml
  ```

### Template 3: Pseudobulk Differential Expression (Hybrid R/Python)
- **Notebook:** [notebooks/03_differential_expression.ipynb](file:///home/surj/Workspace/SC_RNA_workbench_2026/notebooks/03_differential_expression.ipynb)
- **Workflow:** Aggregates single-cell counts per biological sample (`prepare_pseudobulk_data`), executes the reference R `DESeq2` engine (`deseq2_r`), provides explicit pure-Python `deseq2_python` comparison, generates Nature-styled volcano plots, and exports CSV results to `reports/`.
- **Run & Export Report:**
  ```bash
  make report NB=03_differential_expression CONFIG=configs/pseudobulk_de_analysis.yaml
  ```

---

## 5. Config-Driven Differential Expression (FR-5, FR-11)

Differential expression settings are controlled in `configs/*.yaml`:

```yaml
differential_expression:
  engine: "deseq2_r"          # Options: "deseq2_r" (default), "pydeseq2", "limma_voom_r"
  r_seed: 42                  # R random seed for reproducibility
  design_factor: "condition"  # Experimental factor to test
  sample_key: "sample_id"     # Metadata column for biological sample aggregation
  fdr_cutoff: 0.05            # Benjamini-Hochberg FDR threshold
  log2fc_cutoff: 1.0          # Absolute log2 fold change threshold
  contrast:
    - "condition"
    - "Treated"
    - "Control"
```

To switch to pure-Python execution without modifying notebook code, simply set `engine: "pydeseq2"` in your config YAML.

---

## 6. Git Hygiene & Data Versioning (FR-3, FR-4, NFR-6)

- **Clean Git Diffs:** `nbstripout` automatically strips raw cell outputs before git commit.
- **Jupytext Pairing:** Every notebook is paired as `ipynb,py:percent`. Code diffs are reviewable line-by-line in git.
- **Git LFS:** Large single-cell files (`*.h5ad`, `*.loom`, `*.mtx.gz`, `*.bam`) are tracked via Git LFS outside the git commit tree. See [docs/data_versioning.md](file:///home/surj/Workspace/SC_RNA_workbench_2026/docs/data_versioning.md) for workflow details.

Install pre-commit hooks:
```bash
pre-commit install
```

---

## 7. Quality Assurance & Tests

Run all unit and reproducibility tests:

```bash
# Run 42 unit & reproducibility tests covering I/O, QC, R-interop, PyDESeq2, Plotting, and Config
pixi run pytest tests/ -v

# Run code style & formatting checks
pixi run ruff check src/ tests/ notebooks/
pixi run ruff format --check src/ tests/ notebooks/
```
