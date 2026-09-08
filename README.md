# Downstream Single-Cell Transcriptomics Analysis Workbench

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)

A reproducible, production-grade downstream analysis workbench for single-cell RNA-seq (scRNA-seq) data in Python, built around **Scanpy**, **PyDESeq2**, and **AnnData**.

Complies with the technical and scientific specifications in [scRNAseq_Workbench_Requirements.md](file:///home/surj/Workspace/SC_RNA_workbench_2026/scRNAseq_Workbench_Requirements.md).

---

## 1. Architectural Philosophy (NFR-5)

> **No Vendor Lock-in / No Custom API Wrapper:**
> The workbench provides infrastructure and tooling that sits **around** your analysis code — not a replacement for it. Analysts write native Scanpy (`sc.pp.*`, `sc.tl.*`, `sc.pl.*`) and PyDESeq2 code directly in Jupyter Notebooks. The `workbench_utils` module supplies non-invasive utilities for plotting themes, QC calculation, disk-backed AnnData loading, and config-driven parameters.

---

## 2. Directory Conventions (FR-1)

```text
SC_RNA_workbench_2026/
├── configs/               # YAML configuration files (FR-5, NFR-1)
│   ├── default_analysis.yaml
│   ├── strict_qc_analysis.yaml
│   └── pseudobulk_de_analysis.yaml
├── data/
│   ├── raw/               # Raw upstream matrices / 10x folders (Read-only)
│   └── processed/         # Normalized/clustered .h5ad files (Git LFS tracked)
├── docs/                  # Architectural & data versioning documentation
│   └── data_versioning.md
├── notebooks/             # Paired notebooks (Jupytext .ipynb + .py:percent)
│   ├── 01_qc_and_filtering.ipynb
│   ├── 02_clustering_and_annotation.ipynb
│   └── 03_differential_expression_pydeseq2.ipynb
├── reports/               # Executed notebook HTML reports & marker tables (FR-8)
│   ├── 01_qc_and_filtering.html
│   ├── 02_clustering_and_annotation.html
│   ├── 03_differential_expression_pydeseq2.html
│   └── significant_differential_expression_results.csv
├── scripts/               # Automation and scaffolding scripts
│   ├── init_workbench.sh
│   └── run_notebook.sh
├── src/
│   └── workbench_utils/   # Shared utility module (FR-6, NFR-4)
│       ├── __init__.py
│       ├── config.py      # YAML config loader & validation
│       ├── de.py          # PyDESeq2 data preparation & execution
│       ├── io.py          # Backed-mode AnnData I/O & Seurat export (NFR-3)
│       ├── plotting.py    # Nature/Cell themes, violins & volcano plots
│       └── qc.py          # QC metrics, filtering & Scrublet doublets
├── tests/                 # Unit test suite (pytest)
├── .gitattributes         # Git LFS data versioning tracking (FR-4)
├── .pre-commit-config.yaml# Git hygiene hooks: nbstripout + jupytext + ruff (FR-3)
├── environment.yml        # Conda / Mamba environment definition (FR-2)
├── pixi.toml              # Pixi reproducible environment definition
├── pyproject.toml         # Python package metadata & dependencies
└── Makefile               # CLI automation targets (FR-7, FR-8)
```

---

## 3. Quickstart & Installation (NFR-2)

Setting up from a clean machine requires **≤ 2 commands**:

### Option A: Using Pixi (Recommended)
```bash
pixi install
make setup
```

### Option B: Using Mamba / Conda
```bash
mamba env create -f environment.yml
conda activate sc_workbench
pip install -e .
```

---

## 4. End-to-End Analysis Workflow (FR-10)

The workbench includes 3 reproducible template notebooks demonstrating the standard single-cell analysis workflow:

### Step 1: QC & Doublet Filtering
- **Notebook:** [notebooks/01_qc_and_filtering.ipynb](file:///home/surj/Workspace/SC_RNA_workbench_2026/notebooks/01_qc_and_filtering.ipynb)
- **Actions:** Computes library depth, detects MT/ribosomal fractions, simulates doublets with Scrublet, and applies threshold filtering.
- **Run & Export Report:**
  ```bash
  make report NB=01_qc_and_filtering CONFIG=configs/default_analysis.yaml
  ```
- **Output:** `data/processed/01_qc_filtered.h5ad` and `reports/01_qc_and_filtering.html`

### Step 2: Normalization, Clustering & UMAP
- **Notebook:** [notebooks/02_clustering_and_annotation.ipynb](file:///home/surj/Workspace/SC_RNA_workbench_2026/notebooks/02_clustering_and_annotation.ipynb)
- **Actions:** Library scaling (`target_sum=10k`), log1p transformation, HVG selection, PCA, k-NN graph, UMAP projection, Leiden community detection, and cluster biomarker ranking.
- **Run & Export Report:**
  ```bash
  make report NB=02_clustering_and_annotation CONFIG=configs/default_analysis.yaml
  ```
- **Output:** `data/processed/02_clustered.h5ad` and `reports/02_clustering_and_annotation.html`

### Step 3: Differential Expression (PyDESeq2)
- **Notebook:** [notebooks/03_differential_expression_pydeseq2.ipynb](file:///home/surj/Workspace/SC_RNA_workbench_2026/notebooks/03_differential_expression_pydeseq2.ipynb)
- **Actions:** Formats raw integer counts, fits Negative Binomial GLM with PyDESeq2, evaluates biological contrasts, plots publication Volcano Plots, and exports significant markers.
- **Run & Export Report:**
  ```bash
  make report NB=03_differential_expression_pydeseq2 CONFIG=configs/default_analysis.yaml
  ```
- **Output:** `reports/03_differential_expression_pydeseq2.html` and `reports/significant_differential_expression_results.csv`

---

## 5. Config-Driven Parameters (FR-5, NFR-1)

All analysis parameters live in [configs/default_analysis.yaml](file:///home/surj/Workspace/SC_RNA_workbench_2026/configs/default_analysis.yaml):
- **Fixed Random Seed:** `project.random_seed: 42` guarantees identical results across runs.
- **QC Thresholds:** `min_genes`, `max_pct_mito`, `expected_doublet_rate`.
- **Dimensionality & Clustering:** `n_pcs: 30`, `resolution: 0.8`.
- **PyDESeq2 Parameters:** `design_factor: "condition"`, `fdr_cutoff: 0.05`, `log2fc_cutoff: 1.0`.

To validate a configuration file before execution:
```bash
make validate-config CONFIG=configs/strict_qc_analysis.yaml
```

---

## 6. Git Hygiene & Data Versioning (FR-3, FR-4, NFR-6)

- **Clean Git Diffs:** `nbstripout` automatically strips raw cell outputs before git commit.
- **Jupytext Pairing:** Every notebook is paired as `ipynb,py:percent`. Code diffs are reviewable line-by-line in git.
- **Git LFS:** Large single-cell files (`*.h5ad`, `*.loom`, `*.mtx.gz`, `*.bam`) are tracked via Git LFS outside the git commit tree. See [docs/data_versioning.md](file:///home/surj/Workspace/SC_RNA_workbench_2026/docs/data_versioning.md) for workflow details.

To install pre-commit hooks locally:
```bash
pre-commit install
```

---

## 7. How to Add a New Notebook

1. Create a new notebook in [notebooks/](file:///home/surj/Workspace/SC_RNA_workbench_2026/notebooks/) (e.g. `04_pathway_enrichment.ipynb`).
2. Pair it with Jupytext:
   ```bash
   pixi run jupytext --set-formats ipynb,py:percent notebooks/04_pathway_enrichment.ipynb
   ```
3. Tag the second cell with `tags=["parameters"]` and define `config_path = "configs/default_analysis.yaml"`.
4. Load config parameters:
   ```python
   from workbench_utils.config import load_config
   cfg = load_config(config_path)
   ```
5. Execute and generate reports:
   ```bash
   make report NB=04_pathway_enrichment
   ```

---

## 8. Testing & Quality Assurance

Run the test suite and verify code quality:

```bash
# Run 29 unit tests covering I/O, QC, PyDESeq2, and Plotting
make test

# Run Ruff linter
make lint

# Run Ruff code formatter
make format
```
