# Software Requirements Document
## Downstream Single-Cell Transcriptomics Analysis Workbench (Hybrid R/Python)

**Document version:** 2.0 — revised to adopt hybrid R/Python design on a single Jupyter notebook
**Date:** 2026-09-09
**Prepared for:** Software Engineer (Infrastructure/Tooling)
**Prepared by:** [Project Owner]

---

## 1. Background & Objective

The research team currently performs downstream analysis of single-cell RNA-seq
(scRNA-seq) data using Python (numpy, pandas, scikit-learn, PyDESeq2) inside
Jupyter notebooks. There is currently no standardized project structure,
environment management, or reproducibility tooling.

**Objective:** Build a Jupyter-based workbench project that makes it easy and
reliable to use the standard scRNA-seq stack for downstream analysis —
**without wrapping existing libraries behind a custom API**, and **using the
original reference implementation for each statistical method, whether that
implementation lives in Python or R.** The deliverable is infrastructure and
tooling that sits *around* the analysis code, not a replacement for it.

### 0. Language Selection Principle (new in v2)

| Rule | Use R | Use Python |
|---|---|---|
| An R package is the reference implementation of the statistical method (e.g., differential expression via negative-binomial GLM) | ✅ | — |
| The task is core single-cell processing (QC, normalization, clustering, embedding) where scanpy is the field-standard, actively-maintained, GPU-capable tool | — | ✅ |
| The task benefits from Python-first tooling (deep-learning-based integration, GPU acceleration) | — | ✅ |
| Infra/glue/config/reporting/notebook tooling | — | ✅ |

The notebook kernel stays **Python (Jupyter)**. R is invoked only through
`rpy2`, so analysts never need to open RStudio or write R by hand — but the
statistical engine underneath differential expression calls is the original R
implementation, not a Python port.

This is a deliberate change from v1: PyDESeq2 (a Python re-implementation of
DESeq2) is replaced as the *default* DE engine by a thin `rpy2` wrapper that
calls R's `DESeq2` directly, since DESeq2 is the reference implementation the
field cites. PyDESeq2 remains available as an explicit, documented
alternative for teams that prefer a pure-Python dependency chain (see FR-11).

## 2. Scope

### In scope
- Project scaffolding and folder structure
- Environment/dependency management setup (Python + R/Bioconductor in one environment)
- Notebook version-control hygiene (git-friendly notebooks)
- Data version control for large intermediate files (`.h5ad`, `.loom`, etc.)
- Reproducible/parameterized notebook execution
- Reusable helper module (`src/`) for shared utility functions (plotting
  themes, QC config, I/O helpers) — imported into notebooks, not wrapping
  analysis calls
- A thin, isolated R interop layer (`rpy2`) for statistics whose reference
  implementation is R-only (DESeq2, edgeR, limma, MAST-in-R if needed)
- CI checks for code quality (lint/format on `src/`, notebook output
  stripping)
- Report generation from executed notebooks
- Documentation of the project setup for onboarding new analysts

### Out of scope
- Building a custom analysis API/wrapper around scanpy, DESeq2, etc.
- A web UI or dashboard (analysts work directly in Jupyter)
- Managing wet-lab / sample metadata systems
- Cluster/HPC job scheduling (unless explicitly added in a later phase)
- GPU-accelerated tooling (`rapids-singlecell`) — deferred; R has no
  comparable GPU-accelerated ecosystem, so this stays a pure-Python question
  (see Open Questions)

## 3. Module-by-Module Language Assignment

| Stage | Recommended language | Rationale |
|---|---|---|
| Data loading / AnnData I/O | Python (`anndata`, `scanpy`) | Field-standard object model; no R equivalent needed |
| QC metrics & filtering | Python (`scanpy`) | Mature, actively maintained, GPU-ready if needed later |
| Normalization | Python (`scanpy`) | Standard, well-validated |
| Doublet detection | Python (`scrublet`/`scanpy`) | No compelling R alternative |
| Dimensionality reduction (PCA/UMAP) | Python (`scanpy`) | Standard |
| Clustering (Leiden/Louvain) | Python (`scanpy`, `leidenalg`) | Standard |
| **Differential expression (pseudobulk)** | **R (`DESeq2` via rpy2)** | DESeq2 is the reference implementation; avoids relying on a Python port |
| **Differential expression (alternative: limma-voom)** | **R (`limma`)** | No Python equivalent exists |
| Marker gene detection (Wilcoxon rank-sum, within scanpy) | Python (`scanpy.tl.rank_genes_groups`) | Fine as-is for exploratory marker detection; not a substitute for pseudobulk DE |
| Cell type annotation | Python (`scanpy`, reference-based tools) | Field-standard |
| Trajectory/pseudotime (if added later) | Python (`scanpy`, `scVelo`) | Python-first ecosystem |
| Plotting themes / QC plots | Python (`matplotlib`/`scanpy.pl`) | Presentation layer only |

## 4. Functional Requirements

| # | Requirement | Detail |
|---|---|---|
| FR-1 | Project scaffold | Standard folder layout: `data/raw/`, `data/processed/`, `notebooks/`, `src/`, `reports/`, `configs/`. Provide as a cookiecutter template or setup script. |
| FR-2 | Environment management | A single `environment.yml` (mamba/conda) pinning versions for numpy, pandas, scikit-learn, scanpy, anndata, **plus R and the Bioconductor packages listed in Section 6 (DESeq2, limma, edgeR)**, so Python and R dependencies are resolved together in one command. |
| FR-3 | Notebook–git integration | Configure `jupytext` to pair every `.ipynb` with a `.py:percent` file. Configure `nbstripout` as a pre-commit hook so outputs never get committed raw. |
| FR-4 | Data version control | Set up `dvc` (or `git-lfs` if simpler) to track large data files outside git, with a documented workflow for pulling/pushing data. |
| FR-5 | Config-driven parameters | QC thresholds, clustering resolution, DE cutoffs, **and DE engine choice (`deseq2_r` vs `pydeseq2`)** should live in a `configs/*.yaml` file, not hardcoded in notebooks. |
| FR-6 | Shared utility module | `src/workbench_utils/` package containing plotting themes, QC helper functions, common I/O wrappers, **and the `rpy2` DE interop layer (`workbench_utils.de.deseq2_r()`)**. Must be pip-installable in editable mode (`pip install -e .`). |
| FR-7 | Parameterized execution | Support running any notebook via `papermill` with a parameters cell, so the same notebook can be re-run against different samples/configs without duplication. |
| FR-8 | Report export | Provide a script/make-target that converts an executed notebook into an HTML report and saves it to `reports/`. |
| FR-9 | Onboarding documentation | A `README.md` (or docs site) explaining: environment setup, folder conventions, how to add a new notebook, how to run the pipeline end-to-end, **and how the R interop layer works (so analysts never need to write R directly, but can read the thin wrapper if they want to verify it).** |
| FR-10 | Sample/template notebooks | Provide 2–3 template notebooks demonstrating the intended pattern (e.g., QC → clustering → pseudobulk DE via `DESeq2`) so analysts have a reference to copy from. |
| FR-11 | R↔Python DE interop layer (new) | A function that: (1) accepts a `pandas.DataFrame` pseudobulk count matrix + sample metadata, (2) converts to R objects internally via `rpy2`, (3) runs `DESeqDataSetFromMatrix` → `DESeq` → `results`, (4) converts the result back to a `pandas.DataFrame`. Analysts never touch R syntax. `pydeseq2` remains available as a documented, explicit alternative for teams that want to avoid the R dependency entirely — but it must never be silently substituted for the R engine. |
| FR-12 | R dependency isolation (new) | The `rpy2` + Bioconductor dependency chain must be an optional extra (`pip install workbench_utils[r-stats]` or a separate conda environment file), so the core Python-only workflow (scanpy, QC, clustering) never fails to install if R is unavailable. |
| FR-13 | R interop error handling (new) | If a required R/Bioconductor package is missing, the wrapper must raise a clear, actionable error message (with install instructions), never a raw `rpy2` traceback. |

## 5. Non-Functional Requirements

| # | Requirement | Detail |
|---|---|---|
| NFR-1 | Reproducibility | Given the same config and data version, results must be reproducible (fixed random seeds documented in config, applies to both Python and R RNG where relevant). |
| NFR-2 | Low friction | Setting up the environment from a clean machine should take ≤ 2 commands (e.g., `mamba env create -f environment.yml && pip install -e .`), **including R/Bioconductor.** |
| NFR-3 | Scalability (data size) | Utility functions should support AnnData `backed='r'` mode for datasets that don't fit in RAM; document when to switch. Pseudobulk aggregation before the R DE step keeps the `rpy2` payload small regardless of raw cell count. |
| NFR-4 | Maintainability | `src/` Python code must have type hints, docstrings, and pass `ruff`/`black` checks in CI. The R glue script(s) called via `rpy2` should be kept minimal (thin function calls only, no complex control flow) and version-pinned. |
| NFR-5 | No vendor lock-in / no wrapper | Analysts must always be able to call scanpy/DESeq2/etc. APIs directly; the tooling layer must never intercept or modify library behavior — this applies equally to the R interop layer, which is a pass-through, not an abstraction. |
| NFR-6 | Version control cleanliness | Git history must not accumulate notebook output diffs or large binary data files. |
| NFR-7 | R/Bioconductor version pinning (new) | Pin exact versions of `DESeq2`, `limma`, `edgeR` (and their R version) in `environment.yml` to prevent results drifting across environments or over time. |

## 6. Proposed Tech Stack

- **Language:** Python 3.11+ (kernel), R 4.x (invoked via `rpy2`, never as the notebook kernel)
- **Environment:** mamba/conda + `pyproject.toml` + one shared `environment.yml`
- **Notebook tooling:** JupyterLab, `jupytext`, `nbstripout`, `papermill`, `nbconvert`
- **Data versioning:** `dvc` (or `git-lfs`)
- **Python analysis libraries (already approved):** numpy, pandas, scikit-learn, scanpy, anndata, `rpy2`
- **Python DE (optional, explicit alternative to R engine):** pydeseq2
- **R/Bioconductor (optional extra, via rpy2):** `DESeq2`, `limma`, `edgeR`
- **Code quality:** `ruff`, `black`, `pre-commit`
- **Testing:** `pytest` (for `src/` utility functions), plus regression tests for the `rpy2` DE wrapper against a bare R script run
- **CI:** GitHub Actions (lint + test on PR; R interop tests run in a job with Bioconductor pre-installed)

## 7. Deliverables

1. Git repository with the scaffolded project structure (FR-1)
2. Working environment file(s) covering Python **and R/Bioconductor**, and setup instructions (FR-2)
3. Pre-commit configuration (jupytext + nbstripout + ruff/black) (FR-3)
4. DVC/git-lfs configuration and documented data workflow (FR-4)
5. `src/workbench_utils` installable package, including the R interop layer as an optional extra (FR-6, FR-11, FR-12)
6. 2–3 template notebooks with example config files, showing both the R-backed and Python-only DE path (FR-5, FR-10)
7. Makefile or CLI scripts for: env setup, notebook execution via papermill, report export (FR-7, FR-8)
8. README / onboarding doc, including an explanation of the R interop layer (FR-9)
9. CI pipeline (lint + test, including R interop regression test) on pull requests

## 8. Acceptance Criteria

- A new team member can clone the repo, run the setup command, and execute
  the template notebook end-to-end without manual troubleshooting — including
  the R-backed DE step, without writing any R code themselves.
- `git diff` on a notebook after editing shows meaningful code changes only
  (no cell-output noise).
- Running the same notebook twice with the same config produces identical
  results, for both the Python-only and R-backed paths.
- CI fails the PR if `src/` code fails lint or tests, or if the `rpy2`
  wrapper's output diverges from a direct R script run beyond a documented
  numerical tolerance.
- Data files larger than [X MB — to be defined] are never committed directly
  to git; they are tracked via DVC/LFS.
- The Python-only core workflow (QC → clustering) installs and runs
  successfully even in an environment where R is not installed (R extra is
  truly optional per FR-12).

## 9. Open Questions for Engineering Kickoff

- Confirm whether GPU-accelerated tooling (e.g., `rapids-singlecell`) is
  needed now or deferred to a later phase. (No R equivalent exists, so this
  stays Python-only regardless of the hybrid decision.)
- Confirm whether experiment tracking (MLflow/W&B) is in scope for v1.
- Confirm data size expectations (drives the backed-mode / DVC storage
  backend decision — local disk vs. S3/GCS).
- Confirm CI/CD platform constraints (GitHub Actions vs. internal CI), and
  whether the CI runner image can host both Python and Bioconductor without
  excessive build time.
- Confirm whether `pydeseq2` should remain a first-class supported path
  (FR-11) or be deprecated in favor of the R engine once the interop layer
  is validated.
