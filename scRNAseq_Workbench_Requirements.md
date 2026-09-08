# Software Requirements Document
## Downstream Single-Cell Transcriptomics Analysis Workbench

**Document version:** 1.0
**Date:** 2026-09-08
**Prepared for:** Software Engineer (Infrastructure/Tooling)
**Prepared by:** [Project Owner]

---

## 1. Background & Objective

The research team currently performs downstream analysis of single-cell RNA-seq
(scRNA-seq) data using Python (numpy, pandas, scikit-learn, PyDESeq2) inside
Jupyter notebooks. There is currently no standardized project structure,
environment management, or reproducibility tooling.

**Objective:** Build a Jupyter-based workbench project that makes it easy and
reliable to use the standard Python scientific/bioinformatics stack for
scRNA-seq downstream analysis — **without wrapping existing libraries behind
a custom API**. The deliverable is infrastructure and tooling that sits
*around* the analysis code, not a replacement for it.

## 2. Scope

### In scope
- Project scaffolding and folder structure
- Environment/dependency management setup
- Notebook version-control hygiene (git-friendly notebooks)
- Data version control for large intermediate files (`.h5ad`, `.loom`, etc.)
- Reproducible/parameterized notebook execution
- Reusable helper module (`src/`) for shared utility functions (plotting
  themes, QC config, I/O helpers) — imported into notebooks, not wrapping
  analysis calls
- CI checks for code quality (lint/format on `src/`, notebook output
  stripping)
- Report generation from executed notebooks
- Documentation of the project setup for onboarding new analysts

### Out of scope
- Building a custom analysis API/wrapper around scanpy, PyDESeq2, etc.
- A web UI or dashboard (analysts work directly in Jupyter)
- Managing wet-lab / sample metadata systems
- Cluster/HPC job scheduling (unless explicitly added in a later phase)

## 3. Functional Requirements

| # | Requirement | Detail |
|---|---|---|
| FR-1 | Project scaffold | Standard folder layout: `data/raw/`, `data/processed/`, `notebooks/`, `src/`, `reports/`, `configs/`. Provide as a cookiecutter template or setup script. |
| FR-2 | Environment management | `pyproject.toml` (or `environment.yml`) with pinned versions for numpy, pandas, scikit-learn, pydeseq2, scanpy, anndata, and any approved additions. Must support `mamba`/`conda` install in one command. |
| FR-3 | Notebook–git integration | Configure `jupytext` to pair every `.ipynb` with a `.py:percent` file. Configure `nbstripout` as a pre-commit hook so outputs never get committed raw. |
| FR-4 | Data version control | Set up `dvc` (or `git-lfs` if simpler) to track large data files outside git, with a documented workflow for pulling/pushing data. |
| FR-5 | Config-driven parameters | QC thresholds, clustering resolution, DE cutoffs, etc. should live in a `configs/*.yaml` file, not hardcoded in notebooks. |
| FR-6 | Shared utility module | `src/workbench_utils/` package containing plotting themes, QC helper functions, and common I/O wrappers (e.g., load/save AnnData with standard naming). Must be pip-installable in editable mode (`pip install -e .`) so notebooks can `import workbench_utils`. |
| FR-7 | Parameterized execution | Support running any notebook via `papermill` with a parameters cell, so the same notebook can be re-run against different samples/configs without duplication. |
| FR-8 | Report export | Provide a script/make-target that converts an executed notebook into an HTML report and saves it to `reports/`. |
| FR-9 | Onboarding documentation | A `README.md` (or docs site) explaining: environment setup, folder conventions, how to add a new notebook, how to run the pipeline end-to-end. |
| FR-10 | Sample/template notebooks | Provide 2–3 template notebooks demonstrating the intended pattern (e.g., QC → clustering → PyDESeq2 DE) so analysts have a reference to copy from. |

## 4. Non-Functional Requirements

| # | Requirement | Detail |
|---|---|---|
| NFR-1 | Reproducibility | Given the same config and data version, results must be reproducible (fixed random seeds documented in config). |
| NFR-2 | Low friction | Setting up the environment from a clean machine should take ≤ 2 commands (e.g., `mamba env create -f environment.yml && pip install -e .`). |
| NFR-3 | Scalability (data size) | Utility functions should support AnnData `backed='r'` mode for datasets that don't fit in RAM; document when to switch. |
| NFR-4 | Maintainability | `src/` code must have type hints, docstrings, and pass `ruff`/`black` checks in CI. |
| NFR-5 | No vendor lock-in / no wrapper | Analysts must always be able to call scanpy/PyDESeq2/etc. APIs directly; the tooling layer must never intercept or modify library behavior. |
| NFR-6 | Version control cleanliness | Git history must not accumulate notebook output diffs or large binary data files. |

## 5. Proposed Tech Stack

- **Language:** Python 3.11+
- **Environment:** mamba/conda + `pyproject.toml`
- **Notebook tooling:** JupyterLab, `jupytext`, `nbstripout`, `papermill`, `nbconvert`
- **Data versioning:** `dvc` (or `git-lfs`)
- **Analysis libraries (already approved):** numpy, pandas, scikit-learn, pydeseq2, scanpy, anndata
- **Code quality:** `ruff`, `black`, `pre-commit`
- **Testing:** `pytest` (for `src/` utility functions only)
- **CI:** GitHub Actions (lint + test on PR)

## 6. Deliverables

1. Git repository with the scaffolded project structure (FR-1)
2. Working environment file(s) and setup instructions (FR-2)
3. Pre-commit configuration (jupytext + nbstripout + ruff/black) (FR-3)
4. DVC/git-lfs configuration and documented data workflow (FR-4)
5. `src/workbench_utils` installable package (FR-6)
6. 2–3 template notebooks with example config files (FR-5, FR-10)
7. Makefile or CLI scripts for: env setup, notebook execution via papermill, report export (FR-7, FR-8)
8. README / onboarding doc (FR-9)
9. CI pipeline (lint + test) on pull requests

## 7. Acceptance Criteria

- A new team member can clone the repo, run the setup command, and execute
  the template notebook end-to-end without manual troubleshooting.
- `git diff` on a notebook after editing shows meaningful code changes only
  (no cell-output noise).
- Running the same notebook twice with the same config produces identical
  results.
- CI fails the PR if `src/` code fails lint or tests.
- Data files larger than [X MB — to be defined] are never committed directly
  to git; they are tracked via DVC/LFS.

## 8. Open Questions for Engineering Kickoff

- Confirm whether GPU-accelerated tooling (e.g., `rapids-singlecell`) is
  needed now or deferred to a later phase.
- Confirm whether experiment tracking (MLflow/W&B) is in scope for v1.
- Confirm data size expectations (drives the backed-mode / DVC storage
  backend decision — local disk vs. S3/GCS).
- Confirm CI/CD platform constraints (GitHub Actions vs. internal CI).
