# %% [markdown]
# # 🔬 Template 1: Quality Control & Filtering
# **Workflow:** Data Ingestion ➔ QC Metrics ➔ Doublet Detection ➔ Filtering ➔ H5AD Export
#
# Part of the Single-Cell Transcriptomics Analysis Workbench (FR-10, NFR-5).
# Calls standard **Scanpy**, **Scrublet**, and **workbench_utils** helper functions.

# %%
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from workbench_utils.config import load_config
from workbench_utils.io import load_upstream_matrix, save_h5ad
from workbench_utils.plotting import plot_qc_violins, set_publication_style
from workbench_utils.qc import calculate_qc_metrics, detect_doublets, filter_cells, filter_genes

set_publication_style()

# %%
# ==============================================================================
# Parameters Cell for Papermill (FR-7)
# Override via CLI: papermill -p config_path "configs/strict_qc_analysis.yaml"
# ==============================================================================
config_path = "configs/default_analysis.yaml"
input_path = None
output_h5ad = None

# %%
# Load and validate YAML configuration (FR-5, NFR-1)
cfg = load_config(config_path if Path(config_path).exists() else None)
np.random.seed(cfg["project"]["random_seed"])
print(f"Loaded Pipeline Configuration: {cfg['project']['name']}")

# %% [markdown]
# ## 1. Data Ingestion
# Load upstream count matrix or generate reproducible synthetic dataset.

# %%
resolved_input = input_path or cfg["data"].get("input_path")

if resolved_input and Path(resolved_input).exists():
    adata = load_upstream_matrix(resolved_input)
else:
    print("Upstream matrix not found on disk. Generating synthetic demonstration dataset...")
    n_cells = 180
    n_genes = 300
    genes = (
        [f"GENE_{i:03d}" for i in range(260)]
        + [f"MT-ND{j}" for j in range(20)]
        + [f"RPS{k}" for k in range(20)]
    )
    cells = [f"cell_{i:03d}" for i in range(n_cells)]

    # Generate count matrix with biological condition
    rng = np.random.default_rng(cfg["project"]["random_seed"])
    counts = rng.poisson(lam=12.0, size=(n_cells, n_genes)).astype(np.float32)

    # Add higher counts to first 10 genes in Group A
    counts[:90, :15] += rng.poisson(lam=25.0, size=(90, 15))

    obs = pd.DataFrame(index=cells)
    obs["condition"] = ["Treated"] * 90 + ["Control"] * 90
    obs["sample_id"] = [f"sample_{(i % 6) + 1}" for i in range(n_cells)]

    var = pd.DataFrame(index=genes)
    var["gene_name"] = genes

    adata = ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()

print(f"Initial Dataset: {adata.n_obs} cells × {adata.n_vars} genes")

# %% [markdown]
# ## 2. Calculate QC Metrics & Visualize Violins
# Compute UMI counts, detected genes, and mitochondrial/ribosomal fractions.

# %%
adata = calculate_qc_metrics(
    adata,
    mito_prefix=tuple(cfg["qc"]["mito_prefix"]),
    ribo_prefix=tuple(cfg["qc"]["ribo_prefix"]),
)

fig = plot_qc_violins(
    adata,
    keys=("n_genes_by_counts", "total_counts", "pct_counts_mito", "pct_counts_ribo"),
)

# %% [markdown]
# ## 3. Scrublet Doublet Detection
# Estimate multiplet contamination using simulation on the raw count matrix.

# %%
if cfg["qc"].get("run_doublet_detection", True):
    adata = detect_doublets(
        adata,
        expected_doublet_rate=cfg["qc"]["expected_doublet_rate"],
        random_state=cfg["project"]["random_seed"],
    )

# %% [markdown]
# ## 4. Filter Cells and Genes
# Filter low-depth droplets and dead cells according to config thresholds.

# %%
filtered_adata = filter_cells(
    adata,
    min_genes=cfg["qc"]["min_genes"],
    max_genes=cfg["qc"].get("max_genes"),
    min_counts=cfg["qc"]["min_counts"],
    max_counts=cfg["qc"].get("max_counts"),
    max_pct_mito=cfg["qc"]["max_pct_mito"],
    max_doublet_score=cfg["qc"].get("max_doublet_score"),
)

filtered_adata = filter_genes(
    filtered_adata,
    min_cells=cfg["filter_genes"]["min_cells"],
)

print(f"Post-Filtering Dataset: {filtered_adata.n_obs} cells × {filtered_adata.n_vars} genes")

# %% [markdown]
# ## 5. Save Filtered AnnData
# Export clean dataset to `data/processed/` for downstream clustering.

# %%
out_target = output_h5ad or "data/processed/01_qc_filtered.h5ad"
save_h5ad(filtered_adata, out_target)
print(f"QC & Filtering step completed. Output saved to: {out_target}")
