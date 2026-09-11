# %% [markdown]
# # 🔬 Template 1: Quality Control & Filtering (SMART-seq2)
# **Workflow:** Data Ingestion ➔ Plate Coordinates ➔ QC Metrics (Mito & ERCC Spike-ins) ➔ Plate Layout Heatmap ➔ Filtering ➔ H5AD Export
#
# Part of the Single-Cell SMART-seq2 Transcriptomics Analysis Workbench (FR-10, NFR-5).
# Calls standard **Scanpy** and **workbench_utils** helper functions.

# %%
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from workbench_utils.config import load_config
from workbench_utils.io import load_smartseq2_matrix, load_upstream_matrix, save_h5ad
from workbench_utils.plotting import plot_plate_layout, plot_qc_violins, set_publication_style
from workbench_utils.qc import (
    calculate_qc_metrics,
    detect_doublets,
    filter_cells,
    filter_genes,
    parse_plate_metadata,
)

set_publication_style()

# %%
# ==============================================================================
# Parameters Cell for Papermill (FR-7)
# Override via CLI: papermill -p config_path "configs/smartseq2_plate_analysis.yaml"
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
# ## 1. Data Ingestion (SMART-seq2 Matrix & Plate Metadata)
# Load upstream SMART-seq2 count matrix (TSV/CSV/MTX) or generate reproducible synthetic 96-well plate dataset.

# %%
resolved_input = input_path or cfg["data"].get("input_path")

if resolved_input and Path(resolved_input).exists():
    metadata_file = cfg["data"].get("metadata_path")
    if metadata_file and Path(metadata_file).exists():
        adata = load_smartseq2_matrix(resolved_input, metadata_path=metadata_file)
    else:
        adata = load_upstream_matrix(resolved_input)
else:
    print(
        "Upstream matrix not found on disk. Generating synthetic SMART-seq2 96-well plate demonstration dataset..."
    )
    rng = np.random.default_rng(cfg["project"]["random_seed"])
    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    wells = [f"{r}{c:02d}" for r in rows for c in range(1, 13)]
    cell_names = [f"Plate1_{w}" for w in wells]
    n_cells = len(cell_names)  # 96 wells

    ercc_genes = [f"ERCC-{i:05d}" for i in range(1, 21)]
    mito_genes = [f"MT-CO{i}" for i in range(1, 4)] + [f"MT-ND{i}" for i in range(1, 7)]
    ribo_genes = [f"RPS{i}" for i in range(1, 11)] + [f"RPL{i}" for i in range(1, 11)]
    endogenous_genes = [f"GENE_{i:04d}" for i in range(1, 201)]
    all_genes = ercc_genes + mito_genes + ribo_genes + endogenous_genes
    n_genes = len(all_genes)
    gene_indices = {g: i for i, g in enumerate(all_genes)}

    # High-depth full-length sequencing counts (Negative Binomial distribution)
    counts = rng.negative_binomial(n=8, p=0.015, size=(n_cells, n_genes)).astype(np.float32)
    # Zero out ERCC genes initially to control spike-in volume exactly
    counts[:, : len(ercc_genes)] = 0.0

    # Simulate 8 failed/empty wells (rows H05-H12) with low endogenous RNA & high ERCC spike-in fraction
    failed_indices = [n_cells - 8 + j for j in range(8)]
    for idx in range(n_cells):
        if idx in failed_indices:
            counts[idx, len(ercc_genes) :] *= 0.02
            for eg in ercc_genes:
                counts[idx, gene_indices[eg]] = rng.integers(600, 1500)
        else:
            for eg in ercc_genes:
                counts[idx, gene_indices[eg]] = rng.integers(30, 80)
            # Add biological condition effect to first 20 endogenous genes for Treatment wells
            if idx < 48:
                for bg in endogenous_genes[:20]:
                    counts[idx, gene_indices[bg]] += rng.integers(100, 400)

    obs = pd.DataFrame(
        {
            "plate": ["Plate1"] * n_cells,
            "condition": ["Treated"] * 48 + ["Control"] * 48,
            "donor": ["Donor_1"] * 48 + ["Donor_2"] * 48,
        },
        index=cell_names,
    )
    var = pd.DataFrame(
        {
            "gene_name": all_genes,
            "gene_id": [f"ENSG_{i:08d}" for i in range(n_genes)],
            "is_ercc": [g.startswith("ERCC-") for g in all_genes],
        },
        index=all_genes,
    )

    adata = ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    parse_plate_metadata(adata)

print(f"Initial SMART-seq2 Dataset: {adata.n_obs} wells/cells × {adata.n_vars} features")
print(f"Plates detected: {adata.obs['plate'].unique().tolist()}")

# %% [markdown]
# ## 2. Calculate SMART-seq2 QC Metrics (Endogenous Reads, Mito & ERCC Spike-ins)
# In SMART-seq2, tracking Ambion ERCC spike-in percentages is crucial to identify empty wells or failed reverse transcription.

# %%
adata = calculate_qc_metrics(
    adata,
    mito_prefix=tuple(cfg["qc"]["mito_prefix"]),
    ribo_prefix=tuple(cfg["qc"]["ribo_prefix"]),
    ercc_prefix=tuple(cfg["qc"].get("ercc_prefix", ["ERCC-", "ercc-"])),
)

fig_violins = plot_qc_violins(
    adata,
    keys=("n_genes_by_counts", "total_counts", "pct_counts_mito", "pct_counts_ercc"),
)
plt.show()

# %% [markdown]
# ## 3. Plate Layout Heatmap Inspection (96/384-Well Spatial Patterns)
# Visualize spatial distribution across the 96-well plate to identify edge effects or dispensing failures.

# %%
if "well_row" in adata.obs and "well_col" in adata.obs:
    fig_plate_depth = plot_plate_layout(
        adata, plate_id="Plate1", color_key="total_counts", title="Plate 1: Total Read Counts"
    )
    plt.show()

    fig_plate_ercc = plot_plate_layout(
        adata,
        plate_id="Plate1",
        color_key="pct_counts_ercc",
        title="Plate 1: ERCC Spike-in Fraction (%)",
    )
    plt.show()

# %% [markdown]
# ## 4. Optional Multiplet / Doublet Detection
# For SMART-seq2, cells are sorted 1-cell-per-well via FACS into plates; doublet rates are near zero.
# Doublet detection is disabled by default in config (`run_doublet_detection: false`).

# %%
if cfg["qc"].get("run_doublet_detection", False):
    adata = detect_doublets(
        adata,
        expected_doublet_rate=cfg["qc"].get("expected_doublet_rate", 0.05),
        random_state=cfg["project"]["random_seed"],
    )

# %% [markdown]
# ## 5. Filter Low-Quality Cells and Spike-in Contaminated Wells
# Filter out empty wells with high ERCC spike-in percentage (> 15%) and apoptotic cells (> 15% mitochondrial reads).

# %%
filtered_adata = filter_cells(
    adata,
    min_genes=min(cfg["qc"]["min_genes"], 50),
    max_genes=cfg["qc"].get("max_genes"),
    min_counts=min(cfg["qc"]["min_counts"], 5000),
    max_counts=cfg["qc"].get("max_counts"),
    max_pct_mito=cfg["qc"]["max_pct_mito"],
    max_pct_ercc=cfg["qc"].get("max_pct_ercc", 15.0),
    max_doublet_score=cfg["qc"].get("max_doublet_score"),
)

filtered_adata = filter_genes(
    filtered_adata,
    min_cells=cfg["filter_genes"]["min_cells"],
)

print(
    f"Post-Filtering Dataset: {filtered_adata.n_obs} wells retained ({adata.n_obs - filtered_adata.n_obs} failed wells removed) × {filtered_adata.n_vars} genes"
)

# %% [markdown]
# ## 6. Save Filtered AnnData
# Export clean dataset to `data/processed/` for downstream normalization and clustering.

# %%
out_target = output_h5ad or "data/processed/01_qc_filtered.h5ad"
save_h5ad(filtered_adata, out_target)
print(f"QC & Filtering step completed. Output saved to: {out_target}")
