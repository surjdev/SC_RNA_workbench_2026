# %% [markdown]
# # 🧬 Template 2: Normalization, Clustering & Cell Typing
# **Workflow:** Normalization ➔ HVG Selection ➔ PCA ➔ k-NN ➔ UMAP ➔ Leiden ➔ Marker Discovery
#
# Part of the Single-Cell Transcriptomics Analysis Workbench (FR-10, NFR-5).
# Calls native **Scanpy** APIs directly with **workbench_utils** parameterization.

# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from workbench_utils.config import load_config
from workbench_utils.io import load_h5ad, save_h5ad
from workbench_utils.plotting import set_publication_style

set_publication_style()

# %%
# ==============================================================================
# Parameters Cell for Papermill (FR-7)
# ==============================================================================
config_path = "configs/default_analysis.yaml"
input_h5ad = "data/processed/01_qc_filtered.h5ad"
output_h5ad = "data/processed/02_clustered.h5ad"

# %%
cfg = load_config(config_path if Path(config_path).exists() else None)
np.random.seed(cfg["project"]["random_seed"])
print(f"Loaded Pipeline Configuration: {cfg['project']['name']}")

# %% [markdown]
# ## 1. Load Filtered Dataset
# Load AnnData object produced by Template 1.

# %%
input_file = Path(input_h5ad)
if not input_file.exists():
    raise FileNotFoundError(
        f"Input file not found at: {input_file.resolve()}. Please run 01_qc_and_filtering first!"
    )

adata = load_h5ad(input_file)
print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

# %% [markdown]
# ## 2. Normalization & Log Transformation
# Scale counts to target sum and apply natural logarithm.

# %%
sc.pp.normalize_total(adata, target_sum=cfg["normalization"]["target_sum"])
if cfg["normalization"].get("log1p", True):
    sc.pp.log1p(adata)
adata.layers["normalized"] = adata.X.copy()
print("Normalization & log1p complete.")

# %% [markdown]
# ## 3. Highly Variable Gene (HVG) Selection
# Identify genes driving biological heterogeneity.

# %%
n_top = min(cfg["normalization"]["n_top_genes"], adata.n_vars)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=n_top,
    flavor=cfg["normalization"].get("flavor", "seurat"),
    layer="counts" if "counts" in adata.layers else None,
)
print(f"HVG Selection complete: {adata.var['highly_variable'].sum()} variable genes tagged.")

# %% [markdown]
# ## 4. Dimensionality Reduction (PCA & k-NN Graph)
# Compute linear embedding and cell-cell similarity graph.

# %%
n_pcs = min(cfg["reduction"]["n_pcs"], adata.n_obs - 1, adata.n_vars - 1)
sc.tl.pca(
    adata,
    n_comps=n_pcs,
    use_highly_variable=True,
    random_state=cfg["project"]["random_seed"],
)

sc.pp.neighbors(
    adata,
    n_neighbors=cfg["reduction"]["n_neighbors"],
    n_pcs=n_pcs,
    metric=cfg["reduction"].get("metric", "cosine"),
    random_state=cfg["project"]["random_seed"],
)

# %% [markdown]
# ## 5. Manifold Embedding & Community Detection (UMAP & Leiden)
# Project cells onto 2D space and partition into discrete clusters.

# %%
sc.tl.umap(
    adata,
    min_dist=cfg["reduction"].get("umap_min_dist", 0.5),
    spread=cfg["reduction"].get("umap_spread", 1.0),
    random_state=cfg["project"]["random_seed"],
)

sc.tl.leiden(
    adata,
    resolution=cfg["clustering"]["resolution"],
    random_state=cfg["project"]["random_seed"],
)

print(f"Leiden Clustering complete: {adata.obs['leiden'].nunique()} clusters identified.")

# %% [markdown]
# ## 6. Visualizing Clusters & Biological Conditions
# Generate publication-grade UMAP representations.

# %%
color_vars = ["leiden"]
if "condition" in adata.obs:
    color_vars.append("condition")

sc.pl.umap(
    adata,
    color=color_vars,
    frameon=False,
    title=[f"Leiden Clusters (res={cfg['clustering']['resolution']})", "Condition"],
    show=False,
)
plt.show()

# %% [markdown]
# ## 7. Cluster Biomarker Ranking
# Identify differentially expressed marker genes per cluster.

# %%
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method=cfg.get("markers", {}).get("method", "wilcoxon"),
)

marker_df = sc.get.rank_genes_groups_df(adata, group=None)
print("Top Marker Genes Summary:")
top_markers = marker_df.groupby("group").head(3)
print(top_markers[["group", "names", "scores", "logfoldchanges", "pvals_adj"]])

# %% [markdown]
# ## 8. Save Processed Dataset
# Export clustered AnnData to `data/processed/` for differential expression analysis.

# %%
out_target = output_h5ad or "data/processed/02_clustered.h5ad"
save_h5ad(adata, out_target)
print(f"Clustering & Annotation step completed. Output saved to: {out_target}")
