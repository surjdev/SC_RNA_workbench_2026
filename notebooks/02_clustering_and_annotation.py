# %% [markdown]
# # 🧬 Template 2: Normalization, Clustering & Cell Typing (SMART-seq2)
# **Workflow:** Endogenous Filtering ➔ CPM Normalization ➔ HVG Selection ➔ PCA ➔ Plate Batch Check ➔ UMAP ➔ Leiden ➔ Marker Discovery
#
# Part of the Single-Cell SMART-seq2 Transcriptomics Analysis Workbench (FR-10, NFR-5).
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
print(f"Loaded SMART-seq2 AnnData: {adata.n_obs} cells/wells × {adata.n_vars} features")

# %% [markdown]
# ## 2. SMART-seq2 CPM Normalization (Excluding ERCC Spike-ins) & Log Transformation
# For SMART-seq2 full-length sequencing, endogenous gene expression is scaled to Counts Per Million (CPM, target_sum=1,000,000).
# Best practice: synthetic ERCC spike-in transcripts are excluded from the endogenous library size calculation.

# %%
if "is_ercc" in adata.var and adata.var["is_ercc"].any():
    n_ercc = int(adata.var["is_ercc"].sum())
    print(
        f"Filtering out {n_ercc} ERCC spike-in controls from endogenous library size normalization..."
    )
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    adata_endogenous = adata[:, ~adata.var["is_ercc"]].copy()
else:
    adata_endogenous = adata.copy()

target_sum = cfg["normalization"].get("target_sum", 1000000.0)
sc.pp.normalize_total(adata_endogenous, target_sum=target_sum)
if cfg["normalization"].get("log1p", True):
    sc.pp.log1p(adata_endogenous)
adata_endogenous.layers["normalized"] = adata_endogenous.X.copy()
print(
    f"CPM Normalization (target_sum={target_sum:,.0f}) & log1p complete: {adata_endogenous.n_vars} endogenous genes."
)

# %% [markdown]
# ## 3. Highly Variable Gene (HVG) Selection
# Identify genes driving biological heterogeneity.

# %%
n_top = min(cfg["normalization"]["n_top_genes"], adata_endogenous.n_vars)
flavor = cfg["normalization"].get("flavor", "seurat")
layer = "counts" if flavor == "seurat_v3" and "counts" in adata_endogenous.layers else None
sc.pp.highly_variable_genes(
    adata_endogenous,
    n_top_genes=n_top,
    flavor=flavor,
    layer=layer,
)
print(
    f"HVG Selection complete: {adata_endogenous.var['highly_variable'].sum()} variable genes tagged."
)

# %% [markdown]
# ## 4. Dimensionality Reduction (PCA & k-NN Graph)
# Compute linear embedding and cell-cell similarity graph.

# %%
n_pcs = min(cfg["reduction"]["n_pcs"], adata_endogenous.n_obs - 1, adata_endogenous.n_vars - 1)
sc.tl.pca(
    adata_endogenous,
    n_comps=n_pcs,
    use_highly_variable=True,
    random_state=cfg["project"]["random_seed"],
)

sc.pp.neighbors(
    adata_endogenous,
    n_neighbors=min(cfg["reduction"]["n_neighbors"], adata_endogenous.n_obs - 1),
    n_pcs=n_pcs,
    metric=cfg["reduction"].get("metric", "cosine"),
    random_state=cfg["project"]["random_seed"],
)

# %% [markdown]
# ## 5. Manifold Embedding & Community Detection (UMAP & Leiden)
# Project cells onto 2D space and partition into discrete clusters.

# %%
sc.tl.umap(
    adata_endogenous,
    min_dist=cfg["reduction"].get("umap_min_dist", 0.5),
    spread=cfg["reduction"].get("umap_spread", 1.0),
    random_state=cfg["project"]["random_seed"],
)

sc.tl.leiden(
    adata_endogenous,
    resolution=cfg["clustering"]["resolution"],
    random_state=cfg["project"]["random_seed"],
)

print(
    f"Leiden Clustering complete: {adata_endogenous.obs['leiden'].nunique()} clusters identified."
)

# %% [markdown]
# ## 6. Visualizing Clusters, Plates & Biological Conditions
# Generate publication-grade UMAP representations, checking for plate batch effects.

# %%
color_vars = ["leiden"]
if "plate" in adata_endogenous.obs and adata_endogenous.obs["plate"].nunique() > 1:
    color_vars.append("plate")
if "condition" in adata_endogenous.obs:
    color_vars.append("condition")

sc.pl.umap(
    adata_endogenous,
    color=color_vars,
    frameon=False,
    title=[f"Leiden (res={cfg['clustering']['resolution']})"] + color_vars[1:],
    show=False,
)
plt.show()

# %% [markdown]
# ## 7. Cluster Biomarker Ranking
# Identify differentially expressed marker genes per cluster.

# %%
sc.tl.rank_genes_groups(
    adata_endogenous,
    groupby="leiden",
    method=cfg.get("markers", {}).get("method", "wilcoxon"),
)

marker_df = sc.get.rank_genes_groups_df(adata_endogenous, group=None)
print("Top Marker Genes Summary:")
top_markers = marker_df.groupby("group").head(3)
print(top_markers[["group", "names", "scores", "logfoldchanges", "pvals_adj"]])

# %% [markdown]
# ## 8. Save Processed Dataset
# Export clustered AnnData to `data/processed/` for differential expression analysis.

# %%
out_target = output_h5ad or "data/processed/02_clustered.h5ad"
save_h5ad(adata_endogenous, out_target)
print(f"Clustering & Annotation step completed. Output saved to: {out_target}")
