# ---
# jupyter:
#   language_info:
#     name: python
#     version: '3.11'
# ---

# %% [markdown]
# Start with [00_end_to_end.ipynb](00_end_to_end.ipynb) for an executable synthetic walkthrough. This topic notebook requires the indicated input checkpoints and dataset-specific parameters. Open from workbench/notebooks/.

# %% [markdown]
# # 🎯 02. Normalization, Dimensionality Reduction & Clustering
#
# This notebook demonstrates depth scaling, log1p transformation, Highly Variable Gene (HVG) selection, PCA, kNN graph building, 2D UMAP projection, and Leiden clustering.

# %%
import sys
from pathlib import Path

sys.path.insert(0, str(Path("..").resolve()))

import scanpy as sc

import sc_workbench as scw
from sc_workbench import SingleCellWorkbench

wb = SingleCellWorkbench.from_h5ad("../data/qc_filtered.h5ad")
wb

# %% [markdown]
# ### 1. Depth Normalization & Variance Stabilization

# %%
wb.normalize(target_sum=1e4)
wb.select_hvg(n_top_genes=2000)


# %% [markdown]
# ### 2. PCA & Feature Scaling

# %%
wb.scale(max_value=10.0)
wb.run_pca(n_comps=30)
sc.pl.pca_variance_ratio(wb.adata, n_pcs=20, log=True)

# %% [markdown]
# ### 3. k-NN Neighborhood Graph & 2D UMAP Embedding

# %%
wb.compute_neighbors(n_neighbors=15, n_pcs=20)
wb.run_umap(min_dist=0.3)
scw.plotting.plot_umap(wb.adata, color=["total_counts", "n_genes_by_counts"])

# %% [markdown]
# ### 4. Leiden Community Clustering

# %%
wb.cluster(resolution=0.6, method="leiden")
scw.plotting.plot_umap(wb.adata, color="leiden")

# %% [markdown]
# ### 5. Save Clustered Object

# %%
wb.save("../data/clustered.h5ad")
print("✔ Stage 02 completed successfully!")
