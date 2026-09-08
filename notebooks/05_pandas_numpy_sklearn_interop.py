# %% [markdown]
# Start with [00_end_to_end.ipynb](00_end_to_end.ipynb) for an executable synthetic walkthrough. This topic notebook requires the indicated input checkpoints and dataset-specific parameters. Open from workbench/notebooks/.

# %% [markdown]
# # 🐼 Integration Guide: Pandas, NumPy, Scikit-Learn & sc-workbench
#
# **Core Philosophy**:
# - Anything **Pandas, NumPy, and Scikit-Learn** can do natively (data wrangling, matrix manipulation, standard ML algorithms, custom feature engineering), you can do using standard tools.
# - Anything they **cannot** do easily (single-cell manifold topology, Scrublet doublet simulation, Harmony batch integration, Leiden community detection on sparse k-NN graphs, PAGA pseudotime, biomarker ranking, Enrichr), **`sc_workbench`** takes care of seamlessly.
# - You can move back and forth between DataFrames, NumPy arrays, and `SingleCellWorkbench` at any stage.

# %%
import sys
from pathlib import Path

sys.path.insert(0, str(Path("..").resolve()))

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import RandomForestClassifier

from sc_workbench import SingleCellWorkbench

print("Libraries loaded successfully.")

# %% [markdown]
# ## 1. Prepare Data with Pandas & NumPy
# Create or manipulate an expression matrix using standard Pandas DataFrames.

# %%
# Create mock single-cell count matrix in Pandas (100 cells x 200 genes)
np.random.seed(42)
cells = [f"cell_{i:03d}" for i in range(100)]
genes = (
    [f"GENE_{j}" for j in range(180)]
    + [f"MT-ATP{k}" for k in range(10)]
    + [f"RPS{m}" for m in range(10)]
)
counts = np.random.poisson(lam=12, size=(100, 200))

raw_df = pd.DataFrame(counts, index=cells, columns=genes)
print(f"Raw Pandas DataFrame shape: {raw_df.shape}")
raw_df.iloc[:5, :5]

# %% [markdown]
# ## 2. Ingest into sc-workbench for Specialized Single-Cell Processing
# Pass the DataFrame directly into `SingleCellWorkbench.from_df()`.

# %%
# Ingest directly from Pandas
wb = SingleCellWorkbench.from_df(raw_df)

# Run single-cell specific algorithms (Mito QC, Doublet detection, Log1p, HVG, PCA, k-NN graph, UMAP)
wb = (
    wb.calculate_qc()
    .filter_cells(min_genes=20, min_counts=50, max_pct_mito=20.0)
    .normalize(target_sum=1e4)
    .select_hvg(n_top_genes=100)
    .scale()
    .run_pca(n_comps=15)
    .compute_neighbors(n_neighbors=10, n_pcs=10)
    .run_umap()
)
wb

# %% [markdown]
# ## 3. Compare Clustering: Scikit-Learn (KMeans) vs. Single-Cell Graph (Leiden)
# You can run standard machine learning clustering (`KMeans`, `Agglomerative`, `DBSCAN`) or single-cell graph clustering (`Leiden`, `Louvain`) on the same dataset.

# %%
# A. Classical Scikit-Learn KMeans
wb.cluster(method="kmeans", n_clusters=4, key_added="kmeans_clusters")

# B. Specialized Single-Cell Graph Modularity Clustering (Leiden)
wb.cluster(method="leiden", resolution=0.6, key_added="leiden_clusters")

# Inspect results directly via the Pandas .obs property
wb.obs[["total_counts", "n_genes_by_counts", "kmeans_clusters", "leiden_clusters"]].head()

# %% [markdown]
# ## 4. Apply Arbitrary Scikit-Learn Estimators & Supervised Learning
# Use `apply_sklearn()` for decomposition/transformers and `train_classifier()` to evaluate cluster separability.

# %%
# Apply TruncatedSVD from Scikit-Learn onto PCA embeddings
svd = TruncatedSVD(n_components=3, random_state=42)
wb.apply_sklearn(svd, input_source="X_pca", key_added="X_svd")
print("Stored custom SVD in .obsm['X_svd']:", wb.adata.obsm["X_svd"].shape)

# Train a Random Forest classifier using Scikit-Learn to assess cluster separability
rf = RandomForestClassifier(n_estimators=50, random_state=42)
clf_res = wb.train_classifier(rf, target_col="leiden_clusters", use_rep="X_pca", cv=3)
print(f"Random Forest 3-Fold Cross-Validation Accuracy: {clf_res['mean_cv_accuracy']:.2%}")

# %% [markdown]
# ## 5. Seamless Export Back to Pandas DataFrames & NumPy Arrays
# Extract anything back to Pandas or NumPy for customized analysis, Seaborn plotting, or custom file saving.

# %%
# 1. Export normalized expression matrix to Pandas DataFrame
expr_df = wb.to_df(layer="normalized")
print("Expression DataFrame:", type(expr_df), expr_df.shape)

# 2. Export 2D UMAP coordinates to Pandas DataFrame
umap_df = wb.get_embedding("X_umap", as_df=True)
print("UMAP DataFrame:\n", umap_df.head(3))

# 3. Export raw matrix to dense NumPy array
dense_mat = wb.to_numpy(layer="counts")
print("Dense NumPy matrix:", type(dense_mat), dense_mat.shape)

# 4. Direct Pythonic Slicing using Pandas boolean mask
mask = wb.obs["total_counts"] > 1000
high_depth_wb = wb[mask]
print(f"Filtered from {wb.n_cells} cells to {high_depth_wb.n_cells} cells with >1000 counts")
