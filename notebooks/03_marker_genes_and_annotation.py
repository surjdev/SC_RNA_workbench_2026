# ---
# jupyter:
#   language_info:
#     name: python
#     version: '3.11'
# ---

# %% [markdown]
# Start with [00_end_to_end.ipynb](00_end_to_end.ipynb) for an executable synthetic walkthrough. This topic notebook requires the indicated input checkpoints and dataset-specific parameters. Open from workbench/notebooks/.

# %% [markdown]
# # 🔬 03. Marker Gene Discovery & Cell Type Annotation
#
# This notebook demonstrates differential expression analysis across clusters, biomarker ranking, Volcano plot visualization, and cell type assignment.

# %%
import sys
from pathlib import Path

sys.path.insert(0, str(Path("..").resolve()))

import sc_workbench as scw
from sc_workbench import SingleCellWorkbench

wb = SingleCellWorkbench.from_h5ad("../data/clustered.h5ad")
wb

# %% [markdown]
# ### 1. Identify Cluster Marker Genes (Wilcoxon Rank-Sum)

# %%
markers_df = wb.find_markers(groupby="leiden", method="wilcoxon", n_genes=25)
markers_df.head(15)

# %% [markdown]
# ### 2. Top Markers Summary per Cluster

# %%
top_markers = markers_df.groupby("cluster").head(3)
display(top_markers[["cluster", "gene", "logfoldchange", "pvals_adj"]])

# %% [markdown]
# ### 3. Cluster Biomarker DotPlot

# %%
top_genes_dict = markers_df.groupby("cluster")["gene"].apply(lambda x: list(x[:3])).to_dict()
scw.plotting.plot_marker_dotplot(wb.adata, markers=top_genes_dict, groupby="leiden")

# %% [markdown]
# ### 4. Biological Cell Type Assignment

# %%
# Map cluster IDs to biological cell type designations
# Fill only after inspecting tissue-specific markers; cluster IDs carry no cell identity.
cluster_map = {str(k): "Unassigned" for k in wb.obs["leiden"].unique()}
wb.annotate_cell_types(cluster_map, cluster_key="leiden", new_key="cell_type")
scw.plotting.plot_umap(wb.adata, color="cell_type")

# %% [markdown]
# ### 5. Save Annotated Object

# %%
wb.save("../data/annotated.h5ad")
print("✔ Stage 03 completed successfully!")
